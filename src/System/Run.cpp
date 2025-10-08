#include "System.h"
#include "CoordinateManip.h"
#include "Clustering.h"
#include "DoExchanges.h"
#include "Exchange.h"
#include "Evaporation.h"
#include "FileUtils.h"
#include "FindExchanges.h"
#include "GmxUtils.h"
#include "MobileProtons.h"
#include "ParsePKA.h"
#include "Subprocess.h"
#include "Timer.h"
#include "TitratableSites.h"
#include "TopologyUtils.h"
#include "Recentering.h"
#include "Writers.h"

void System::run() {
	auto tInit = now();
	params = Parameters::Get();

    // ----- Initialize -----
    bool continuation = params.isCont();
    if (!continuation) {
        initialize();
    }
    else {
        string fname = to_string(params.getStepCont() - 1);
        state.parseGRO(to_string(params.getStepCont() - 1) + ".gro");
        call_pdb2gmx(fname + ".top", getPdb2GmxInputs(state), state);
    }

    parsePKA(state);
    Writers writers(fs::path("..") / "data", state);

    int numRestarts = 1; // For checkpoint naming
    int numContSteps = 2000;

    // ----- Main loop -----
    int numSteps = static_cast<int>(round(params.getTime() * 250)); // 1 ns = 250 4ps steps
    int initStep = (continuation) ? params.getStepCont() : 1;
    RunFlags flags{ true, false, false }; // Flags for topology generation

    for (int step = initStep; step < numSteps; ++step) {
        doStep(step, writers, flags, tInit, numRestarts, numContSteps);
        flags.reset();
    }

    auto runtime = duration_cast<chrono::duration<double>>(now() - tInit).count() / 3600.0;
	cout << "\n\nProgram completed in " << setprecision(3) << runtime << " hr.\n";
}

void System::doStep(
    int step, 
    Writers& writers, 
    RunFlags& flags,
    const chrono::steady_clock::time_point& tInit, 
    int& numRestarts, 
    int& numContSteps) 
{
    auto tStart = now();
    static int waterCutoff = (params.getSettingsPresent().test(WATERCUTOFF)) 
        ? params.getWaterCutoff() : getProteinMass(state) * 0.0035;
    static double itemp = params.getITemp();
    static double ftemp = params.getFTemp();
    static bool mpm = params.isMPM(); 

    // ----- Water clustering -----
    WaterInfo waterInfo = getWaterInfo(state);
    auto clusters = clusterWaters(state, waterInfo);

    // ----- Temperature ramp/mobile protons-----
    double temp = itemp; 
    auto it = state.residueSet.find("SOL");
    if (it == state.residueSet.end() || it->second.size() < waterCutoff) {
        temp = ftemp;

        // Mobile protons must be done before titratable site assigning
        if (mpm && step % 10 == 0) { // Only do every 40 ps
            mobileProtons(state, flags, writers.exchangeWriter, clusters, step, temp); 
        }
    }

    // ----- Titratable site info -----
    unordered_set<shared_ptr<Atom>> hydrogens, acceptors;
    getTitSites(state, hydrogens, acceptors);
    clusterTitSites(hydrogens, acceptors, state, waterInfo);
    computeClusterPHs(clusters, state);

    if (step % 50 == 0) flags.createRunFile = true; // Regenerate .tpr at least every 50 steps
    if (step % 250 == 0) parsePKA(state);           // Recalculate residue pKa's every 1 ns
    fixDisulfides(state, flags);                    // Fix missed disulfide bug
    unordered_set<shared_ptr<Residue>> skipList;    // For prevention of 'rattling' of proton between exchange pair

    for (int hop = 0; hop < 5; ++hop) { // Grotthuss hops

        // ----- Exchanges -----
        if (hop != 0) {
            waterInfo = getWaterInfo(state);
            getTitSites(state, hydrogens, acceptors);
        }

        unordered_set<shared_ptr<Exchange>> exchanges = 
            findExchanges(state, hydrogens, acceptors, waterInfo, clusters, temp, step, hop);
        pruneExchanges(exchanges, skipList, flags, temp);

        writers.exchangeWriter.write(exchanges);
        doExchanges(exchanges, state, clusters, hydrogens, acceptors);

        if (exchanges.empty()) break;
    }

    // ----- Final system processing -----
    auto gWaters = removeEvaporated(state, flags);
    recenter(state, gWaters, flags);

    // Finalize .gro
    string prevStepStr = to_string(step - 1);
    string stepStr = to_string(step);
    string preGroName = stepStr + "_pre.gro";
    string newGroName = stepStr + ".gro";
    if (flags.createRunFile)
        state.writeGRO(preGroName, params.getBoxSize());
    else
        copyFile(to_string(step - 1) + ".gro", preGroName);
    auto tExchange = now();

    // Finalize .top
    string topName = stepStr + ".top";
    if (flags.protExchanges)  // New protein .itp's 
        call_pdb2gmx(stepStr + ".top", getPdb2GmxInputs(state), state);
    else if (flags.exchanges) // Just modify res counts in .top
        createTOP(stepStr + ".top", state);
    else
        copyFile(prevStepStr + ".top", stepStr + ".top");
    auto tTop = now();

    static const auto& gmxEnv = params.getGMXEnv();
    modifyMDPgrps("prodrun", state, temp);

    chrono::steady_clock::time_point tGrompp;
    chrono::steady_clock::time_point tMD;

    // ---- Molecular Dynamics -----
    if (flags.createRunFile) {
        { // grompp
            string command = gmxEnv + " grompp -f prodrun.mdp -c " + preGroName + " -p " +
                topName + " -n prodrun.ndx -o " + stepStr + ".tpr -maxwarn 100";
            subprocess(command);
            tGrompp = now();
        }

        { // mdrun
            deleteFile(newGroName);
            call_mdrun(stepStr);
            tMD = now();
        }

        // Reset checkpoint vars
        numRestarts = 1;
        numContSteps = 2000;
    }
    else { // Use Checkpoint
        ++numRestarts;
        numContSteps += 2000;

        { // convert-tpr
            copyFile(prevStepStr + ".cpt", stepStr + ".cpt");
            string command = gmxEnv + " convert-tpr -s " + prevStepStr + ".tpr -o " + stepStr + 
                ".tpr -nsteps " + to_string(numContSteps);
            subprocess(command);
            tGrompp = now();
        }

        { // mdrun 
            static const bool hpc = params.isHPC();
            static const bool gpu = params.isGPU();

            string command = gmxEnv + " mdrun -cpi " + stepStr + " -deffnm "
                + stepStr + " -noappend -cpo " + stepStr;
            if (hpc)        
                command += " -ntmpi 1 -nb gpu";
            else if (gpu)   
                command += " -nb gpu";

            subprocess(command);
            tMD = now();
        }

        for (const auto& ftype : { ".gro", ".xtc", ".edr", ".log" })
            renameCheckpointFile(ftype, step, numRestarts);
    }

    if (!fs::exists(newGroName))
        throw runtime_error("MD failed to complete successfully during step " + stepStr +
            ". \nCheck gmx terminal output.");
    else
        state.update(newGroName);

    // ---- MD file clean-up -----
    static const bool save = params.isSave();
    if (!save) {
        deleteFile(stepStr + ".edr");
        deleteFile(preGroName);
        if (step % 10 != 0) {
            deleteFile(stepStr + ".xtc");
            deleteFile(stepStr + ".log");
        }
        if (step > 2) {
            string oldStep = to_string(step - 2);
            for (const auto& ftype : { ".cpt", ".tpr", ".top" })
                deleteFile(oldStep + ftype);
            deleteFile(oldStep + "_prev.cpt");
            deleteFile(oldStep + "_preMD.ndx");
        }
    }

    // ----- Write data -----
    writers.chargeWriter.write(state);
    writers.compositionWriter.write(state, gWaters);
    writers.temperatureWriter.write(temp);
    writers.timerWriter.write(tInit, tStart, tExchange, tTop, tGrompp, tMD);

    cout << "\rStep " << step << " completed in " << setprecision(3) << seconds(now() - tStart) << " s." << flush;
}