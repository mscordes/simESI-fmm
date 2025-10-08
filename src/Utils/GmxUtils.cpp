#include "GmxUtils.h"
#include "Atom.h"
#include "Constants.h"
#include "CoordinateManip.h"
#include "FileUtils.h"
#include "Parameters.h"
#include "Subprocess.h"

#include <fstream>
#include <map>

// Call mdrun, modify for GPU acceleration or HPC mode dependent on user args
void call_mdrun(const string& fname) {
    static bool hpc = Parameters::Get().isHPC();
    static bool gpu = Parameters::Get().isGPU();
    static string gmxEnv = Parameters::Get().getGMXEnv();

    string command;
    if (hpc)      command = gmxEnv + " mdrun -ntmpi 1 -nb gpu -v -deffnm " + fname;
    else if (gpu) command = gmxEnv + " mdrun -nb gpu -v -deffnm " + fname;
    else          command = gmxEnv + " mdrun -v -deffnm " + fname;
    subprocess(command);
}

// Create .ndx and modify .mdp file for COM translational/rotational removal, T-coupling, FMM, and non-bonded cutoffs
void modifyMDPgrps(const string& opfname, const State& state, double temperature) {
    static bool isHumid = Parameters::Get().isHumid();
    static bool isFMM = Parameters::Get().isFMM();
    static double boxD = Parameters::Get().getBoxSize();
    static string gasTemperature = to_string(Parameters::Get().getGasTemp());
    string temperatureStr = to_string(temperature);

    // ---- ndx groups -----
    map<string, vector<int>> ndxGrps;
    int atomCount = 0;

    for (const auto& [chain, protein] : state.proteinMap) 
        for (const auto& residue : protein->residues) 
            for (int i = 0; i < residue->atoms.size(); ++i) 
                ndxGrps["Protein"].push_back(++atomCount);

    ndxGrps["Protein_Droplet"] = ndxGrps["Protein"];

    for (const auto& resType : Constants::topOrder) {
        if (!state.residueSet.contains(resType)) continue;

        if (resType == "NNN" || resType == "OOO") {
            for (const auto& residue : state.residueSet.at(resType)) 
                for (int i = 0; i < residue->atoms.size(); ++i) 
                    ndxGrps["Atmosphere"].push_back(++atomCount);
        }
        else if (resType == "SOL" && isHumid) { // Differentiate gaseous vs. droplet water
            const auto& skeleton = getProteinSkeletonCoords(state);
            const double cutoff_sq = 5.0 * 5.0; // 5 nm cutoff for now

            for (const auto& r : state.residueSet.at(resType)) {
                const auto& rc = r->atoms[0]->coord;
                
                bool inDroplet = false;
                for (const auto& c : skeleton) {
                    if (rc.distance_sq(c) < cutoff_sq) { 
                        for (int i = 0; i < r->atoms.size(); ++i) {
                            ndxGrps["Droplet"].push_back(++atomCount);
                            ndxGrps["Protein_Droplet"].push_back(atomCount);
                        }
                        inDroplet = true;
                        break;
                    }
                }

                if (!inDroplet) {
                    for (int i = 0; i < r->atoms.size(); ++i) {
                        ndxGrps["Atmosphere"].push_back(++atomCount);
                    }
                }
            }
        }
        else {
            for (const auto& residue : state.residueSet.at(resType)) {
                for (int i = 0; i < residue->atoms.size(); ++i) {
                    ndxGrps["Droplet"].push_back(++atomCount);
                    ndxGrps["Protein_Droplet"].push_back(atomCount);
                }
            }
        }
    }

    // ----- write ndx -----
    string ndxName = opfname + ".ndx";
    deleteFile(ndxName);

    ofstream ndx(ndxName);
    if (!ndx.is_open())
        throw runtime_error("Could not open file " + opfname + ".ndx for writing.");

    for (const auto& [grpName, indices] : ndxGrps) {
        ndx << "[ " << grpName << " ]\n";
        int count = 1;
        for (const auto& idx : indices) {
            ndx << right << setw(4) << idx << " ";
            ++count;
            if (count % 16 == 0) {
                ndx << '\n';
                count = 1;
            }
        }
        ndx << '\n';
    }
    ndx.close();

    // ----- modify mdp -----
    string mdpFname = opfname + ".mdp";
    ifstream imdp(mdpFname);
    if (!imdp.is_open())
        throw runtime_error("Could not open " + mdpFname + " for writing.");

    string line;
    vector<string> lines;
    while (getline(imdp, line)) 
        lines.push_back(line);
    imdp.close();

    bool droplet = ndxGrps.contains("Droplet");
    bool atmosphere = ndxGrps.contains("Atmosphere");

    // tc_grps
    string tcGrps = "Protein ";
    string refT = temperatureStr + " ";
    string tauT = "100 ";
    if (droplet) {
        tcGrps += "Droplet ";
        refT += temperatureStr + " ";
        tauT += "5 ";
    }
    if (atmosphere) {
        tcGrps += "Atmosphere ";
        refT += gasTemperature + " ";
        tauT += "5 ";
    }

    static string rdist = to_string((boxD * 0.99) / 3); // Modify non-bonded cutoffs to fit within box

    for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].find("comm_grps") != string::npos) {
            string commGrps = "Protein";
            if (droplet) commGrps += "_Droplet";
            if (atmosphere) commGrps += " Atmosphere ";
            lines[i] = "comm_grps          = " + commGrps;
        }
        else if (lines[i].find("tc_grps") != string::npos) {
            lines[i] = "tc_grps            = " + tcGrps;
        }
        else if (lines[i].find("ref_t") != string::npos) {
            lines[i] = "ref_t              = " + refT;
        }
        else if (lines[i].find("tau_t") != string::npos) {
            lines[i] = "tau_t              = " + tauT;
        }
        else if (lines[i].find("fmm-override-tree-depth") != string::npos && isFMM) {
            string fmmDepth;
            if (atomCount > 1000000)    fmmDepth = "4";
            else if (atomCount > 50000) fmmDepth = "3";
            else if (atomCount > 30000) fmmDepth = "2";
            else if (atomCount > 10000) fmmDepth = "1";
            else                        fmmDepth = "0";
            lines[i] = "fmm-override-tree-depth = " + fmmDepth;
        }
        else if (lines[i].find("rlist") != string::npos) {
            lines[i] = "rlist              = " + rdist;
        }
        else if (lines[i].find("rcoulomb") != string::npos && lines[i].find("rcoulomb-switch") == string::npos) {
            lines[i] = "rcoulomb           = " + rdist;
        }
        else if (lines[i].find("rvdw") != string::npos && lines[i].find("rvdw-switch") == string::npos) {
            lines[i] = "rvdw               = " + rdist;
        }
    }

    ofstream omdp(mdpFname);
    if (!omdp) throw runtime_error("Could not open " + mdpFname + " for writing.");
    for (const auto& l : lines) {
        omdp << l << '\n';
    }
}