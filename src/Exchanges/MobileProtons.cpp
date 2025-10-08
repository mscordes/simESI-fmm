#include "MobileProtons.h"
#include "Charge.h"
#include "Constants.h"
#include "FileUtils.h"
#include "FindExchanges.h"
#include "ForcefieldParser.h"
#include "MathUtils.h"
#include "TopologyUtils.h"
#include "Transforms.h"
#include "RNG.h"

MobileProtonSite::MobileProtonSite(
    const shared_ptr<Residue>& residue_, int idx_,
    bool nterm_, bool cterm_)
    : residue(residue_), idx(idx_), nterm(nterm_), cterm(cterm_)
{    
    protonated = isAA_Protonated(residue, nterm, cterm);

    static const unordered_set<string> basicAminoAcids  = { "LYS", "ARG",  "HIS" };
    static const unordered_set<string> acidicAminoAcids = { "GLU", "ASP" };

    const string& resName = residue->name;
    if (protonated && (nterm || basicAminoAcids.contains(resName)))
        charge = 1.0;
    else if (!protonated && (cterm || acidicAminoAcids.contains(resName)))
        charge = -1.0;
    else
        charge = 0.0; 

    string acceptorName;
    if (nterm)
        acceptorName = *Constants::acceptorsMap.at("NTERM").begin();
    else if (cterm)
        acceptorName = *Constants::acceptorsMap.at("CTERM").begin();
    else {
        auto it = Constants::acceptorsMap.find(resName); 
        if (it == Constants::acceptorsMap.end())
            throw runtime_error("MobileProtonSite: Could not find residue " + resName + " in acceptorsMap.");
        acceptorName = *it->second.begin();
    }

    bool found = false;
    for (const auto& a : residue->atoms) {
        if (acceptorName == a->name) {
            coord = a->coord;
            found = true;
            break;
        }
    }

    if (!found)
        throw runtime_error("MobileProtonSite: Could not find acceptor atoms of residue " + resName);
}

// Assess protein intramolecular proton transfers (mobile protons or salt-bridge annihilation)
void mobileProtons(
    State& state,
    RunFlags& flags,
    ExchangeWriter& exchangeWriter,
    const unordered_map<int, Cluster>& clusters,
    const int step,
    double temp)
{
    // Find all mobile proton donor/acceptors, record sites in pdb2gmx order
    // NOTE: sites must be in the gas-phase for MPM transfer
    unordered_map<char, vector<shared_ptr<MobileProtonSite>>> aminoAcids;
    unordered_set<shared_ptr<MobileProtonSite>> hydrogens, acceptors;

    int resCount = 0;
    for (const auto& [chain, p] : state.proteinMap) {
        for (const auto& resType : Constants::pdb2gmxOrder) {
            for (const auto& r : p->residues) {
                if (r->name == resType) {
                    shared_ptr<MobileProtonSite> site =
                        make_shared<MobileProtonSite>(r, resCount, false, false);

                    aminoAcids[chain].push_back(site);
                    ++resCount;

                    if (clusters.at(r->cluster).water == 0) {
                        if (isAA_Protonated(r)) hydrogens.insert(site);
                        else                    acceptors.insert(site);
                    }
                }
            }
        }

        // N-termini
        const auto& nterm = p->residues.front();
        shared_ptr<MobileProtonSite> nterm_site =
            make_shared<MobileProtonSite>(nterm, resCount, true, false);

        aminoAcids[chain].push_back(nterm_site);
        ++resCount;

        if (clusters.at(nterm->cluster).water == 0) {
            if (isAA_Protonated(nterm, true, false)) hydrogens.insert(nterm_site);
            else                                     acceptors.insert(nterm_site);
        }

        // C-termini
        const auto& cterm = p->residues.back();
        shared_ptr<MobileProtonSite> cterm_site =
            make_shared<MobileProtonSite>(cterm, resCount, false, true);

        aminoAcids[chain].push_back(cterm_site);
        ++resCount;

        if (clusters.at(cterm->cluster).water == 0) {
            if (isAA_Protonated(cterm, false, true)) hydrogens.insert(cterm_site);
            else                                     acceptors.insert(cterm_site);
        }
    }
    if (hydrogens.size() == 0 && acceptors.size() == 0) return;

    // Store a copy of original protonation states
    unordered_map<char, vector<shared_ptr<MobileProtonSite>>> preAminoAcids;

    // Add in charged adducts as well, these will only be considered during ES calc, not for transfer
    vector<shared_ptr<MobileProtonSite>> adducts;
    const unordered_map<string, int> adductCharges = {{"NXH", 1}, {"ATX", -1}};

    size_t totalCount = 0;
    for (const auto& [resType, _] : adductCharges) {
        auto it = state.residueSet.find(resType);
        if (it != state.residueSet.end()) totalCount += it->second.size();
    }
    adducts.reserve(totalCount);

    for (const auto& [resType, charge] : adductCharges) {
        auto it = state.residueSet.find(resType);
        if (it == state.residueSet.end()) continue;

        for (const auto& r : it->second) {
            adducts.emplace_back(make_shared<MobileProtonSite>(r, resCount, false, false));
            ++resCount;
        }
    }

    // Flatten all MobileProtonSite objects
    vector<shared_ptr<MobileProtonSite>> chargedSites;
    chargedSites.reserve(aminoAcids.size() + adducts.size());

    for (const auto& [chain, sites] : aminoAcids) 
        for (const auto& s : sites) 
            chargedSites.push_back(s);

    for (const auto& s : adducts)
        chargedSites.push_back(s);

    // Get coords/charges of all MobileProtonSites
    vector<Vec3D> coords;
    coords.reserve(chargedSites.size());

    vector<double> charges;
    charges.reserve(chargedSites.size());

    for (const auto& s : chargedSites) {
        coords.push_back(s->coord);
        charges.push_back(s->charge);
    }

    // Get all distance between unique pairs, only done once
    const auto distances = cdist(coords);

    // Modify charges/protonation state to "(de)protonate" a residue
    auto modifyResidue = [](
        const shared_ptr<MobileProtonSite>& site, 
        bool protonate, vector<double>& charges) {
            if (protonate) {
                charges[site->idx] += 1.0;
                site->protonated = true; 
            }
            else {
                charges[site->idx] -= 1.0;
                site->protonated = false;
            }
        };


    // Compute change the apparaent GPB of a pattern of protonation
    auto computeGPBapp = [](
        const auto& aminoAcids,
        const vector<vector<double>>& distances, 
        const vector<double>& charges) {
            static const auto& GPBs = Constants::GPBs;
            double GPBapp = 0.0;
            for (const auto& pair : aminoAcids) {
                for (const auto& s : pair.second) {
                    if (s->protonated) {
                        if (s->nterm)      GPBapp += GPBs.at("NTERM");
                        else if (s->cterm) GPBapp += GPBs.at("CTERM");
                        else               GPBapp += GPBs.at(s->residue->name);
                    }
                }
            }
            return GPBapp + computeCoulomb(distances, charges);
        };

    // ----- Find exchanges -----
    unordered_set<shared_ptr<Exchange>> exchanges;
    
    double GPBapp = computeGPBapp(aminoAcids, distances, charges);
    vector<double> GPBapps = { GPBapp };

    vector<shared_ptr<MobileProtonSite>> shuffledHydrogens(hydrogens.begin(), hydrogens.end());
    vector<shared_ptr<MobileProtonSite>> shuffledAcceptors(acceptors.begin(), acceptors.end());

    while (true) {
        shuffle(shuffledHydrogens.begin(), shuffledHydrogens.end(), RNG::gen);
        shuffle(shuffledAcceptors.begin(), shuffledAcceptors.end(), RNG::gen);

        for (const auto& hsite : shuffledHydrogens) {
            modifyResidue(hsite, false, charges);

            for (const auto& asite : shuffledAcceptors) {
                modifyResidue(asite, true, charges);

                double final_GPBapp = computeGPBapp(aminoAcids, distances, charges);
                double delta = final_GPBapp - GPBapp;

                if (MCMC(delta, temp)) {
                    hydrogens.insert(asite);
                    acceptors.insert(hsite);

                    hydrogens.erase(hsite);
                    acceptors.erase(asite);

                    GPBapp = final_GPBapp;

                    // Just use first atom
                    shared_ptr<Atom> h_atom = hsite->residue->atoms[0];
                    shared_ptr<Atom> a_atom = asite->residue->atoms[0];

                    exchanges.emplace(
                        make_shared<Exchange>(state, h_atom, a_atom, delta, step, temp));

                    goto endloop;
                }

                modifyResidue(asite, false, charges);
            }

            modifyResidue(hsite, true, charges);
        }
        endloop:


        // Get that last 25% of iterations
        GPBapps.push_back(GPBapp);
        size_t n = GPBapps.size();
        auto quartile = vector<double>(GPBapps.begin() + (n * 3 / 4), GPBapps.end());

        // Convergence checks
        const double KT = (Constants::R / 1000) * temp;
        if (stdev(quartile) >= 6.0 * KT) continue;
        if (GPBapp - quartile.front() <= 0.1 * KT) continue;
        break;
    }

    // ----- Facilitate exchanges -----
    if (!exchanges.empty()) {
        flags.createRunFile = true; 
        exchangeWriter.write(exchanges);

        // Extract pdb2gmx input from a MobileProtonSite object
        auto getPdb2GmxInput = [](const shared_ptr<MobileProtonSite>& site) -> string {
            const auto& r = site->residue;
            const bool protonated = site->protonated;

            string input;
            if (site->nterm) {
                const auto& codes = getNtermPdb2GmxCodes(r);
                input = protonated ? codes.first : codes.second;
            }
            else if (site->cterm) {
                const auto& codes = Constants::pdb2gmxCodes.at("CTERM");
                input = protonated ? codes.first : codes.second;
            }
            else if (r->name == "HIS") {
                input = protonated ? "2" : "0";
            }
            else {
                const auto& codes = Constants::pdb2gmxCodes.at(r->name);
                input = protonated ? codes.first : codes.second;
            }

            return input;
        };

        // Create .top/.gro with final protonation states
        map<char, vector<string>> pdb2gmxInputs;
        for (const auto& [chain, sites] : aminoAcids)
            for (const auto& s : sites) 
                pdb2gmxInputs[chain].push_back(getPdb2GmxInput(s));

        call_pdb2gmx(to_string(step) + ".top", pdb2gmxInputs, state, true); 
        State postState;
        postState.parseGRO(to_string(step) + ".gro");
        deleteFile(to_string(step) + ".gro");


        // Begin updating coordinates of State given differing protonation states
        auto updateAtoms = [&](
            const shared_ptr<Residue>& pre_r,
            const shared_ptr<Residue>& post_r,
            bool nterm, bool cterm)
            {
                if (isAA_Protonated(pre_r, nterm, cterm) != isAA_Protonated(post_r, nterm, cterm)) {
                    unordered_map<string, pair<Vec3D, Vec3D>> coordVels;
                    for (const auto& a : pre_r->atoms)
                        coordVels[a->name] = { a->coord, a->velocity };

                    const auto& titHydrogens = Constants::hydrogensMap.at(pre_r->name);
                    for (const auto& a : post_r->atoms) {
                        if (titHydrogens.contains(a->name)) {
                            a->velocity = sampleMaxwell(1, temp, a->element)[0];

                            // pdb2gmx screws up Ctermini new H placement so fix that one edge case here
                            if (cterm && a->name == "HT2") {
                                shared_ptr<Atom> OT2, C;
                                for (const auto& a : post_r->atoms) {
                                    if (a->name == "OT2") OT2 = a;
                                    else if (a->name == "C") C = a;
                                }
                                a->coord = rotateBond(C->coord, OT2->coord, a->coord, 107.0);
                                a->coord = setBondLength(OT2->coord, a->coord, 0.1);
                            }
                        }
                        else {
                            const auto& cv = coordVels.at(a->name);
                            a->coord = cv.first;
                            a->velocity = cv.second;
                        }
                    }

                    pre_r->atoms = post_r->atoms;
                }
        };

        for (const auto& [chain, protein] : state.proteinMap) {
            const auto& pre_residues = protein->residues;
            const auto& post_residues = postState.proteinMap.at(chain)->residues;

            for (size_t i = 1; i < pre_residues.size()-1; ++i) {
                const auto& pre_r = pre_residues[i]; 
                const auto& post_r = post_residues[i];
                if (!Constants::titAminoAcids.contains(pre_r->name)) continue;
                updateAtoms(pre_r, post_r, false, false);
            }

            // Termini
            const auto&  pre_nterm = pre_residues.front();
            const auto& post_nterm = post_residues.front();
            updateAtoms(pre_nterm, post_nterm, true, false);

            const auto& pre_cterm = pre_residues.back();
            const auto& post_cterm = post_residues.back();
            updateAtoms(pre_cterm, post_cterm, false, true);
        }
    }
}   