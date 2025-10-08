#include "ForcefieldParser.h"
#include "Constants.h"
#include "StringUtils.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

namespace fs = filesystem;

static void parseRTP(
    const fs::path& filename,
    unordered_map<string, unordered_map<string, double>>& chargeMap)
{
    ifstream in(filename);
    if (!in) throw runtime_error("Could not open " + filename.string());

    string currentResidue;
    bool keep = false;

    string line;
    while (getline(in, line)) {
        if (isOnlyWhitespace(line) || line[0] == ';') continue;

        trimRef(line);
        if (line.front() == '[') {
            string name = trim(line.substr(1, line.size() - 2));

            if (name == "atoms") {
                keep = true;
            }
            else if (name == "bonds" || name == "impropers" || name == "cmap") {
                keep = false;
            }
            else {
                currentResidue = name;
                chargeMap[currentResidue] = {};
            }
            continue;
        }

        if (!currentResidue.empty() && keep) {
            istringstream iss(line);
            string atomName, atomType;
            double charge;
            int group;

            if (iss >> atomName >> atomType >> charge >> group) {
                chargeMap[currentResidue][atomName] = charge;
            }
        }
    }
}

static void parseTDB(
    const fs::path& filename,
    unordered_map<string, unordered_map<string, double>>& chargeMap)
{
    ifstream in(filename);
    if (!in) throw runtime_error("Could not open " + filename.string());

    string currentResidue;
    bool keep = false; 
    bool replace = false;
    bool add = false;

    string line;
    while (getline(in, line)) {
        if (isOnlyWhitespace(line) || line[0] == ';') continue;

        trimRef(line);
        if (line.front() == '[') {
            string name = trim(line.substr(1, line.size() - 2));

            if (name == "replace") {
                keep = true;
                replace = true;
                add = false;
            }
            else if (name == "add") {
                keep = true;
                replace = false;
                add = true;
            }
            else if (name == "delete" || name == "impropers" || name == "None" || name == "add") {
                keep = false;
                replace = false; 
                add = false;
            }
            else {
                currentResidue = name;
                chargeMap[currentResidue] = {};
                keep = false; 
                replace = false;
                add = false;
            }
            continue;
        }

        if (!currentResidue.empty() && keep) {
            bool passed = false;
            string atomName; 
            double charge;
            if (replace) {
                istringstream iss(line);
                string atomType;
                double mass;
                if (iss >> atomName >> atomType >> mass >> charge) 
                    passed = true;
            }
            else if (add) {
                istringstream iss(line);
                double mass;
                int group;
                if (iss >> atomName >> mass >> charge >> group)
                    passed = true;
            }

            if (!passed) continue;

            if (currentResidue.find("NH") != string::npos && (atomName == "HC" || atomName == "H")) {
                for (const string& hName : { "H1", "H2", "H3" })
                    chargeMap[currentResidue][hName] = charge;
            }
            else if (currentResidue == "COO-" && atomName == "OC") {
                for (const string& hName : { "OT1", "OT2", "OXT"})
                    chargeMap[currentResidue][hName] = charge;
            }
            else if (currentResidue == "COOH") {
                if (atomName == "OB")
                    chargeMap[currentResidue]["OT1"] = charge;
                else if (atomName == "OT2") {
                    chargeMap[currentResidue][atomName] = charge;
                    chargeMap[currentResidue]["OXT"] = charge;
                }
                else if (atomName == "H")
                    chargeMap[currentResidue]["HT2"] = charge;
                else 
                    chargeMap[currentResidue][atomName] = charge;
            }
            else {
                chargeMap[currentResidue][atomName] = charge;
            }
        }
    }
}

unordered_map<string, unordered_map<string, double>> getForceFieldCharges() {
    unordered_map<string, unordered_map<string, double>> chargeMap;
    fs::path ff_dir = "charmm36.ff";
    parseRTP(ff_dir / "aminoacids.rtp", chargeMap);
    parseTDB(ff_dir / "aminoacids.n.tdb", chargeMap);
    parseTDB(ff_dir / "aminoacids.c.tdb", chargeMap);
    return chargeMap;
}