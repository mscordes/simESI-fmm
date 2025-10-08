#include "State.h"
#include "Atom.h"
#include "Constants.h"
#include "StringUtils.h"

#include <charconv>
#include <string>
#include <iomanip>
#include <fstream>
#include <stdexcept>

// Only updates coordinates and velocities
void State::update(const string& groFile) {
    ifstream file(groFile);
    if (!file.is_open())
        throw runtime_error("Error opening .gro file: " + groFile);

    string line;
    getline(file, line);
    getline(file, line);

    auto parseDouble = [](string_view sv) -> double {
        while (!sv.empty() && isspace(static_cast<unsigned char>(sv.front())))
            sv.remove_prefix(1);
        while (!sv.empty() && isspace(static_cast<unsigned char>(sv.back())))
            sv.remove_suffix(1);

        double val{};
        auto res = from_chars(sv.data(), sv.data() + sv.size(), val);
        if (res.ec != errc{}) {
            throw runtime_error("Failed to parse number from: '" + string(sv) + "'");
        }
        return val;
        };

    auto updateAtom = [&](const shared_ptr<Atom>& a, const string& line) {
        string_view view(line);

        // Column ranges (0-based, [start, start+len))
        string_view atomName = view.substr(9, 6);
        if (a->name != trim(string(atomName)))
            throw runtime_error("State::updateAtom: Atom mismatch for atom " + a->name +
                " at line " + line);

        a->coord = {
            parseDouble(view.substr(20, 8)),
            parseDouble(view.substr(28, 8)),
            parseDouble(view.substr(36, 8))
        };
        a->velocity = {
            parseDouble(view.substr(44, 8)),
            parseDouble(view.substr(52, 8)),
            parseDouble(view.substr(60, 8))
        };
    };

    for (const auto& [chain, p] : proteinMap) {
        for (const auto& r : p->residues) {
            for (const auto& a : r->atoms) {
                getline(file, line);
                updateAtom(a, line);
            }
        }
    }

    for (const auto& resType : Constants::topOrder) {
        auto it = residueSet.find(resType);
        if (it == residueSet.end()) continue;

        for (const auto& r : it->second) {
            for (const auto& a : r->atoms) {
                getline(file, line);
                updateAtom(a, line);
            }
        }
    }

    file.close();
}