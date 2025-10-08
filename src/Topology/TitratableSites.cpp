#include "TitratableSites.h"
#include "Constants.h"
#include "TopologyUtils.h"

#include <map>

// All non-water titratable sites
void getTitSites(
	const State& state, 
	unordered_set<shared_ptr<Atom>>& hydrogens,
	unordered_set<shared_ptr<Atom>>& acceptors)
{
	hydrogens.clear();
	acceptors.clear();

	// Size estimate
	int aminoAcidCount = 0;
	for (const auto& [_, p] : state.proteinMap)
		for (const auto& r : p->residues)
			aminoAcidCount += 1;

	hydrogens.reserve(aminoAcidCount);
	acceptors.reserve(aminoAcidCount);

	for (const auto& [resType, residues] : state.residueSet) {

		if (Constants::titAminoAcids.contains(resType)) {
			const auto& hydrogenNames = Constants::hydrogensMap.at(resType);
			const auto& acceptorNames = Constants::acceptorsMap.at(resType);

			for (const auto& residue : residues) {
				if (isAA_Protonated(residue)) {
					for (const auto& atom : residue->atoms)
						if (hydrogenNames.contains(atom->name))
							hydrogens.insert(atom);
				}
				else {
					if (resType != "HIS") {
						for (const auto& atom : residue->atoms)
							if (acceptorNames.contains(atom->name))
								acceptors.insert(atom);
					}
					else {
						if (isHISD(residue)) { // HISD
							for (const auto& atom : residue->atoms)
								if (atom->name == "NE2")
									acceptors.insert(atom);
						}
						else { // HISE
							for (const auto& atom : residue->atoms)
								if (atom->name == "ND1")
									acceptors.insert(atom);
						}
					}
				}
			}
		}
		else if (Constants::hydrogensMap.contains(resType)) {
			const auto& hydrogenNames = Constants::hydrogensMap.at(resType);
			for (const auto& residue : residues) {
				for (const auto& atom : residue->atoms)
					if (hydrogenNames.contains(atom->name))
						hydrogens.insert(atom);
			}
		}
		else if (Constants::acceptorsMap.contains(resType)) {
			const auto& acceptorNames = Constants::acceptorsMap.at(resType);
			for (const auto& residue : residues) {
				for (const auto& atom : residue->atoms)
					if (acceptorNames.contains(atom->name))
						acceptors.insert(atom);
			}
		}
	}

	for (const auto& [chain, protein] : state.proteinMap) { // Termini
		const auto& nterm = protein->residues.front();
		if (isAA_Protonated(nterm, true, false)) {
			const auto& hydrogenNames = Constants::hydrogensMap.at("NTERM");
			for (const auto& atom : nterm->atoms)
				if (hydrogenNames.contains(atom->name))
					hydrogens.insert(atom);
		}
		else {
			const auto& acceptorNames = Constants::acceptorsMap.at("NTERM");
			for (const auto& atom : nterm->atoms)
				if (acceptorNames.contains(atom->name))
					acceptors.insert(atom);
		}

		const auto& cterm = protein->residues.back();
		if (isAA_Protonated(cterm, false, true)) {
			const auto& hydrogenNames = Constants::hydrogensMap.at("CTERM");
			for (const auto& atom : cterm->atoms)
				if (hydrogenNames.contains(atom->name))
					hydrogens.insert(atom);
		}
		else {
			const auto& acceptorNames = Constants::acceptorsMap.at("CTERM");
			for (const auto& atom : cterm->atoms)
				if (acceptorNames.contains(atom->name))
					acceptors.insert(atom);
		}
	}
}