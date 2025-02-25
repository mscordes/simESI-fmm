#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <array>
#include <set>
#include <filesystem>
#include <memory>
#include <Classes.cpp>

namespace Core {

    void printHelp();

    void validateConfig(const Config& config);

    Config getArgs(int argc, char* argv[]);

    std::vector<std::string> readFile(const std::string& fileName);

    std::vector<std::string> splitLine(const std::string& line);

    std::vector<std::vector<std::string>> splitLines(const std::vector<std::string>& lines);

    std::filesystem::path makeDirs(Config& config);

    std::unordered_map<int, float> getPkaVals(const CoordInfo& coordInfo, const std::string& pdb);

    void updatePkavals(std::unordered_map<int, float>& pkaMap, CoordInfo& coordInfo);

    std::vector<std::shared_ptr<Atom>> readGROatoms(const std::string& gro);

    std::vector<std::shared_ptr<Residue>> binAtoms(std::vector<std::shared_ptr<Atom>>& atoms);

    std::vector<std::shared_ptr<Protein>> binResidues(const std::vector<std::shared_ptr<Residue>>& residues);

    std::vector<std::array<float, 3>> extractCoordinates(const std::vector<std::shared_ptr<Atom>>& atoms);

    std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>
        getResidueMap(const std::vector<std::shared_ptr<Residue>>& residues);

    std::unordered_map<std::string, int> getNumResidues(
        const std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>& residueMap);

    std::map<std::string, std::vector<std::shared_ptr<Atom>>> getProteinAtoms(const std::vector<std::shared_ptr<Protein>>& proteins);

    std::vector<std::array<float, 3>> getProteinCoords(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms);

    CoordInfo parseAtoms(std::vector<std::shared_ptr<Atom>> atoms, const float& boxSize);

    CoordInfo buildCoordInfo(const std::string& file, const float& boxSize);

    void writePDB(const std::string& filename, std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms,
        const std::array<float, 3>& box_vectors);

    void auto_gmx_input(const std::string& command, const std::vector<std::string>& inputs);

    void runMD(const std::string& fname, const Config& config);

    void equilibriation(const Config& config, CoordInfo& coordInfo, const std::unordered_map<int, float>& pkaMap, const float& init_boxSize, 
        const std::vector<std::string>& topOrder);

    std::unordered_map<std::string, std::vector<std::string>> setProtStates(const CoordInfo& coordInfo, const std::unordered_map<int, float>& pkaVals,
        const Config& config);

    std::unordered_map<std::string, std::vector<std::string>> getProtStates(const CoordInfo& coordInfo);

    void writeGRO(const std::string& filename, const std::vector<Atom>& atoms, const std::array<float, 3>& box_vectors);

    void writeGRO(const std::string& filename, const std::vector<std::shared_ptr<Atom>>& atoms, const std::array<float, 3>& box_vectors);

    void writeTOP(const std::string& topName, const CoordInfo& coordInfo, const std::unordered_map<std::string,
        std::vector<std::string>>&pdb2gmxMap, const bool keepGro, const Config& config, const std::vector<std::string>& topOrder);

    void createTOP(const std::string& fname, const CoordInfo& coordInfo, const std::vector<std::string>& topOrder);

    void reorderAtoms(CoordInfo& coordInfo, const std::vector<std::string>& topOrder);

    float getProteinMass(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms);

    void copyFile(const std::filesystem::path& ifname, const std::filesystem::path& ofname);

    void deleteFile(const std::filesystem::path& fname);

    void copyFiles(const std::filesystem::path& idir, const std::filesystem::path& odir);

    void copyRecursive(const std::filesystem::path& src, const std::filesystem::path& target) noexcept;

    std::chrono::time_point<std::chrono::high_resolution_clock> timer();

    std::string duration(const std::chrono::time_point<std::chrono::high_resolution_clock>& start,
        const std::chrono::time_point<std::chrono::high_resolution_clock>& end);

    void renameCheckpointFile(const std::string& ftype, const int& step, const int& num_restarts);

    std::vector<std::shared_ptr<Atom>> seedAtmosphere(const Config& config, std::vector<std::shared_ptr<Atom>>& atoms,
        std::vector<std::array<float, 3>>& coords, const std::array<float, 3>& boxVectors, const float& dropletRadius);

    CoordInfo formDroplet(const Config& config, CoordInfo& coordInfo, const std::vector<std::string>& topOrder);

    void modifyMDPgrps(const std::string& ndx_fname, const std::string& mdp_fname, const Config& config,
        const std::vector<std::string>& topOrder, CoordInfo& coordInfo, const float& temperature);

    void topFromCoord(const std::string& groFname, const std::string& topFname, const float& boxSize, const Config& config, 
        const std::vector<std::string>& topOrder);

    void write_output(const std::filesystem::path& target_dir, const std::string& input, const std::string& output_fname);

    std::vector<std::array<float, 3>> sampleMaxwell(int n, float temperature, std::string element);

    std::vector<float> getProtCharges(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms);

    std::vector<float> getCharges(const std::map<std::string, std::vector<std::shared_ptr<Atom>>>& proteinAtoms,
        const std::unordered_map<std::string, int>& numResidues, const std::vector<std::string>& topOrder);
    
    int getNetCharge(const std::vector<float>& charges);

    std::array<float, 3> getCOM(const std::vector<std::array<float, 3>>& coords, const int& numCoords);

    float getDistance(const std::array<float, 3>& coord1, const std::array<float, 3>& coord2); \

    std::array<float, 3> setBondLength(const std::array<float, 3>& coord1,
            const std::array<float, 3>& coord2, const float& length);

    std::vector<std::shared_ptr<Atom>> fixDisulfides(const std::unordered_map<std::string, std::vector<std::shared_ptr<Residue>>>& residueMap,
        std::vector<std::shared_ptr<Atom>>& atoms);

    std::vector<float> norm(const std::vector<std::array<float, 3>>& coords,
        const std::array<float, 3>& refCoord);

    std::vector<std::vector<float>> cdist(const std::vector<std::array<float, 3>>& coords1,
        const std::vector<std::array<float, 3>>& coords2);

    bool removeEvaporated(const Config& config, CoordInfo& coordInfo, const std::vector<std::string>& topOrder,
        const std::string& top_fname, std::vector<std::array<float, 3>>& proteinCarbons, std::unordered_map<std::string, int>& vaporDict);

    std::vector<int> computeClusters(const std::vector<std::array<float, 3>>& points);
    
    std::unordered_map<int, int> findClusterWaters(const std::vector<int>& clusterIDs);

    std::vector<int> pinCluster(const std::vector<std::shared_ptr<Atom>>& titAtoms, const std::vector<std::array<float, 3>>& waterCoords,
        const std::vector<int>& clusterIDs, const std::unordered_map<int, int>& clusterWaters, const bool& donor);

    std::vector<int> pinClusterWaters(const std::vector<int>& clusterIDs, const std::unordered_map<int, Cluster>& clusters);

    std::unordered_map<int, Cluster> binClusters(const std::vector<int>& clusterIDs, const CoordInfo& coordInfo, const Config& config);
    
    TitratableSites getTitratableSites(const CoordInfo& coordInfo, const std::vector<int>& clusterIDs,
        const std::unordered_map<int, Cluster>& clusters, const std::unordered_map<int, float>& pkaMap);

    std::vector<std::string> getProteinResNames(const CoordInfo& coordInfo);

    float grotthussEnergy(const CoordInfo& coordInfo, const std::vector<float>& charges, const TitratableAtom& donor,
        const TitratableAtom& acceptor, const int& hop);

    std::vector<Exchange> findPairs(const CoordInfo& coordInfo, const TitratableSites& titSites,
        const std::vector<int>& clusterIDs, const std::unordered_map<int, Cluster>& clusters, const float& temperature,
        const int& step, const int& hop, const std::vector<float>& charges, std::vector<std::array <float, 3>>& skipCoords,
        bool& pairs, bool& prot);

    float getAngle(const std::array<float, 3>& p1, const std::array<float, 3>& p2, const std::array<float, 3>& p3);

    float convertToRadians(const float& degrees);

    std::array<float, 3> addVectors(const std::array<float, 3>& u, const std::array<float, 3>& v);

    std::array<float, 3> subtractVectors(const std::array<float, 3>& u, const std::array<float, 3>& v);

	std::array<float, 3> scaleVector(const std::array<float, 3>& u, const float& scalar);

	std::array<float, 3> crossProduct(const std::array<float, 3>& u, const std::array<float, 3>& v);

    std::array<float, 3> normalizeVector(const std::array<float, 3>& vec);

    float dotProduct(const std::array<float, 3>& u, const std::array<float, 3>& v);

    std::array<float, 3> planeNormalVector(const std::array<float, 3>& p1, const std::array<float, 3>& p2,
        const std::array<float, 3>& p3, const std::array<float, 3>& p4);

    std::array<float, 3> rotateBond(const std::array<float, 3>& p1, const std::array<float, 3>& p2,
        std::array<float, 3>& p3, const float& theta);

    std::vector<std::array<float, 3>> rotationMatrix(const std::array<float, 3>& axis, const float& theta);

    std::array<float, 3> dotProductforRot(const std::vector<std::array<float, 3>>& rotationMatrix, const std::array<float, 3>& v);

    void H3O_transform(std::vector<std::shared_ptr<Atom>>& atoms);

    void water_transform(std::vector<std::shared_ptr<Atom>>& atoms);

    void deleteAtoms(std::vector<std::shared_ptr<Atom>>& atoms, const std::vector<int>& indicesToDelete);

    std::vector<int> isGasPhase(const CoordInfo& coordInfo, const std::vector<std::array<float, 3>>& proteinCarbons,
        const std::string resType);

    void doExchanges(CoordInfo& coordInfo, const std::vector<Exchange>& exchanges, const bool& pairs, const bool& prot,
        const std::vector<std::string>& topOrder, const Config& config);

    bool recenterDroplet(CoordInfo& coordInfo, const std::vector<std::array<float, 3>>& proteinCoords);

    bool correctGasVelocities(std::vector<std::shared_ptr<Atom>>& atoms, float temperature);

    void fixAmmoniaGas(const Config& config, CoordInfo& coordInfo, bool& gasCorrect);

    void fixAceticGas(const Config& config, CoordInfo& coordInfo, bool& gasCorrect);

    void simulation(const Config& config, const std::filesystem::path& trialPath);
}