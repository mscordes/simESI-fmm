#include "Core.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>

namespace Core {

    // Read a file line by line
    std::vector<std::string> readFile(const std::string& fileName) {
        std::vector<std::string> lines;  // Vector to store the lines
        std::ifstream file(fileName);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open the file " << fileName << std::endl;
            return lines;  // Return an empty vector if the file can't be opened
        }

        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);  // Add each line to the vector
        }

        file.close();
        return lines;  // Return the vector containing all lines
    }

    // Function to split a single line into words
    std::vector<std::string> splitLine(const std::string& line) {
        std::vector<std::string> words;
        std::stringstream ss(line);  // Create a stringstream from the line
        std::string word;

        // Split the line by spaces and add words to the vector
        while (ss >> word) {
            words.push_back(word);
        }

        return words;
    }

    // Function to split each line in the vector into a vector of words
    std::vector<std::vector<std::string>> splitLines(const std::vector<std::string>& lines) {
        std::vector<std::vector<std::string>> allWords;

        for (const std::string& line : lines) {
            allWords.push_back(splitLine(line));  // Split each line into words
        }

        return allWords;
    }

    void copyFile(const std::filesystem::path& ifname, const std::filesystem::path& ofname) {
        std::filesystem::copy_file(ifname, ofname, std::filesystem::copy_options::overwrite_existing);
    }

    // Copies all files from one dir to another, does not overwrite!
    void copyFiles(const std::filesystem::path& idir, const std::filesystem::path& odir) {
		if (!std::filesystem::exists(odir)) {
			std::cerr << "Output dir for copying does not exist." << std::endl;
            std::exit(1);
		}
        if (!std::filesystem::exists(idir)) {
            std::cerr << "Input dir for copying does not exist." << std::endl;
            std::exit(1);
        }

		for (const auto& entry : std::filesystem::directory_iterator(idir)) {
            if (!std::filesystem::exists(odir / entry.path().filename())) {
                std::filesystem::copy_file(entry.path(), odir / entry.path().filename());
            }
		}
    }

    // Recursively copies all files and folders from src to target and overwrites existing files in target.
    void copyRecursive(const std::filesystem::path& src, const std::filesystem::path& target) noexcept
    {
        try
        {
            std::filesystem::copy(src, target, std::filesystem::copy_options::overwrite_existing | std::filesystem::copy_options::recursive);
        }
        catch (std::exception& e)
        {
            std::cout << e.what();
        }
    }

	// Deletes a file
    void deleteFile(const std::filesystem::path& fname) {
        try {
            if (std::filesystem::exists(fname)) { 
                std::filesystem::remove(fname); 
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Filesystem error when deleteing file. " << e.what() << std::endl;
            std::exit(1);
        }
        catch (const std::exception& e) {
            std::cerr << "General error when deleteing file. " << e.what() << std::endl;
        }
    }

    // Write to output file in data dir
    void write_output(const std::filesystem::path& target_dir, const std::string& input, const std::string& output_fname) {
        // Combine target directory and file name to form the path
        std::filesystem::path file_path = target_dir / output_fname;

        // Open the file in append mode
        std::ofstream file(file_path, std::ios_base::app);

        // Check if the file was opened successfully
        if (file.is_open()) {
            file << input << '\n';
            file.close();
        }
        else {
            std::cerr << "Failed to open file: " << file_path << std::endl;
        }
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> timer() {
        return std::chrono::high_resolution_clock::now();
    }

    std::string duration(const std::chrono::time_point<std::chrono::high_resolution_clock>& start,
        const std::chrono::time_point<std::chrono::high_resolution_clock>& end) {
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4) << elapsed_seconds.count();
        return oss.str();
    }

    // Renames files to account for GROMACS checkpoint naming
	void renameCheckpointFile(const std::string& ftype, const int& step, const int& num_restarts) {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << num_restarts;

        // Construct source and destination file paths
        std::string source = std::to_string(step) + ".part" + oss.str() + ftype;
        std::string destination = std::to_string(step) + ftype;

        try {
            if (std::filesystem::exists(source)) {
                std::filesystem::rename(source, destination);
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Filesystem error when renaming file. " << e.what() << std::endl;
			std::exit(1);
        }
	}


}