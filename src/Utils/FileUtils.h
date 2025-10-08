#pragma once
#include <filesystem>
#include <string>

using namespace std;
namespace fs = filesystem;

void deleteFile(const fs::path& fname);
void createDir(const fs::path& dest);
void copyFile(const fs::path& source, const fs::path& dest);
void copyFiles(const fs::path& idir, const fs::path& odir);
void copyRecursive(const fs::path& source, const fs::path& dest);
void populateOutputDirs();
void renameCheckpointFile(const string& ftype, int step, int num_restarts);