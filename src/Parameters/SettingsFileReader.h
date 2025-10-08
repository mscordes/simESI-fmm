#pragma once
#include <filesystem>

using namespace std;
namespace fs = filesystem;


class SettingsFileReader {
    private:
        fs::path fileName;

    public:
        SettingsFileReader(const fs::path&);
        virtual bool Read();
};