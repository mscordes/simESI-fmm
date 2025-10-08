#pragma once
#include <filesystem>
#include <fstream>
#include <string>

using namespace std;
namespace fs = filesystem;

class ParamInfoWriter {
    private:
        const fs::path fileName;

    public:
        explicit        ParamInfoWriter(const fs::path& fileName);
        virtual void    write() const;
};