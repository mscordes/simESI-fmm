#pragma once
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <stdexcept>    

using namespace std;
namespace fs = filesystem;

class Writer {
    protected:
        fs::path filename;
        ofstream outFile;

    public:
        Writer(fs::path filename);
        ~Writer();
};