#include "Writer.h"  

Writer::Writer(fs::path filename_) : filename(move(filename_)) {
    outFile = ofstream(filename, ios::app);
    if (!outFile.is_open())
        throw runtime_error("Could not open Writer for " + filename.string() + ".");
}

Writer::~Writer() {
    if (outFile.is_open()) {
        outFile.flush();
        outFile.close();
    }
}