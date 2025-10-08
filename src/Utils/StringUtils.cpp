#include "StringUtils.h"

#include <algorithm> 
#include <cctype> 
#include <cstring> 
#include <iostream>
#include <sstream>

void trimRef(string& s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(),
        [](unsigned char ch) { return !isspace(ch); }));
    s.erase(find_if(s.rbegin(), s.rend(),
        [](unsigned char ch) { return !isspace(ch); }).base(),
        s.end());
}

string trim(const string& s) {
    auto start = find_if(s.begin(), s.end(),
        [](unsigned char ch) { return !isspace(ch); });

    auto end = find_if(s.rbegin(), s.rend(),
        [](unsigned char ch) { return !isspace(ch); }).base();

    if (start >= end) {
        return ""; 
    }
    return string(start, end);
}

void stripComment(string& line) {
    size_t commentPos = line.find('#');
    if (commentPos != string::npos) {
        line = line.substr(0, commentPos);
    }
}

string toLower(string line) {
    for (int i = 0; i < line.size(); ++i) {
        line[i] = tolower(line[i]);
    }
    return line;
}

void cleanString(string& line) {
    line = toLower(line);
    stripComment(line);
    trimRef(line);
}

bool isOnlyWhitespace(const string& s) {
    return s.empty() &&
        all_of(s.begin(), s.end(),
            [](unsigned char c) { return isspace(c); });
}

bool isNumeric(const char* string) {
    for (int i = 0; i < strlen(string); ++i) {
        if (!isdigit(string[i]) && string[i] != '.' && string[i] != '-') {
            // Check scientific notation, #e# or #E#
            if (string[i] == 'e' || string[i] == 'E') {
                if (i == 0 || i == strlen(string) - 1) return false;
                if (isdigit(string[i - 1]) && isdigit(string[i + 1])) continue;
            }
            return false;
        }
    }
    return true;
}

bool isInteger(const string& s) {
    stringstream ss(s);
    int x;
    return (ss >> x) && (ss.eof());
}

bool isAlphabetic(const string& str) {
    for (char c : str) {
        if (!isalpha(static_cast<unsigned char>(c))) return false;
    }
    return true; 
}

bool extractDoubleAfter(string line, double& data) {
    size_t delim = line.find("=");
    if (delim != string::npos) {
        data = stod(line.substr(delim + 1, string::npos).c_str());
        return true;
    }
    return false;
}

bool extractIntAfter(string line, int& data) {
    size_t delim = line.find("=");
    if (delim != string::npos) {
        data = stoi(line.substr(delim + 1, string::npos).c_str());
        return true;
    }
    return false;
}

bool parseStringBefore(string& input, string& output) {
    size_t split = input.find('=');
    if (split != string::npos) {
        output = input.substr(0, split);
		trimRef(output);
        return true;
    }
    return false;
}

bool parseStringAfter(string& input, string& output) {
    size_t split = input.find('=');
    if (split != string::npos) {
        output = input.substr(split + 1, string::npos);
		trimRef(output);
        return true;
    }
    return false;
}

// Bool containing string must be 'yes' or 'no'
bool parseBoolAfter(string& input, bool& output) {
    string stringWithBool;
    if (parseStringAfter(input, stringWithBool)) {
        if (stringWithBool.find("yes") != string::npos) {
            output = true;
            return true;
        }
        else if (stringWithBool.find("no") != string::npos) {
            output = false;
            return true;
        }
        else {
            return false;
        }
    }
    return false;
}

vector<string> splitLine(const string& line) {
    vector<string> words;
    stringstream ss(line);
    string word;
    while (ss >> word) words.push_back(word);
    return words;
}