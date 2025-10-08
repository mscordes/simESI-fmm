#pragma once
#include <string>
#include <vector>

using namespace std;

void	trimRef(string& s);
string	trim(const string& s);
void	stripComment(string& line);
string	toLower(string);
void	cleanString(string& line);
bool	isOnlyWhitespace(const string&);
bool	isNumeric(const char*);
bool	isInteger(const string& s);
bool	isAlphabetic(const string& str);
bool	extractDoubleAfter(string, double&);
bool	extractIntAfter(string, int&);
bool	parseStringBefore(string&, string&);
bool	parseStringAfter(string&, string&);
bool	parseBoolAfter(string& input, bool& output);
vector<string> splitLine(const string& line);