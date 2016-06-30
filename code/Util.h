#ifndef Util_H
#define Util_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

namespace Util {

    extern std::string         getLineFromFile(std::string fileName, int lineNum);
    extern int                 flip(int x);
    extern std::string         getTime(void);
    
	extern std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
	extern std::vector<std::string> split(const std::string &s, char delim);
	extern std::string getTime(void);

	extern std::string boolToString(bool);
	extern std::string doubleToString(double);
	extern std::string intToString(int);

	extern int stringToInt(std::string);
	extern double stringToDouble(std::string);
	extern bool stringToBool(std::string);

    extern std::istream& safeGetline(std::istream& is, std::string& t);

}

#endif
