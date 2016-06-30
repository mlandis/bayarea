#include "Util.h"

std::vector<std::string>& Util::split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> Util::split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	return split(s, delim, elems);
}


std::string Util::getLineFromFile(std::string fileName, int lineNum) {

	/* open the file */
	std::ifstream fileStream(fileName.c_str());
	if (!fileStream)
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	//while( Util::safeGetline(fileStream, linestring).good() )
    while( !Util::safeGetline(fileStream, linestring).eof() )
		{
		line++;
		if (line == lineNum)
			break;
		}

	/* close the file */
	fileStream.close();

	if (line != lineNum)
		{
		std::cerr << "The file \"" + fileName + "\" has " << line << " lines. Could not find line " << lineNum << std::endl;
		exit(1);
		}

	return linestring;
}

int Util::flip(int x) {

	if (x == 0)
		return 1;
	else
		return 0;

}


std::string Util::getTime(void)
{
	std::string s;

	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	// two-digit representation
	// YYMMDDhhmmss
	strftime (buffer,80,"%y%m%d%H%M%S",timeinfo);

	s = buffer;

	return s;
}

bool Util::stringToBool(std::string s)
{
	if (s == "True" || s == "T" || s == "1" || s == "TRUE")
		return true;
	else return false;
}


std::string Util::boolToString(bool tf)
{
	if (tf) return "True";
	else return "False";
}

std::string Util::intToString(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

std::string Util::doubleToString(double number)
{
	std::stringstream ss;
	ss << number;
	return ss.str();
}

int Util::stringToInt(std::string s)
{
	return atoi(s.c_str());
}

double Util::stringToDouble(std::string s)
{
	return atof(s.c_str());
}

std::istream& Util::safeGetline(std::istream& is, std::string& t)
{
    t.clear();
    
    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.
    
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();
    
    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}
