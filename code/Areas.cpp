#include "Areas.h"
#include "Util.h"

Areas::Areas(std::string fileName) {

	/* open the file */
	std::ifstream seqStream(fileName.c_str());
	if (!seqStream) 
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		std::exit(1);
		}

	std::string linestring = "";
	int line = 0;
	std::string theSequence = "";
	int taxonNum = 0;
	matrix = NULL;
	numTaxa = numAreas = 0;
	while( !Util::safeGetline(seqStream, linestring).eof() )
		{
		std::istringstream linestream(linestring);
		int ch;
		std::string word = "";
		int wordNum = 0;
		int siteNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			//cout << "word(" << wordNum << ") = " << word << endl;
			if (line == 0)
				{
				/* read the number of taxa/chars from the first line */
				int x;
				std::istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numAreas = x;
				if (numTaxa > 0 && numAreas > 0 && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numAreas];
					for (int i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numAreas;
					for (int i=0; i<numTaxa; i++)
						for (int j=0; j<numAreas; j++)
							matrix[i][j] = 0;
					}
				}
			else
				{
				if (wordNum == 1 && word != "")
					{
                    taxonNames.push_back(word);
                    taxonNum++;
					}
				else
					{
                    for (int i=0; i<word.length(); i++)
                        {
                        char area = word.at(i);
                        if (area == '0')
                            matrix[taxonNum-1][siteNum++] = 0;
                        else if (area == '1')
                            matrix[taxonNum-1][siteNum++] = 1;
                        }
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		//cout << linestring << endl;
		line++;
		}	
	
	//numTaxa = 0;
	//numTaxa++;

	/* close the file */
	seqStream.close();

}

Areas::~Areas(void) {

	delete [] matrix[0];
	delete [] matrix;
}

bool Areas::isInArea(int tIdx, int aIdx) {
    
    if (matrix[tIdx][aIdx] == 1)
        return true;
    return false;
}

void Areas::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << std::endl;

}

int Areas::getNumOccuppiedInArea(int aIdx) {
    
    int n = 0;
    for (int i=0; i<numTaxa; i++)
        {
        if ( isInArea(i, aIdx) == true )
            n++;
        }
    return n;
}

std::string Areas::getTaxonName(int i) {

	return taxonNames[i];

}

int Areas::getTaxonIndex(std::string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;

}

void Areas::print(void) {

    listTaxa();
    std::cout << std::endl;
	int **x;
    x = matrix;
	std::cout << "        ";
	for (int i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << std::endl;
	std::cout << "--------";
	for (int i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << std::endl;	
	for (int j=0; j<numAreas; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (int i=0; i<numTaxa; i++)
			{
			std::cout << std::setw(3) << x[i][j];
			}
		std::cout << std::endl;
		}
}


