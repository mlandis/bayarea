#ifndef Areas_H
#define Areas_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

class Areas {

	public:
                                    Areas(std::string fileName);  
                                   ~Areas(void);
        int                         getNumTaxa(void) { return numTaxa; }
        int                         getNumAreas(void) { return numAreas; }
        int                         getNumOccuppiedInArea(int aIdx);
        int                         getTaxonIndex(std::string ns);
        std::vector<std::string>&   getTaxonList(void) { return taxonNames; }
        bool                        isInArea(int tIdx, int aIdx);
        void                        listTaxa(void);
        void                        print(void);
        std::string                 getTaxonName(int i);
        int                         areaState(int tIdx, int aIdx) { return matrix[tIdx][aIdx]; }

    private:
        int                         numTaxa;
        int                         numAreas;
        std::vector<std::string>    taxonNames;
        bool                        compressedData;
        int**                       matrix;
};

#endif
