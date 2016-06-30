#ifndef Mcmc_H
#define Mcmc_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Atlas.h"
#include "Areas.h"
#include "BranchHistory.h"
#include "GeoCoords.h"
#include "MbRandom.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"

class Areas;
class BranchHistory;
class GeoCoords;
class MbRandom;
class Model;
class Node;
class Settings;
class Tree;

class Mcmc {
    
    public:
                                Mcmc(Areas* a, GeoCoords* g, MbRandom* r, Model* m, Tree* t, Settings* s);
                               ~Mcmc(void);
    
    private:
        void                    runChain(void);
        void                    sampleParms(int n);
        void                    samplePaths(int n);
        void                    updateAncestralStateCounts(int n);
        void                    writeAncestralStateFreqs(void);
        void                    writeNhxFile(void);
        std::string             addNodeNhxToString(Node* p, std::string s);
        GeoCoords*              geoPtr;
        Areas*                  areasPtr;
        MbRandom*               ranPtr;
        Settings*               settingsPtr;
        Tree*                   treePtr;
        Model*                  modelPtr;
        std::vector<double>     proposalProbs;
        std::ofstream           outStrm;
        std::ofstream           histStrm;
        std::ofstream           areaPosteriorStrm;
        std::ofstream           nhxStrm;
        int                     burnIn;
        int                     numSamples;
};

#endif
