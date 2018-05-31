#ifndef Model_H
#define Model_H

#include <set>
#include <vector>
class Areas;
class GeoCoords;
class MbRandom;
class Node;
class Settings;
class Tree;


// NOTE:
// sumOfRates() and transitionRate() allow different dispersal mechanisms to be modeled.



class Model {
    
    public:
                                Model(Areas* a, GeoCoords* g, Tree* t, MbRandom* r, Settings* s);
                               ~Model(void);
        double                  lnLikelihood(int space);
        double					lnLikelihoodForPath(Node* p, int space);
        double                  getGainRate(void) { return areaRates[1]; }
        double                  getLossRate(void) { return areaRates[0]; }
        double					getGainPathRate(void) { return pathRates[1]; }
        double					getLossPathRate(void) { return pathRates[0]; }
        double                  getDistancePower(void) { return distancePower; }
        double                  getCarryingCapacity(void) { return carryingCapacity; }
        double                  getCarryingPower(void) { return carryingPower; }
        double					getRangeSize(void)	{ return rangeSize; }
        int                     getNumAreas(void) { return numAreas; }
        int                     getNumChanges(int space);
        int                     getNumGains(int space);
        int                     getNumLosses(int space);
        int                     getNumGainsTree(int space);
        int                     getNumLossesTree(int space);
        Tree*                   getTreePtr(void) { return treePtr; }
        double                  getLnLike(void) { return lnLike; }
        void                    showConditionalLikelihoods(void);
        void                    showTransitionProbabilities(void);
        double                  sumOfRates(std::vector<bool>& x, int srModelType);
        double                  sumOfRates(int srModelType, double eventTime, std::vector<double>& rates, std::vector<bool>& x);
        double                  transitionRate(int trModelType, double eventTime, std::vector<double>& rates, std::vector<bool>& fromState, std::vector<bool>& toState, int toArea);
        double                  transitionRate(std::vector<bool>& fromState, std::vector<bool>& toState, int toArea, int trModelType, double t);

        void                    updateNode(int mcmcCycle, int na);
        void                    updateAreaRate(int mcmcCycle, double tuning);
        void					updatePathRate(int mcmcCycle, double tuning);
        void                    updateDistancePower(int mcmcCycle, double tuning);
        void                    updateCarryingPower(int mcmcCycle, double tuning);
        void                    updateCarryingCapacity(int mcmcCycle, double tuning);
        void					updateConsensus(int mcmcCycle, int na);
    
    private:
        double                  calcLnProposalRatio(void);
        void                    conditionalLikelhoodDown(void);
        void                    conditionalLikelhoodDown(std::set<int>& areaSet);
        bool                    drawAncestralAreas(Node* p, std::set<int>& areaSet);
        bool                    drawChanges(Node* p, std::set<int>& areaSet);
        void					drawConsensus(std::set<int>& areaSet);
        void                    inititalizeConditionalLikelihoods(void);
        void                    initializeHistory(void);
        double                  lnBinomialCoefficient(int n, int k);
        double                  lnFactorial(int n);
        double                  lnProbK(int n, int k);
        double                  lnProposedStateProb(int space);
        double					lnProposedStateProbForPath(Node* p, int space);
        double					lnProposedConsensusProb(std::set<int>& areaSet, int space);
        int                     numOnStates(std::vector<bool>& x);
        void                    initializeParameters(void);
    
        Areas*                  areasPtr;
        MbRandom*               ranPtr;
        Settings*               settingsPtr;
        Tree*                   treePtr;
        GeoCoords*              geoPtr;
        std::vector<double>		areaRateParm;
        std::vector<double>		pathRateParm;
        double					distancePowerParm;
        double                  carryingPowerParm;
        double                  carryingCapacityAlphaParm;
        double                  carryingCapacityBetaParm;
        double                  scaleFactor;
        double					rangeSize;
        std::vector<double>     areaRates;
        std::vector<double>		pathRates;
        double					distancePower;
        double					carryingCapacity;
        double					carryingPower;
        double                  lnLike;
        int                     numAreas;
        std::vector<std::vector<double> > distances;
        std::string             geoDistancePowerStr;

        bool					useSteppingStone;
        double					betaSteppingStone;
    
        bool                    useDistanceThreshold;
        double                  distanceThreshold;

        int                     modelType;
        bool                    useScaleFactor;
};

#endif
