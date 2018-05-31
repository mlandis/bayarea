#include "Areas.h"
#include "AreaChange.h"
#include "BranchHistory.h"
#include "ConditionalLikelihood.h"
#include "GeoCoords.h"
#include "MbRandom.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "Settings.h"
#include "TransitionProbability.h"
#include "Tree.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <pthread.h>

#define INDEPENDENCE 1
#define CARRYING_CAPACITY 2
#define DISTANCE_NORM 3
#define CARRYING_CAPACITY_AND_DISTANCE_NORM 4

#define MAX_EVENTS 10000

#define DEBUG_MH 0
#define DEBUG_PATH 0
#define DEBUG_PROPOSE_RATE 0
#define DEBUG_PROPOSE_NODE_PATH 0
#define DEBUG_PROPOSE_NODE 0
#define DEBUG_PROPOSE_BRANCH 0
#define DEBUG_PROPOSE_DISTANCE_POWER 0
#define DEBUG_PROPOSE_CARRYING_CAPACITY 0
#define DEBUG_PROPOSE_CARRYING_POWER 0
#define DEBUG_DRAW_CHANGES 0
#define DEBUG_DRAW_CHANGES2 0
#define DEBUG_DISTANCE_TR_RATE 0
#define DEBUG_LIKELIHOOD 0
#define DEBUG_LIKELIHOOD2 0
#define DEBUG_DISTANCE_SUM_RATE 0

Model::Model(Areas* a, GeoCoords* g, Tree* t, MbRandom* r, Settings* s) {
    
    // remember the address of important objects
    areasPtr    = a;
    ranPtr      = r;
    settingsPtr = s;
    treePtr     = t;
    geoPtr      = g;
    
    // initialize area rate values (per-sequence)
    numAreas     = areasPtr->getNumAreas();
    areaRateParm.resize(2);
    areaRateParm[0] = settingsPtr->getAreaLossPrior();
    areaRateParm[1] = settingsPtr->getAreaGainPrior();
    areaRates.resize(2);
    areaRates[0] = ranPtr->exponentialRv(areaRateParm[0]);
    areaRates[1] = ranPtr->exponentialRv(areaRateParm[1]);
    //areaRates[0] = 0.05;
    //areaRates[1] = 0.5;
	//areaRates[0] = settingsPtr->getAreaLossRate(); //0.1;	// rate of change 1->0
	//areaRates[1] = settingsPtr->getAreaGainRate(); //0.1;	// rate of change 0->1

    // initialize path proposal rate values (per-site)
    pathRateParm.resize(2);
    pathRateParm[0] = settingsPtr->getAreaLossPrior();
    pathRateParm[1] = settingsPtr->getAreaGainPrior();
    pathRates.resize(2);
    pathRates[0] = ranPtr->exponentialRv(pathRateParm[0]);
    pathRates[1] = ranPtr->exponentialRv(pathRateParm[1]);
    //pathRates[0] = 0.05;
    //pathRates[1] = 0.5;
    //pathRates[0] = areaRates[0];
    //pathRates[1] = areaRates[1];

    // initialize carrying capacity values
    // will want to make this a free parameter later
    carryingPowerParm = 1.0;
    carryingCapacityAlphaParm = 2.0;
    carryingCapacityBetaParm = 2.0;
    carryingCapacity = ranPtr->betaRv(1, 1);
    carryingPower = ranPtr->exponentialRv(carryingPowerParm);
    //carryingCapacity = settingsPtr->getCarryingCapacity();
    //carryingPower = settingsPtr->getCarryingPower();


    // initialize distance power parameter
    distancePowerParm = 1.0;
    distancePower = fabs(ranPtr->normalRv(0.0, distancePowerParm)) / 10.0;
    //distancePower = 0.0;
    //distancePower = settingsPtr->getGeoDistancePower();

    // initialize stepping stone settings
    useSteppingStone = settingsPtr->getUseSteppingStone();
    betaSteppingStone = settingsPtr->getBetaSteppingStone();

    // initialize distance threshold settings
    useDistanceThreshold = false;
    distanceThreshold = 1e-3;
    
    // initialize model settings
    modelType    = settingsPtr->getModelType();
    
    geoDistancePowerStr = (settingsPtr->getGeoDistanceTruncate() ? "haversinePowerInverseTruncated" : "haversinePowerInverse");

    // load pre-calculated distance map
    //geoPtr->changeDistancePower("haversinePowerInverse", 0.0, distancePower);
    //distances = geoPtr->getDistance("haversinePowerInverse", 0.0);
    //geoPtr->changeDistancePower("haversinePowerInverseTruncated", 0.0, distancePower);
    //distances = geoPtr->getDistance("haversinePowerInverseTruncated", 0.0);
    geoPtr->changeDistancePower(geoDistancePowerStr, 0.0, distancePower);
    distances = geoPtr->getDistance(geoDistancePowerStr, 0.0);
    
    
    // use parismony-style approximation to initalize parameters
    if (settingsPtr->getGuessInitialRates())
        initializeParameters();
    
    //initializeOrderedDistanceIdx();
    // std::cout << geoPtr->getDistanceStr("haversinePowerInverse", 0.0) << "\n";
    
    // attach a vector of conditional likelihoods to each node in the tree
    inititalizeConditionalLikelihoods();
    
    // initialize history
    initializeHistory();
    
    // initialize the likelihood
    lnLike = lnLikelihood(0);
    if (useSteppingStone == true)
    	lnLike *= betaSteppingStone;
    
    //std::cout << "treeLength = " << treePtr->getTreeLength() << "\n";
    std::cout << "Initial log likelihood:\t" << lnLike << std::endl;
}

Model::~Model(void) {
}

double Model::calcLnProposalRatio(void) {
    
    // we assume that we go from 0 -> 1 and that 1 is always the proposed state
    double lnProb0 = lnProposedStateProb(0);
    double lnProb1 = lnProposedStateProb(1);
    return lnProb0 - lnProb1;
}

void Model::conditionalLikelhoodDown(void) {
    
    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->isTip() == false )
            {
            double* clP = p->getConditionalLikelihood()->getCls();
            std::vector<double*> clD;
            std::vector<double**> tiD;
            for (int i=0; i<p->getNumDescendants(); i++)
                {
                double* x = p->getDescendantIndexed(i)->getConditionalLikelihood()->getCls();
                clD.push_back(x);
                double** y = p->getDescendantIndexed(i)->getTransitionProbability()->getTiProb();
                tiD.push_back(y);
                }
            
            for (int a=0; a<areasPtr->getNumAreas(); a++)
                {
                for (int i=0; i<2; i++)
                    {
                    double prob = 1.0;
                    for (int d=0; d<p->getNumDescendants(); d++)
                        {
                        double linProb = 0.0;
                        for (int j=0; j<2; j++)
                            {
                            linProb += clD[d][j] * tiD[d][i][j];
                            }
                        prob *= linProb;
                        }
                    clP[i] = prob;
                    }
                
                clP += 2;
                for (int d=0; d<p->getNumDescendants(); d++)
                    clD[d] += 2;
                }
            
            }
        }
}

void Model::conditionalLikelhoodDown(std::set<int>& areaSet) {
    
    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->isTip() == false )
            {
            double* clP_start = p->getConditionalLikelihood()->getCls();
            double* clP = p->getConditionalLikelihood()->getCls();
            std::vector<double*> clD;
            std::vector<double*> clD_start;
            std::vector<double**> tiD;
            for (int i=0; i<p->getNumDescendants(); i++)
                {
                double* x = p->getDescendantIndexed(i)->getConditionalLikelihood()->getCls();
                clD.push_back(x);
                clD_start.push_back(x);
                double** y = p->getDescendantIndexed(i)->getTransitionProbability()->getTiProb();
                tiD.push_back(y);
                }
            
            for (std::set<int>::iterator it = areaSet.begin(); it != areaSet.end(); it++)
                {
                int a = (*it);
                for (int i=0; i<2; i++)
                    {
                    double prob = 1.0;
                    for (int d=0; d<p->getNumDescendants(); d++)
                        {
                        double linProb = 0.0;
                        for (int j=0; j<2; j++)
                            {
                            linProb += clD[d][j] * tiD[d][i][j];
                            }
                        prob *= linProb;
                        }
                    clP[i] = prob;
                    }
                
                clP = clP_start + 2*a;
                for (int d=0; d<p->getNumDescendants(); d++)
                    clD[d] = clD_start[d] + 2*a;
                }
            
            }
        }
}

bool Model::drawAncestralAreas(Node* p, std::set<int>& areaSet) {
    
    // get the vector of descendants of node p
    std::vector<Node*> descendants = p->getDescendants();

    // loop over the areas updating the ancestral state a node p for each
    for (std::set<int>::iterator a = areaSet.begin(); a != areaSet.end(); a++)
    {
    	// initialize the probability for the area
        double prob[2];
        BranchHistory* h = p->getBranchHistory(1);
        int startState = h->getAncestralStateForArea(*a);
        double** t = p->getTransitionProbability()->getTiProb();
        prob[0] = t[startState][0];
        prob[1] = t[startState][1];
        
        // account for the factors that occur along the branches leading to the descendants
        for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
            // get the transition probability matrix for this branch
            double** t = (*it)->getTransitionProbability()->getTiProb();
            
            // also, get the ending state for the branch
            int endState = 0;
            if ( (*it)->getNumDescendants() == 0 )
            {
                endState = areasPtr->areaState( (*it)->getIndex(), (*a) );
            }
            else
            {
                std::vector<bool> s = (*it)->getDescendantIndexed(0)->getBranchHistory(1)->getAncestralState();
                endState = s[(*a)];
            }
            
            // and deal with the transition probability along the branch
            prob[0] *= t[0][endState];
            prob[1] *= t[1][endState];
        }
        
        // calculate the probability of being in each of the two states for node p using Bayes formula and
        // choose a state using these probabilities
        double sum = prob[0] + prob[1];
        int ancState = 0;
        if ( ranPtr->uniformRv() > prob[0] / sum )
        //if ( ranPtr->uniformRv() > carryingCapacity * sum / prob[0])
            ancState = 1;

        // set the state remembering that it is stored in the histories for the descendants
        for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
            BranchHistory* h = (*it)->getBranchHistory(1);
            h->setAncestralArea( (*a), ancState );
        }
    }

    // resample consensus sequence if root
    if (p->getAncestor() == NULL) {
        for (std::set<int>::iterator a = areaSet.begin(); a != areaSet.end(); a++)
        {
            double u = ranPtr->uniformRv();
            if (u < areaRates[0] / (areaRates[0] + areaRates[1]))
                p->getBranchHistory(1)->setAncestralArea(*a,0);
            else
                p->getBranchHistory(1)->setAncestralArea(*a,1);
        }
    }

    if (p->getNumDescendants() == 0)
       ;// std::cout << p->getIndex() << "\t has 0 descendants\n";

    std::vector<bool> s1 = p->getDescendantIndexed(0)->getBranchHistory(1)->getAncestralState();
    if (numOnStates(s1) == 0)
    {
        //std::cout << "redraw node\n";
        //std::cout << numOnStates(s1);

        // restore node
        std::vector<bool> s0 = p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralState();
        for (std::set<int>::iterator a = areaSet.begin(); a != areaSet.end(); a++)
        {
            int ancState = s0[*a];
            for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
            {
                BranchHistory* h = (*it)->getBranchHistory(1);
                h->setAncestralArea( (*a), ancState );
            }
        }
        s1 = p->getDescendantIndexed(0)->getBranchHistory(1)->getAncestralState();

        // sample failed
        return false;
    }
    
    // sample succeeded
    return true;
}

bool Model::drawChanges(Node* p, std::set<int>& areaSet) {
    
   // std::cout << "start drawChanges\n";
    // get the length of the branch in absolute time
//	double pathRatesSum = pathRates[0] + pathRates[1];
    double len = p->getLen();

#if DEBUG_DRAW_CHANGES
    std::cout << "\ndrawChanges()\t" << p->getIndex() << "\t" << p->getLen() << "\n";
    std::cout << pathRates[0] << "\t" << pathRates[1] << "\n";
    std::cout << "areaSet\t";
    for (std::set<int>::iterator a = areaSet.begin(); a != areaSet.end(); a++)
        std::cout << *a << "\t";
    std::cout << "\n";

#endif

    // get the history
    BranchHistory* bh = p->getBranchHistory(1);
    
    // get the ancestor area configuration for the branch (for determining the start state)
    std::vector<bool> ancState = bh->getAncestralState();

    // loop over areas to update
    for (std::set<int>::iterator a = areaSet.begin(); a != areaSet.end(); a++)
    {
        // get the beginning and ending states for the branch
        int startState = (int)ancState[*a];
        int endState;
        if (p->getNumDescendants() == 0)
        {
            endState = areasPtr->areaState(p->getIndex(), *a);
        }
        else
        {
            endState = p->getDescendantIndexed(0)->getBranchHistory(1)->getAncestralState()[*a];
        }
        
        // sample the proposed branch history, reject improper histories
        int curState = startState;
//        int sampleCount = 0;
        std::vector<double> times;

        do
        {
            times.clear(); // why would I comment this?
            double t = 0.0;
            curState = startState;

#if DEBUG_DRAW_CHANGES
            sampleCount++;
#endif

            while (t < len)
            {

                // calculate rate
                double lambda = pathRates[1];
                if (curState == 1)
                    lambda = pathRates[0];
                
                // exponential waiting time
                t += ranPtr->exponentialRv(lambda);
                
                // potential change
                if (t < len)
                {
                    times.push_back(t / len); // unit time (0,1)
                    if (curState == 0)
                        curState = 1;
                    else
                        curState = 0;

#if DEBUG_DRAW_CHANGES

                    std::cout << "\t\tt: " << t * len << "\tlen: " << len << "\tsampleCount: " << sampleCount << "\tlambda: " << lambda << "\t" << *a << "\t" << (curState == 0 ? "1->0" : "0->1") << "\n";
#endif
                }
            }

            //if (curState != endState) std::cout << "rejected: curState != endState\n";

        } while(curState != endState);
        
        // add the changes to the history
        for (std::vector<double>::iterator it = times.begin(); it != times.end(); it++)
        {
            //std::cout << "addMe\t" << *it * p->getLen() << "\t" << *a << "\n";
            AreaChange* c = new AreaChange(*a, *it);
            bh->addAreaChange( c );
        }
    }
    //std::cout << "drawChanges check\n";
    // final check to make sure path does not include extinction configuration
    //bh = p->getBranchHistory(1);
    bool redrawChanges = false;

#if 0
    if (p->getAncestor() != NULL)
    {
        bh = p->getBranchHistory(1);
        std::set<AreaChange*,comp_history> changes = bh->getBranchChanges();
        std::vector<bool> s = bh->getAncestralState();

/*
        std::cout << "0.0000* \t";
        for (int i = 0; i < s.size(); i++)
        {
            std::cout << s[i];
        }
        std::cout << "\n";
*/

        for (std::set<AreaChange*,comp_history>::iterator it = changes.begin(); it != changes.end(); it++)
        {
            s[(*it)->getArea()].flip();


/*
            std::cout << (*it)->getPosition() * p->getLen() << "\t";
            for (int i = 0; i < s.size(); i++)
            {
                std::cout << s[i];
            }
            std::cout << "\n";
*/

            //std::cout << (*it)->getPosition() * p->getLen() << "\t" << (*it)->getArea() << "\t" << numAreas << "\t" << numOnStates(s) << "\n";
            if (numOnStates(s) == 0)
            {
                //std::cout << "numOnStates(s) == 0\n";
                redrawChanges = true;
                //it = changes.end();
                //break; // completely unexpected behavior! break should only exit out of immediate loop, but this exits out of random routines downstream!!!
                it = changes.end();
                it++;
            }
        }

/*
        if (p->getNumDescendants() > 0)
        {
            std::cout << p->getLen() << "* \t";
            std::vector<bool> sd = p->getDescendantIndexed(0)->getBranchHistory(1)->getAncestralState();
            for (int i = 0; i < sd.size(); i++)
            {
                std::cout << sd[i];
            }
            std::cout << "\n";

        }
        else
        {
            std::cout << p->getLen() << "* \t";
            for (int i = 0; i < numAreas; i++)
                std::cout << areasPtr->areaState(p->getIndex(), i);
            std::cout << "\n";
        }
*/

        /*
        if (redrawChanges == true)
        {
           // std::cout << "redraw path\n";
            bh->removeChanges(areaSet);
           // std::cout << "remove changes\n";
            drawChanges(p, areaSet);
        }*/

    }
#endif

#   if DEBUG_DRAW_CHANGES2
    std::vector<Node*> descendants = p->getDescendants();
    if (p->getAncestor() == NULL)
        {
        BranchHistory* h0 = p->getBranchHistory(0);
        BranchHistory* h1 = p->getBranchHistory(1);
        std::cout << "History along root branch before: " << std::endl;
        h0->print();
        std::cout << "History along root branch after: " << std::endl;
        h1->print();
        }
    /*
    for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        BranchHistory* h0 = (*it)->getBranchHistory(0);
        BranchHistory* h1 = (*it)->getBranchHistory(1);
        std::cout << "History along root branch before: " << std::endl;
        h0->print();
        std::cout << "History along root branch after: " << std::endl;
        h1->print();
        }
    */
#   endif

    return redrawChanges;
}

void Model::drawConsensus(std::set<int>& areaSet) {

	Node* p = treePtr->getRoot();
	std::vector<bool> consensusAncestor;
	int nOn = 0;
	do
	{

		// consensus ancestral area is occupied w.p. of tip occupancy frequency
		consensusAncestor.clear();
		for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
		{
			int nOccuppied = areasPtr->getNumOccuppiedInArea(aIdx);
			double prob = (double)nOccuppied / areasPtr->getNumTaxa();
			if (ranPtr->uniformRv() < prob)
				consensusAncestor.push_back( true );
			else
				consensusAncestor.push_back( false );
		}
		nOn = 0;
		for (int i=0; i<consensusAncestor.size(); i++)
		{
			if (consensusAncestor[i] == true)
				nOn++;
		}
	} while (nOn == 0);

	for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
	{
        int startState = 0;
//        int endState = 0;
		if (consensusAncestor[aIdx] == true)
			startState = 1;
		//if ( p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralStateForArea(aIdx) == true )
		//    endState = 1;

		p->getBranchHistory(1)->setAncestralArea(aIdx, startState);
	}
}

double Model::lnProposedConsensusProb(std::set<int>& areaSet, int space) {

	Node* p = treePtr->getRoot();
	BranchHistory* bh = p->getBranchHistory(space);
	std::vector<bool> s = bh->getAncestralState();

	double lnP = 0.0;

	for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
	{
		int nOccuppied = areasPtr->getNumOccuppiedInArea(aIdx);
		double prob = (double)nOccuppied / areasPtr->getNumTaxa();
		if (s[aIdx] == true)
			lnP += log(prob);
		else
			lnP += log(1 - prob);
	}

	return lnP;
}

int Model::getNumChanges(int space) {
    
    int nChanges = 0;
    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->getAncestor() != NULL )
            {
            BranchHistory* h = p->getBranchHistory(space);
            nChanges += h->getNumChanges();
            }
        }
    return nChanges;
}

int Model::getNumGains(int space) {

    int numGains = 0;
    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->getAncestor() != NULL )
            {
            BranchHistory* h = p->getBranchHistory(space);
            numGains += h->getNumGains();
            }
        }
    return numGains;
}

int Model::getNumLosses(int space) {
    
    int numLosses = 0;
    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->getAncestor() != NULL )
            {
            BranchHistory* h = p->getBranchHistory(space);
            numLosses += h->getNumLosses();
            }
        }
    return numLosses;
}


int Model::getNumGainsTree(int space){
	int numGains = 0;
	for (int n=0; n<treePtr->getNumNodes(); n++)
	{
		Node* p = treePtr->getDownPassNode(n);
		if ( p->getAncestor() != NULL )
		{
			BranchHistory* h = p->getBranchHistory(space);
			numGains += h->getNumGains();
		}
	}
	return numGains;
}

int Model::getNumLossesTree(int space) {
	int numLosses = 0;
	for (int n=0; n<treePtr->getNumNodes(); n++)
	{
		Node* p = treePtr->getDownPassNode(n);
		if ( p->getAncestor() != NULL )
		{
				BranchHistory* h = p->getBranchHistory(space);
				numLosses += h->getNumLosses();
		}
	}
	return numLosses;
}

void Model::initializeParameters(void) {

    // Current slapdash parameter initialization will underestimate the rates because:
    // 1. does not account for repeated substitutions (could probably estimate this)
    // 2. does not account for asymmetry of rates when estimating sum of rates
    // ... MLE implementation would certainly perform better, but requires more code.

    // 1. approx the ratio of rates
    std::vector<bool> s;
    s.resize(areasPtr->getNumAreas());
    int numOn = 0;
    for (int i = 0; i < areasPtr->getNumTaxa(); i++)
    {
        for (int j = 0; j < areasPtr->getNumAreas(); j++)
        {
            if (areasPtr->isInArea(i,j) == true)
            {
                numOn++;
            }
        }
    }
    double f = (double)numOn / (areasPtr->getNumTaxa() * areasPtr->getNumAreas());
    //std::cout << "numOn:\t" << numOn << "\n";
    //std::cout << "numTaxa:\t" << areasPtr->getNumTaxa() << "\tnumAreas:\t" << areasPtr->getNumAreas() << "\n";
    //std::cout << "f:\t" << f << "\n";

    // 2. rough approx the minimum number of changes observed from independent contrasts
    int count = 0;
    double v = 0.0;

    // examine tips only
    for (int i = 0; i < treePtr->getNumTaxa(); i++)
    {
        Node* p = treePtr->getDownPassNode(i);
        Node* q = p->getAncestor()->getDescendants()[1];
        if (p != q && p->isTip() && q->isTip())
        {
            v += p->getLen() + q->getLen();

            // get count of changes between tips
            for (int j = 0; j < areasPtr->getNumAreas(); j++)
            {
                if (areasPtr->isInArea(p->getIndex(), j) != areasPtr->isInArea(q->getIndex(), j))
                {
                    count++;
                }
            }
        }
    }
    double r = (double)count / (v * areasPtr->getNumAreas());
    //std::cout << "v:\t" << v << "\n";
    //std::cout << "r:\t" << r << "\n";

    pathRates[1] = f * 2 * r;
    areaRates[1] = f * 2 * r;
    pathRates[0] = (1 - f) * 2 * r;
    areaRates[0] = (1 - f) * 2 * r;
    distancePower = fabs(ranPtr->normalRv(0.0, distancePowerParm));
    
    std::cout << "Heuristic parameter initialization\n";
    std::cout << "\tArea gain rate = " << pathRates[1] << "\n";
    std::cout << "\tArea loss rate = " << pathRates[0] << "\n";
    std::cout << "\tDistance power = " << distancePower << "\n";
    
}

void Model::inititalizeConditionalLikelihoods(void) {
    
    for (int i=0; i<treePtr->getNumNodes(); i++)
        {
        Node* p = treePtr->getDownPassNode(i);
        
        // for each Node, assign an empty ConditionalLikelihood
        ConditionalLikelihood* c = new ConditionalLikelihood( p, areasPtr->getNumAreas(), 2 );
        p->attachConditionalLikelihoods(c);
        
        // for each Node, assign TransitionProbability according to branch length in units of expected subst. per Area
        // (note: root Node gets a long tail to approx. stationary distr.)
        TransitionProbability* tp = new TransitionProbability( p );
        p->attachTransitionProbability(tp);
       // if ( p->getAncestor() != NULL)
        tp->setProbs( p->getLen(), pathRates[1], pathRates[0], numAreas );
 //       else
   //     	tp->setProbs( treePtr->getTreeLength(), pathRates[1], pathRates[0], numAreas );
            //tp->setProbs( 10.0  / (pathRates[1] + pathRates[0]), pathRates[1], pathRates[0], numAreas );
        
        // for each Node, assign an empty BranchHistory
        BranchHistory* bh0 = new BranchHistory( p, areasPtr->getNumAreas() );
        BranchHistory* bh1 = new BranchHistory( p, areasPtr->getNumAreas() );
        (*bh1) = (*bh0);
        p->attachBranchHistory(bh0, 0);
        p->attachBranchHistory(bh1, 1);
        
        // for each tip Node, update ConditionalLikelihood according to observed Area data
        if ( p->isTip() == true )
            {
            double* clPtr = c->getCls();
            for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
                {
                bool inArea = areasPtr->isInArea( p->getIndex(), aIdx ); 
                if (inArea == true)
                    clPtr[1] = 1.0;
                else
                    clPtr[0] = 1.0;
                clPtr += 2;
                }
            }
        }
}

void Model::initializeHistory(void) {
    
	//std::cout << "MODEL: Initializing branch histories\n";
    // fill in ancestral states for each character

    // step 1: fill in conditional likelihoods down the tree
	// for each internal Node, for each Area, update ConditionalLikelihood via pruning algorithm
    for (int n=0; n<treePtr->getNumNodes(); n++)
    {
        Node* p = treePtr->getDownPassNode(n);


        if ( p->isTip() == false )
        {
            double* clP = p->getConditionalLikelihood()->getCls();
            std::vector<double*> clD;
            std::vector<double**> tiD;
            for (int i=0; i<p->getNumDescendants(); i++)
            {
                double* x = p->getDescendantIndexed(i)->getConditionalLikelihood()->getCls();
                clD.push_back(x);
                double** y = p->getDescendantIndexed(i)->getTransitionProbability()->getTiProb();
                tiD.push_back(y);
            }
            
            // for each Area for a given Node
            for (int a=0; a<areasPtr->getNumAreas(); a++)
            {
                for (int i=0; i<2; i++)
                {
                    double prob = 1.0;

                    // computes the cond. like. given the descendants of the Node (e.g. 2 for bifurcating tree)
                    for (int d=0; d<p->getNumDescendants(); d++)
                    {
                        double linProb = 0.0;
                        for (int j=0; j<2; j++)
                        {
                            linProb += clD[d][j] * tiD[d][i][j];
                        }
                        prob *= linProb;
                    }
                    clP[i] = prob;
                }

                // update pointers for the next Area
                clP += 2;
                for (int d=0; d<p->getNumDescendants(); d++)
                    clD[d] += 2;
                }
            }
        }

    //showConditionalLikelihoods();
    
    // step 2: sample a state for each interior node
    double stationaryFreqs[2];
    double pathRatesSum = pathRates[0] + pathRates[1];
    stationaryFreqs[0] = pathRates[0] / pathRatesSum; // 1->0
    stationaryFreqs[1] = pathRates[1] / pathRatesSum; // 0->1


    if (modelType == CARRYING_CAPACITY || modelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
    {
    	// I think the correct stationary frequencies will require solving the equilibrium values for a BD process w/ carrying capacity...
        // or approximating via simulation
    	stationaryFreqs[0] = 1.0 - carryingCapacity;
    	stationaryFreqs[1] = carryingCapacity;
    }

    int minumumNumChanges = 0;
    for (int n=treePtr->getNumNodes()-1; n>=0; n--)
    {
        Node* p = treePtr->getDownPassNode(n);
        if ( p->isTip() == false )
        {

            // used to ensure initial sampled cfgs exclude the extinction cfg
            bool hasOne = false;
            double bestClP = 0.0;
            double* bestClPPtr = NULL;


            //std::cout << "node " << p->getIndex() << "\n";
            if (p->getAncestor() == NULL)
            {
                // at the root
                double* clP = p->getConditionalLikelihood()->getCls();
                for (int a=0; a<areasPtr->getNumAreas(); a++)
                {

                    double sum = 0.0;
                    for (int i=0; i<2; i++)
                    {
                        sum += clP[i];
                    }
                    if (clP[1] / sum > bestClP)
                    {
                        bestClPPtr = clP;
                        bestClP = clP[1] / sum;
                    }

                    double u = ranPtr->uniformRv();
                    if (u < clP[0] / sum)//stationaryFreqs[0])// * clP[0] / sum)
                    {
                        clP[0] = 1.0;
                        clP[1] = 0.0;
                    }
                    else
                    {
                        clP[0] = 0.0;
                        clP[1] = 1.0;
                        hasOne = true;
                    }
                    clP += 2;
                }
            }
            else
            {
                // interior node (besides root)
                double* clP = p->getConditionalLikelihood()->getCls();
                double* clA = p->getAncestor()->getConditionalLikelihood()->getCls();
//                double** tiP = p->getTransitionProbability()->getTiProb();
                for (int a=0; a<areasPtr->getNumAreas(); a++)
                {
                    int ancState = 0;
                    if (clA[1] > 0.0)
                        ancState = 1;

                    double sum = 0.0;
                    for (int i=0; i<2; i++) {
                        sum += clP[i];
                    }
                    double u = ranPtr->uniformRv();
                    int desState;
                  //  std::cout << u << " < " << clP[0] << " / " << sum << " = " << clP[0] / sum << "\n";

                    if (clP[1] / sum > bestClP)
                    {
                        bestClPPtr = clP;
                        bestClP = clP[1] / sum;
                    }

                    if (u < clP[0] / sum)
                    {
                        desState = 0;
                        clP[0] = 1.0;
                        clP[1] = 0.0;
                    }
                    else
                    {
                        desState = 1;
                        clP[0] = 0.0;
                        clP[1] = 1.0;
                        hasOne = true;
                    }
                    if (ancState != desState)
                        minumumNumChanges++;
                    clP += 2;
                    clA += 2;
                }
            }

            // if a node is initialized to an extinction cfg, flip the area with the highest clP to on
            if (hasOne == false)
            {
                //std::cout << "resurrect " << p->getIndex() << "\n";
                bestClPPtr[0] = 0.0;
                bestClPPtr[1] = 1.0;
            }
		}
	}
    //std::cout << "Minimum number of changes = " << minumumNumChanges << std::endl;
    //showConditionalLikelihoods();
    
    // step 3: simulate histories for each branch
    std::set<int> areaSet;
    for (int i = 0; i < areasPtr->getNumAreas(); i++)
        areaSet.insert(i);

    for (int n=0; n<treePtr->getNumNodes(); n++)
        {
        Node* p = treePtr->getDownPassNode(n);
        //std::cout << "MODEL: Initializing node " << p->getIndex() << "\n";

        // non-root
        if (p->getAncestor() != NULL)
        {
            bool resampleNode = true;

            double* clA = p->getAncestor()->getConditionalLikelihood()->getCls();
            double* clP = p->getConditionalLikelihood()->getCls();
            for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
            {
                int startState = 0, endState = 0;
                if (clA[1] > 0.0)
                {
                    startState = 1;
                    resampleNode = false;
                }
                if (clP[1] > 0.0)
                    endState = 1;
                clA += 2;
                clP += 2;

                p->getBranchHistory(1)->setAncestralArea(aIdx, startState);
            }
            while(drawChanges(p, areaSet)) p->getBranchHistory(1)->removeChanges(areaSet);
            //drawChanges(p, areaSet);

            (*p->getBranchHistory(0)) = (*p->getBranchHistory(1));
        }

        // root
        else
        {

            std::vector<bool> consensusAncestor;
            int nOn = 0;
            do 
            {
            	// consensus ancestral area is occupied w.p. of tip occupancy frequency
                consensusAncestor.clear();
                for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
                {
                    int nOccuppied = areasPtr->getNumOccuppiedInArea(aIdx);
                    double prob = (double)nOccuppied / areasPtr->getNumTaxa();
                    if (ranPtr->uniformRv() < prob)
                        consensusAncestor.push_back( true );
                    else
                        consensusAncestor.push_back( false );
                }
                nOn = 0;
                for (int i=0; i<consensusAncestor.size(); i++)
                {
                    if (consensusAncestor[i] == true)
                        nOn++;
                }
            } while (nOn == 0);
            
            for (int aIdx=0; aIdx<areasPtr->getNumAreas(); aIdx++)
            {
                int startState = 0;
//                int endState = 0;
                if (consensusAncestor[aIdx] == true)
                    startState = 1;
                //if ( p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralStateForArea(aIdx) == true )
                //    endState = 1;
                
                p->getBranchHistory(1)->setAncestralArea(aIdx, startState);
                
            }
            drawChanges(p, areaSet);


            /*
            std::vector<bool> s;
            s = p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralState();
            for (int i = 0; i < s.size(); i++)
            	std::cout << s[i];
            std::cout << "\n";
            */

            (*p->getBranchHistory(0)) = (*p->getBranchHistory(1));
        }
    }
    
    // copy histories from 0 -> 1
    for (int n=0; n<treePtr->getNumNodes(); n++)
    {
        Node* p = treePtr->getDownPassNode(n);
        BranchHistory* bh0 = p->getBranchHistory(0);
        BranchHistory* bh1 = p->getBranchHistory(1);
        (*bh1) = (*bh0);
    }
}

double Model::lnFactorial(int n) {
    
    double lnF = 0.0;
    for (int i=1; i<=n; i++)
        lnF += log((double)i);
    return lnF;
}

double Model::lnBinomialCoefficient(int n, int k) {
    
    double lnN = lnFactorial(n) - lnFactorial(k) - lnFactorial(n-k);
    return lnN;
}

double Model::lnProbK(int n, int k) {

    double lnP = lnFactorial(n) - lnFactorial(k) - lnFactorial(n-k) - (double)n * log(2.0);
    return lnP;
}

double Model::lnLikelihood(int space) {
    
	//return 0.0;

    double lnL = 0.0;
#if DEBUG_LIKELIHOOD2
    std::cout << "lnLikelihoodForTree\n";
#endif

    for (int n=0; n<treePtr->getNumNodes(); n++)
    {
        Node* p = treePtr->getDownPassNode(n);
#if DEBUG_LIKELIHOOD
        std::stringstream debugSs;
#endif
        double v = 0.0;
//        double dLnL = lnL;
        if (p->getAncestor() != NULL)
        {
            //std::cout << "\tNode " << p->getIndex() << "\n";
            v = p->getLen();
            BranchHistory* h = p->getBranchHistory(space);
            std::set<AreaChange*,comp_history> changes = h->getBranchChanges();
            std::vector<bool> fromState = h->getAncestralState();
//            double sumTransitionRate = 0.0;
            double waitingTimeSum = 0.0;
            double pos = 0.0;
            double eventTime = 0.0; // for dynamic maps, currently disabled
            double lambda = sumOfRates(modelType, eventTime, areaRates, fromState); // sum of rates away from fromState
            for (std::set<AreaChange*,comp_history>::iterator it = changes.begin(); it != changes.end(); it++)
            {

                // get the configuration for the "to" state
                std::vector<bool> toState = fromState;
                toState[(*it)->getArea()].flip();

                // how long have we waited for the change?
                double waitingTime = ((*it)->getPosition() - pos) * v;
                waitingTimeSum += waitingTime;

                // get the rate of transition from fromState to toState
                double rateFromTo = transitionRate(modelType, eventTime, areaRates, fromState, toState, (*it)->getArea());

                // deal with the factors added by this change in area
                lnL += -(lambda * waitingTime) + log(rateFromTo);

	#if DEBUG_LIKELIHOOD
                debugSs << "\t\t\t" << p->getIndex() << "," << (*it)->getArea() << ":(" << fromState[(*it)->getArea()] << "->" << toState[(*it)->getArea()];
	#endif
                pos = (*it)->getPosition();
                fromState = toState;
                lambda = sumOfRates(modelType, eventTime, areaRates, fromState);

	#if DEBUG_LIKELIHOOD

               debugSs << "\t" << waitingTime << ")\n";
               debugSs << "\t\t\t\tt:   " << waitingTime << "\tr:\t" << rateFromTo << "\n";
               debugSs << "\t\t\t\tnumOn: " << numOnStates(fromState) << " / " << numAreas << "\n";
               debugSs << "\t\t\t\tl_10:\t" << areaRates[0] << "\tl_01:\t" << areaRates[1] << "\tl:\t" << lambda << "\n";

               debugSs << "\t\t\t\tlnP: " << lnL << " += " << log(rateFromTo) << " - " << lambda << " * " << waitingTime << "\n";
            //   debugSs << "\t\t\t\tlnLinc: " << lnL << "\trft " << rateFromTo << "\n";
               debugSs << "\t\t\t\twaitingTimeSum:\t" << waitingTimeSum << " / " << v << "\n";
               std::cout << debugSs.str();
               debugSs.clear();
	#endif
            }
            // one more factor for this branch that is the probability of no change along a branch of this length
            lnL += -lambda * ( (1.0-pos)*v );
        }
#if DEBUG_LIKELIHOOD2
        std::cout << "\t\t\tlnProb?\t" << p->getIndex() << "\t" << dLnL - lnL << "\n";
#endif
    }

    return lnL;
}


double Model::lnLikelihoodForPath(Node* p, int space) {

	//return 0.0;

    double lnL = 0.0;

#if DEBUG_LIKELIHOOD2
        std::stringstream debugSs;
        std::cout << "lnLikelihoodForPath\n";
#endif
    double v = 0.0;
    if (p->getAncestor() != NULL)
    {
        v = p->getLen();
        BranchHistory* h = p->getBranchHistory(space);
        std::set<AreaChange*,comp_history> changes = h->getBranchChanges();
        std::vector<bool> fromState = h->getAncestralState();
//        double sumTransitionRate = 0.0;
        double waitingTimeSum = 0.0;
        double pos = 0.0;
        double eventTime = 0.0; // for dynamic maps, currently disabled
        double lambda = sumOfRates(modelType, eventTime, areaRates, fromState); // sum of rates away from fromState
        for (std::set<AreaChange*,comp_history>::iterator it = changes.begin(); it != changes.end(); it++)
        {

            // get the configuration for the "to" state
            std::vector<bool> toState = fromState;
            toState[(*it)->getArea()].flip();

            // how long have we waited for the change?
            double waitingTime = ((*it)->getPosition() - pos) * v;
            waitingTimeSum += waitingTime;

            // get the rate of transition from fromState to toState
            double rateFromTo = transitionRate(modelType, eventTime, areaRates, fromState, toState, (*it)->getArea());

            // deal with the factors added by this change in area
            lnL += -(lambda * waitingTime) + log(rateFromTo);

#if DEBUG_LIKELIHOOD
            debugSs << "\t\t\t" << p->getIndex() << "," << (*it)->getArea() << ":(" << fromState[(*it)->getArea()] << "->" << toState[(*it)->getArea()];
#endif
            pos = (*it)->getPosition();
            fromState = toState;
            lambda = sumOfRates(modelType, eventTime, areaRates, fromState);

#if DEBUG_LIKELIHOOD

           debugSs << "\t" << waitingTime << ")\n";
           debugSs << "\t\t\t\tt:   " << waitingTime << " = (" << (*it)->getPosition() << " - " << pos << ") * " << v << "\n";
           debugSs << "\t\t\t\tl10:\t" << areaRates[0] << "\tl_01:\t" << areaRates[1] << "\tl:\t" << lambda << "\n";

           debugSs << "\t\t\t\tlnP: " << lnL << " += " << log(rateFromTo) << " - " << lambda << " * " << waitingTime << "\n";
        //   debugSs << "\t\t\t\tlnLinc: " << lnL << "\trft " << rateFromTo << "\n";
           debugSs << "\t\t\t\twaitingTimeSum:\t" << waitingTimeSum << " / " << v << "\n";
#endif
        }
        // one more factor for this branch that is the probability of no change along a branch of this length
        lnL += -lambda * ( (1.0-pos)*v );
    }
#if DEBUG_LIKELIHOOD2
        std::cout << debugSs.str() << "\n";
        std::cout << "\t\t\tlnProb!\t" << p->getIndex() << "\t" << lnL << "\n\n";
#endif

    return lnL;

}


double Model::lnProposedStateProb(int space) {

    //return 0.0;

#if DEBUG_PROPOSE_BRANCH
    std::cout << "\t\tBranch proposal\n";
#endif
    double lnProb = 0.0;
    for (int n=0; n<treePtr->getNumNodes(); n++)
    {
    	Node* p = treePtr->getDownPassNode(n);

        BranchHistory* h = p->getBranchHistory(space);
        double v = 0.0;
		std::stringstream debugSs;
        if (p->getAncestor() != NULL)
        {
        //	if (p->getAncestor()->getAncestor() != NULL)
        //	{
				v = p->getLen(); // absolute time
				std::set<AreaChange*,comp_history> changes = h->getBranchChanges();
				std::vector<bool> fromState = h->getAncestralState();
				double pos = 0.0;
				double eventTime = 0.0;
				double lambda = sumOfRates(INDEPENDENCE, 0.0, pathRates, fromState);
				double waitingTimeSum = 0.0;
		#if DEBUG_PROPOSE_BRANCH

				debugSs << "\t\tn" << p->getIndex() << "\tv" << v << "\n";
		#endif
				for (std::set<AreaChange*,comp_history>::iterator it = changes.begin(); it != changes.end(); it++)
				{

					// get the configuration for the "to" state
					std::vector<bool> toState = fromState;
					toState[(*it)->getArea()].flip();

					// how long have we waited for the change?
					double waitingTime = ((*it)->getPosition() - pos) * v;
					waitingTimeSum += waitingTime;

					// when did the event occur in absolute time? (this informs which map to load)
					// eventTime += waitingTime; // currently disabled

					// get the rate of transition from the "from" state to the "to" state
					double rateFromTo = transitionRate(INDEPENDENCE, eventTime, pathRates, fromState, toState, (*it)->getArea());

					// deal with the factors added by this change in area
					lnProb += log(rateFromTo) - (waitingTime * lambda);

	#if DEBUG_PROPOSE_BRANCH
				debugSs << "\t\t\t" << (*it)->getArea() << ":(" << fromState[(*it)->getArea()] << "->" << toState[(*it)->getArea()];
				debugSs << "\t" << waitingTime << ")\n";
				debugSs << "\t\t\t\tt:   " << waitingTime << " = (" << (*it)->getPosition() << " - " << pos << ") * " << v << "\n";
				debugSs << "\t\t\t\tl10:\t" << pathRates[0] << "\tl_01:\t" << pathRates[1] << "\tl:\t" << lambda << "\tr:\t" << rateFromTo << "\n";
				debugSs << "\t\t\t\tlnP: " << lnProb << " += " << log(rateFromTo) << " - " << lambda << " * " << waitingTime << "\n";
				debugSs << "\t\t\t\twaitingTimeSum:\t" << waitingTimeSum << " / " << v << "\n";
	#endif

					pos = (*it)->getPosition();
					fromState = toState;
					lambda = sumOfRates(INDEPENDENCE, eventTime, pathRates, fromState);
				}
				// one more factor for this branch that is the probability of no change along a branch of this length
				lnProb += -lambda * ( (1.0-pos)*v );
        //    }
#if DEBUG_PROPOSE_BRANCH
        std::cout << debugSs.str() << "\n";
        std::cout << "\t\t\tlnProb\t" << lnProb << "\n";
#endif
        }
	}

    // prior probability of state at the root of the tree
    lnProb += 0.0; // probability of observing root is 1
    return lnProb;
}

double Model::lnProposedStateProbForPath(Node* p, int space) {

    //return 0.0;

#if DEBUG_PROPOSE_BRANCH
    std::cout << "\t\tBranch proposal\n";
#endif
    double lnProb = 0.0;

	BranchHistory* h = p->getBranchHistory(space);
	double v = 0.0;
	std::stringstream debugSs;
	if (p->getAncestor() != NULL)
	{
		//if (p->getAncestor()->getAncestor() != NULL)

	//	{
			v = p->getLen(); // absolute time
			std::set<AreaChange*,comp_history> changes = h->getBranchChanges();
			std::vector<bool> fromState = h->getAncestralState();
			double pos = 0.0;
			double eventTime = 0.0;
			double lambda = sumOfRates(INDEPENDENCE, 0.0, pathRates, fromState);
			double waitingTimeSum = 0.0;
	#if DEBUG_PROPOSE_BRANCH

			debugSs << "\t\tn" << p->getIndex() << "\tv" << v << "\n";
	#endif
			for (std::set<AreaChange*,comp_history>::iterator it = changes.begin(); it != changes.end(); it++)
			{

				// get the configuration for the "to" state
				std::vector<bool> toState = fromState;
				toState[(*it)->getArea()].flip();

				// how long have we waited for the change?
				double waitingTime = ((*it)->getPosition() - pos) * v;
				waitingTimeSum += waitingTime;

				// when did the event occur in absolute time? (this informs which map to load)
				// eventTime += waitingTime; // currently disabled

				// get the rate of transition from the "from" state to the "to" state
				double rateFromTo = transitionRate(INDEPENDENCE, eventTime, pathRates, fromState, toState, (*it)->getArea());

				// deal with the factors added by this change in area
				lnProb += log(rateFromTo) - (waitingTime * lambda);

		#if DEBUG_PROPOSE_BRANCH
				debugSs << "\t\t\t" << (*it)->getArea() << ":(" << fromState[(*it)->getArea()] << "->" << toState[(*it)->getArea()];
				debugSs << "\t" << waitingTime << ")\n";
				debugSs << "\t\t\t\tt:   " << waitingTime << " = (" << (*it)->getPosition() << " - " << pos << ") * " << v << "\n";
				debugSs << "\t\t\t\tl10:\t" << pathRates[0] << "\tl_01:\t" << pathRates[1] << "\tl:\t" << lambda << "\tr:\t" << rateFromTo << "\n";
				debugSs << "\t\t\t\tlnP: " << lnProb << " += " << log(rateFromTo) << " - " << lambda << " * " << waitingTime << "\n";
				debugSs << "\t\t\t\twaitingTimeSum:\t" << waitingTimeSum << " / " << v << "\n";
		#endif

				pos = (*it)->getPosition();
				fromState = toState;
				lambda = sumOfRates(INDEPENDENCE, eventTime, pathRates, fromState);
			}
			// one more factor for this branch that is the probability of no change along a branch of this length
			lnProb += -lambda * ( (1.0-pos)*v );

		#if DEBUG_PROPOSE_BRANCH
			std::cout << debugSs.str() << "\n";
			std::cout << "\t\t\tlnProb\t" << lnProb << "\n";
		#endif

	//	}
	}

	// prior probability of state at the root of the tree
	lnProb += 0.0; // probability of observing root is 1
	return lnProb;
}

int Model::numOnStates(std::vector<bool>& x) {

    int n = 0;
    for (int i=0; i<x.size(); i++)
        {
        if (x[i] == true)
            n++;
        }
    return n;
}

void Model::showConditionalLikelihoods(void) {

    for (int i=0; i<treePtr->getNumNodes(); i++)
        {
        Node* p = treePtr->getDownPassNode(i);
        std::cout << std::setw(4) << p->getIndex() << " -- ";
        ConditionalLikelihood* c = p->getConditionalLikelihood();
        c->print();
        std::cout << std::endl;
        }
}

void Model::showTransitionProbabilities(void) {

    for (int i=0; i<treePtr->getNumNodes(); i++)
        {
        Node* p = treePtr->getDownPassNode(i);
        std::cout << std::setw(4) << p->getIndex() << ":" << std::endl;
        TransitionProbability* tp = p->getTransitionProbability();
        tp->print();
        std::cout << std::endl;
        }
}

double Model::sumOfRates(int srModelType, double eventTime, std::vector<double>& rates, std::vector<bool>& x) {

	// returns the sum of rates away from the current state
	double r = 0.0;
	int numOn = numOnStates(x);
	int numOff = (int)x.size() - numOn;

	if (srModelType == INDEPENDENCE || srModelType == DISTANCE_NORM)
	{
		// area gain
		r = numOff * rates[1];

		// area loss
		if ( numOn != 1 )
			r += numOn * rates[0];
	}
	else if (srModelType == CARRYING_CAPACITY || srModelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
	{
		// area gain
		r = numOff * rates[1] * pow(numAreas * carryingCapacity / (double)numOn, carryingPower);
		//std::cout << r << "\t" << carryingCapacity << "\t" << numOn << "\t" << carryingPower << "\n";
		// area loss
		if (numOn != 1)
			r += numOn * rates[0] * pow((double)numOn / (numAreas * carryingCapacity), carryingPower);
		//std::cout << r << "\t" << carryingCapacity << "\t" << numOn << "\t" << carryingPower << "\n";
	}

    return r;
}

double Model::transitionRate(int trModelType, double eventTime, std::vector<double>& rates, std::vector<bool>& fromState, std::vector<bool>& toState, int toArea) {
    
	double transitionRate = 0.0;
	int numOn = numOnStates(fromState);
	int numOff = numAreas - numOn;

	if (trModelType == CARRYING_CAPACITY)
	{

		// area loss
		if (numOn - numOnStates(toState) > 0)
			transitionRate = rates[0] * pow((double)numOn / (numAreas * carryingCapacity), carryingPower);

		// area gain
		else
		    transitionRate = rates[1] * pow(numAreas * carryingCapacity / (double)numOn, carryingPower);
	}

	else if (trModelType == INDEPENDENCE)
	{
		// area loss
		if (numOn - numOnStates(toState) > 0)
			transitionRate = rates[0];

		// area gain
		else
			transitionRate = rates[1];
	}

	else if (trModelType == DISTANCE_NORM || trModelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
	{

		// area loss
		if (numOn - numOnStates(toState) > 0) {

			transitionRate = rates[0];

			// add carrying capacity effects
			if (modelType == CARRYING_CAPACITY_AND_DISTANCE_NORM) {
				transitionRate *= pow((double)numOn / (numAreas * carryingCapacity), carryingPower);

			}
		}


		// area gain
		else
		{
			double toAreaImmigrationRate = 0.0;
			double sumImmigrationRate = 0.0;
			double meanImmigrationRate = 0.0;

			// generate vector of indices for on/off areas
			std::vector<int> onIdx;
			std::vector<int> offIdx;
			for (int i = 0; i < numAreas; i++)
			{
				if (fromState[i] == 1)
					onIdx.push_back(i);
				else
					offIdx.push_back(i);
			}

			// calculate sum of distances-weights to toArea and to all areas
			std::vector<int>::iterator onIdxIt;
			std::vector<int>::iterator offIdxIt;
			for (onIdxIt = onIdx.begin(); onIdxIt != onIdx.end(); onIdxIt++)
			{
			    //std::cout << *onIdxIt << "\n";
				// for toArea
				toAreaImmigrationRate += distances[*onIdxIt][toArea];

				// for all areas
				for (offIdxIt = offIdx.begin(); offIdxIt != offIdx.end(); offIdxIt++)
				{
				   // std::cout << "\t" << *offIdxIt << "\n";
					sumImmigrationRate += distances[*onIdxIt][*offIdxIt];
				}
			}
			meanImmigrationRate = sumImmigrationRate / numOff;

			transitionRate = rates[1] * toAreaImmigrationRate / meanImmigrationRate;
			//std::cout << "\t\tTR1\t" << transitionRate << "\n";
			if (modelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
				transitionRate *= pow(numAreas * carryingCapacity / (double)numOn, carryingPower);
			//std::cout << "\t\tTR2\t" << transitionRate << "\n";

#if DEBUG_DISTANCE_TR_RATE
			double sumTransitionRate = rates[1] * numOff;
			std::cout << "\t\tmeanIR:\t" << meanImmigrationRate << "\t=\t" << sumImmigrationRate << "\t/\t" << numOff << "\n";
			std::cout << "\t\tratioR:\t" << rates[1] / rates[0] << "\t=\t" << rates[1] << "\t/\t" << rates[0] << "\n";
			std::cout << "\t\t    TR:\t" << transitionRate << "\t=\t" << rates[1] << "\t*\t" << toAreaImmigrationRate << "\t/\t" << meanImmigrationRate << "\n";
			std::cout << "\t\t sumTR:\t" << sumTransitionRate << "\t+\t" << numOn * rates[1] << "\n";
			std::cout << "\n";
#endif
		}
	}
	//std::cout << "\t\tTR\t" << transitionRate << "\n";
	if (std::isnan(transitionRate)) {
		transitionRate = 10e-300;
	}
	return transitionRate;
}

void Model::updateNode(int mcmcCycle, int na) {
    
    // pick a random node to update
    Node* p = treePtr->randomInteriorNode(ranPtr);
    //std::cout << "update node:\t" << p->getIndex() << "\n";

    // pick areas to update
    std::set<int> areaSet;
    while ( areaSet.size() < na )
    {
        //std::cout << areaSet.size() << " " << na << "\n";
        areaSet.insert( (int)(ranPtr->uniformRv()*areasPtr->getNumAreas()) );
    }

    //std::cout << "\tlnProp\n";
    // compute numerator of lnProposalRatio paths incident to node
    double lnProposalRatio = 0.0;
    if (p->getAncestor() != NULL) {
    	lnProposalRatio += lnProposedStateProbForPath(p, 0);
    }
  	for (int i = 0; i < p->getDescendants().size(); i++) {
  		Node* q = p->getDescendants()[i];
  		lnProposalRatio += lnProposedStateProbForPath(q, 0);
  	}

#   if DEBUG_PROPOSE_NODE_PATH
    std::cout << "Areas to update: ";
    for (std::set<int>::iterator it = areaSet.begin(); it != areaSet.end(); it++)
        std::cout << (*it) << " ";
    std::cout << std::endl;
    std::cout << "Node to update: " << p->getIndex() << std::endl;
#   endif
    
    //std::cout << "\tremChanges\n";
    // get a reference to the descendant nodes
    std::vector<Node*> descendants = p->getDescendants();
        
    // erase the changes on the branches incident to p
    BranchHistory* h = p->getBranchHistory(1);
    h->removeChanges(areaSet);

    for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
    {
        BranchHistory* h = (*it)->getBranchHistory(1);
        h->removeChanges(areaSet);
    }

#   if DEBUG_PROPOSE_NODE_PATH
    if (p->getAncestor() != NULL)
        {
        BranchHistory* h0 = p->getBranchHistory(0);
        BranchHistory* h1 = p->getBranchHistory(1);
        std::cout << "History along root branch before changes removed: " << std::endl;
        h0->print();
        std::cout << "History along root branch after changes removed: " << std::endl;
        h1->print();
        }
    for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        BranchHistory* h0 = (*it)->getBranchHistory(0);
        BranchHistory* h1 = (*it)->getBranchHistory(1);
        std::cout << "History along root branch before changes removed: " << std::endl;
        h0->print();
        std::cout << "History along root branch after changes removed: " << std::endl;
        h1->print();
        }
#   endif
    

    // make new ancestral areas for node p
    //std::cout << "\tdrawareas\n";
//    bool valid_node_sample = drawAncestralAreas(p, areaSet);
    
    int num_attempts = 1e2;
    bool valid_node_sample = false;
    while (num_attempts >= 0 && valid_node_sample == false) {
        valid_node_sample = drawAncestralAreas(p, areaSet);
        num_attempts -= 1;
    }
    

    if (valid_node_sample) {
        // update histories for each branch incident to node p
        //std::cout << "\taddChanges\n";
        while(drawChanges(p, areaSet)) h->removeChanges(areaSet);
        for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
            while(drawChanges(*it, areaSet)) (*it)->getBranchHistory(1)->removeChanges(areaSet);
    }
    
/*
 * drawChanges(p, areaSet);
    for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        drawChanges(*it, areaSet);
 *
 */

#   if DEBUG_PROPOSE_NODE_PATH
    if (p->getAncestor() != NULL)
        {
        BranchHistory* h0 = p->getBranchHistory(0);
        BranchHistory* h1 = p->getBranchHistory(1);
        std::cout << "History along root branch before: " << std::endl;
        h0->print();
        std::cout << "History along root branch after: " << std::endl;
        h1->print();
        }
    for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        BranchHistory* h0 = (*it)->getBranchHistory(0);
        BranchHistory* h1 = (*it)->getBranchHistory(1);
        std::cout << "History along root branch before: " << std::endl;
        h0->print();
        std::cout << "History along root branch after: " << std::endl;
        h1->print();
        }
#   endif
    

    // calculate the acceptance probability
    //std::cout << "\tacceptanceProb\n";
    double newLnL = lnLike;
    if (p->getAncestor() != NULL)
    	newLnL += lnLikelihoodForPath(p, 1) - lnLikelihoodForPath(p, 0);
    for (int i = 0; i < p->getDescendants().size(); i++) {
    	Node* q = p->getDescendants()[i];
    	newLnL += lnLikelihoodForPath(q, 1) - lnLikelihoodForPath(q, 0);
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " -- " << std::fixed << std::setprecision(8) << lnLike << " -> " << newLnL;
    double lnLikelihoodRatio = newLnL - lnLike;
    if (useSteppingStone == true)
    	lnLikelihoodRatio *= betaSteppingStone;

    if (p->getAncestor() != NULL)
    	lnProposalRatio -= lnProposedStateProbForPath(p, 1);
	for (int i = 0; i < p->getDescendants().size(); i++) {
		Node* q = p->getDescendants()[i];
		lnProposalRatio -= lnProposedStateProbForPath(q, 1);
	}

	//lnProposalRatio = calcLnProposalRatio();

    double lnR = lnLikelihoodRatio + lnProposalRatio;
    // Reject samples that cannot sample valid histories after a reasonable number of attempts.
    // This is valid because we could, in theory, have taken the first invalid proposed state
    // as valid, then let the likelihood function reject the invalid sample (prob=0).
    if (!valid_node_sample) {
        lnR = -500;
    }
#if DEBUG_PROPOSE_NODE
 //   double slowLnProposalRatio = lnProposedStateProb(0) - lnProposedStateProb(1);
 //   double slowLnL0 = lnLikelihood(0);
 //   double slowLnL1 = lnLikelihood(1);

    std::cout << "\n";
    std::cout << "\t\tnewLnL\t" << newLnL << "\toldLnL\t" << lnLike << "\n";//<< "\tslowLnL1\t" << slowLnL1  << "\tslowLnL0\t" << slowLnL0 << "\n";
   // std::cout << "\t\tslowLnProposalRatio\t" << slowLnProposalRatio << "\n";
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
   // std::cout << "\t\tnumChgCrwnGrp\t" << numChangesInCrownGroup << "\n";
    std::cout << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
    	r = exp(lnR);

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;

    /*
    if (newLnL != newLnL) std::cout << "\t\t\t\tnewLnL == nan\n";
    if (lnLike != lnLike) std::cout << "\t\t\t\toldLnL == nan\n";
    if (newLnL < -pow(10,300)) std::cout << "\t\t\t\tnewLnL == -Inf\n";
    if (lnLike < -pow(10,300)) std::cout << "\t\t\t\toldLnL == -Inf\n";
    */
    // update the state of the chain
    std::string result;
    if (accept == true)
    {
    	(*p->getBranchHistory(0)) = (*p->getBranchHistory(1));
        for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
            (*(*it)->getBranchHistory(0)) = (*(*it)->getBranchHistory(1));
        lnLike = newLnL;
        result = "Accept new history around node #";
    }
    else
    {
        (*p->getBranchHistory(1)) = (*p->getBranchHistory(0));
        for (std::vector<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
            (*(*it)->getBranchHistory(1)) = (*(*it)->getBranchHistory(0));
        result = "Reject new history around node #";
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
    {
        std::cout << " " << lnR << " " << result << p->getIndex() << " (" << getNumGainsTree(0) << "," << getNumLossesTree(0) << "," << treePtr->getRoot()->getBranchHistory(0)->getNumChanges() << "," <<  getNumGainsTree(0) + getNumLossesTree(0) << ")" << std::endl;
        //treePtr->getRoot()->getBranchHistory(0)->print();
    }
}

void Model::updateAreaRate(int mcmcCycle, double tuning) {

	int m = (ranPtr->uniformRv() < 0.5 ? 1 : 0);
    double oldRate = areaRates[m];
    double newRate = oldRate * exp( -tuning*(ranPtr->uniformRv()-0.5) );
    //double newRate = oldRate + tuning*(ranPtr->uniformRv()-0.5);
    if (newRate < 0) newRate *= -1;
    
    areaRates[m] = newRate;

    double newLnL = lnLikelihood(0);
    double lnLikelihoodRatio = newLnL - lnLike;
    if (useSteppingStone == true)
        lnLikelihoodRatio *= betaSteppingStone;

    //double lnProposalRatio   = 0.0;
    double lnProposalRatio = log(newRate) - log(oldRate); // equiv. to log(u)
    //double lnPriorRatio      = ranPtr->lnExponentialPdf(areaRateParm[m], newRate) - ranPtr->lnExponentialPdf(areaRateParm[m], oldRate);
    double lnPriorRatio      = log( pow(areaRateParm[m], 2) + pow(oldRate, 2) ) - log( pow(areaRateParm[m], 2) + pow(newRate, 2) );
    double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

#if DEBUG_PROPOSE_RATE
    std::cout << "\t\tAreaRate\n";
    std::cout << "\t\trate -> rate\'\t" << oldRate << " -> " << newRate << "\n";
    std::cout << "\t\tnewLnL\t" << newLnL << "\toldLnL\t" << lnLike << "\n";
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
    std::cout << "\t\tlnPriorRatio\t" << lnPriorRatio << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
        r = exp(lnR);
    
    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " -- " << std::fixed << std::setprecision(8) << lnLike << " -> " << newLnL;

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;
    if (newRate > 1000)
        accept = false;

    std::string result;
    if (accept == true)
    {
        areaRates[m] = newRate;
        if (!settingsPtr->getUseAuxiliarySampling())
            pathRates[m] = newRate;
        lnLike = newLnL;
        result = "Accept new areaRates";
    }
    else
    {
        areaRates[m] = oldRate;
        if (!settingsPtr->getUseAuxiliarySampling())
            pathRates[m] = oldRate;
        result = "Reject new areaRates";
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << " " << lnR << " " << result << " (" << m << ": " << areaRates[m] << ")" << std::endl;
}

void Model::updatePathRate(int mcmcCycle, double tuning) {

	double oldPathLnL = lnProposedStateProb(0);

	int m = (ranPtr->uniformRv() < 0.5 ? 1 : 0);
    double oldRate = pathRates[m];
    double newRate = oldRate * exp( -tuning*(ranPtr->uniformRv()-0.5) );
    //double newRate = oldRate + tuning*(ranPtr->uniformRv()-0.5);
    //if (newRate < 0) newRate *= -1;

    pathRates[m] = newRate;

    double newPathLnL = lnProposedStateProb(0);
    double lnLikelihoodRatio = newPathLnL - oldPathLnL;
    if (useSteppingStone == true)
        lnLikelihoodRatio *= betaSteppingStone;

    double lnProposalRatio   = log(newRate) - log(oldRate);
    //double lnProposalRatio   = 0.0;
    //double lnPriorRatio      = ranPtr->lnExponentialPdf(areaRateParm[m], newRate) - ranPtr->lnExponentialPdf(areaRateParm[m], oldRate);
    double lnPriorRatio      = log( pow(pathRateParm[m], 2) + pow(oldRate, 2) ) - log( pow(pathRateParm[m], 2) + pow(newRate, 2) );
    double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

#if DEBUG_PROPOSE_RATE
    std::cout << "\t\tPathRate\n";
    std::cout << "\t\trate -> rate\'\t" << oldRate << " -> " << newRate << "\n";
    std::cout << "\t\tnewPathLnL\t" << newPathLnL << "\toldPathLnL\t" << lnLike << "\n";
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
    std::cout << "\t\tlnPriorRatio\t" << lnPriorRatio << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
        r = exp(lnR);

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " ~p " << std::fixed << std::setprecision(8) << oldPathLnL << " -> " << newPathLnL;

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;
    if (newRate > 1000)
        accept = false;

    std::string result;
    if (accept == true)
    {
        pathRates[m] = newRate;
        for (int n=0; n<treePtr->getNumNodes(); n++)
		{
			Node* p = treePtr->getDownPassNode(n);

            TransitionProbability* tp = p->getTransitionProbability();
            tp->setProbs( p->getLen(), pathRates[1], pathRates[0], numAreas );
		}

        result = "Accept new pathRates";
    }
    else
    {
        pathRates[m] = oldRate;
        result = "Reject new pathRates";
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << " " << lnR << " " << result << " (" << m << ": " << pathRates[m] << ")" << std::endl;
}


void Model::updateDistancePower(int mcmcCycle, double tuning) {


    double oldPower = distancePower;
    double lnProposalRatio = 0;
    
    if (settingsPtr->isGeoDistancePositive())
    {
        double scaler = exp(tuning * (ranPtr->uniformRv() - 0.5));
        distancePower *= scaler;
        lnProposalRatio = log(scaler);
    }
    else
    {
        distancePower += ranPtr->normalRv(0.0, tuning);
        lnProposalRatio = 0; // symmetry of normal pdf
    }

    //geoPtr->changeDistancePower("haversinePowerInverseTruncated", 0.0, distancePower);
    //distances = geoPtr->getDistance("haversinePowerInverseTruncated", 0.0);
    //geoPtr->changeDistancePower("haversinePowerInverse", 0.0, distancePower);
    //distances = geoPtr->getDistance("haversinePowerInverse", 0.0);
    geoPtr->changeDistancePower(geoDistancePowerStr, 0.0, distancePower);
    distances = geoPtr->getDistance(geoDistancePowerStr, 0.0);

    double newLnL = lnLikelihood(0);
    double lnLikelihoodRatio = newLnL - lnLike;
    if (useSteppingStone == true)
        lnLikelihoodRatio *= betaSteppingStone;

    double lnPriorRatio      = log( pow(distancePowerParm, 2) + pow(oldPower, 2) ) - log( pow(distancePowerParm, 2) + pow(distancePower, 2) );
    double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

#if DEBUG_PROPOSE_DISTANCE_POWER
    std::cout << "\t\trate -> rate\'\t" << oldPower << " -> " << distancePower << "\n";
    std::cout << "\t\tlnLike\t" << lnLike << "\t" << newLnL << "\n";
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
    std::cout << "\t\tlnPriorRatio\t" << lnPriorRatio << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
        r = exp(lnR);

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " -- " << std::fixed << std::setprecision(8) << lnLike << " -> " << newLnL;

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;

    std::string result;
    if (accept == true)
    {
        lnLike = newLnL;
        result = "Accept new distancePower";
    }
    else
    {
        //std::cout << "newPower\t" << distancePower << "\toldPower\t" << oldPower << "\n";
        distancePower = oldPower;
        geoPtr->changeDistancePower(geoDistancePowerStr, 0.0, distancePower);
        distances = geoPtr->getDistance(geoDistancePowerStr, 0.0);
        //geoPtr->changeDistancePower("haversinePowerInverseTruncated", 0.0, distancePower);
        //distances = geoPtr->getDistance("haversinePowerInverseTruncated", 0.0);
        //geoPtr->changeDistancePower("haversinePowerInverse", 0.0, distancePower);
        //distances = geoPtr->getDistance("haversinePowerInverse", 0.0);

        result = "Reject new distancePower";
        //std::cout << "\t\tpostLnL\t" << lnLikelihood(0) << "\n";
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << " " << lnR << " " << result << " (" << distancePower << ")" << std::endl;
}

void Model::updateCarryingPower(int mcmcCycle, double tuning) {

//    double oldPathLnL = lnProposedStateProb(0);

    double oldPower = carryingPower;
    carryingPower = oldPower * exp( -tuning*(ranPtr->uniformRv()-0.5) );

    double newLnL = lnLikelihood(0);
    double lnLikelihoodRatio = newLnL - lnLike;
    if (useSteppingStone == true)
        lnLikelihoodRatio *= betaSteppingStone;
    double lnProposalRatio   = log(carryingPower) - log(oldPower);
    double lnPriorRatio      = log( pow(carryingPowerParm, 2) + pow(oldPower, 2) ) - log( pow(carryingPowerParm, 2) + pow(carryingPower, 2) );
    double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

#if DEBUG_PROPOSE_CARRYING_POWER
    std::cout << "\t\trate -> rate\'\t" << oldPower << " -> " << carryingPower << "\n";
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
    std::cout << "\t\tlnPriorRatio\t" << lnPriorRatio << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
        r = exp(lnR);

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " -- " << std::fixed << std::setprecision(8) << lnLike << " -> " << newLnL;

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;

    std::string result;
    if (accept == true)
    {
        lnLike = newLnL;
        result = "Accept new carryingPower";
    }
    else
        {
        carryingPower = oldPower;
        result = "Reject new carryingPower";
        }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << " " << lnR << " " << result << " (" << carryingPower << ")" << std::endl;
}


void Model::updateCarryingCapacity(int mcmcCycle, double tuning) {

//    double oldPathLnL = lnProposedStateProb(0);

    double oldCapacity = carryingCapacity;
    carryingCapacity = ranPtr->betaRv(oldCapacity * tuning / (1 - oldCapacity), tuning);

    double newLnL = lnLikelihood(0);
    double lnLikelihoodRatio = newLnL - lnLike;
    if (useSteppingStone == true)
        lnLikelihoodRatio *= betaSteppingStone;

    // surprised I couldn't find a way to simplify these expressions (even w/ FullSimplify in Mathematica)
    double lnProposalRatio   = ranPtr->lnBetaPdf(carryingCapacity * tuning / (1 - carryingCapacity), tuning, oldCapacity)
                               - ranPtr->lnBetaPdf(oldCapacity * tuning / (1 - oldCapacity), tuning, carryingCapacity);
    double lnPriorRatio      = ranPtr->lnBetaPdf(carryingCapacityAlphaParm, carryingCapacityBetaParm, carryingCapacity)
                               - ranPtr->lnBetaPdf(carryingCapacityAlphaParm, carryingCapacityBetaParm, oldCapacity);

    double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;


#if DEBUG_PROPOSE_CARRYING_CAPACITY
    std::cout << "\t\trate -> rate\'\t" << oldCapacity << " -> " << carryingCapacity << "\n";
    //std::cout << "\t\t"
    std::cout << "\t\tlnLikeRatio\t" << lnLikelihoodRatio << "\n";
    std::cout << "\t\tlnProposalRatio\t" << lnProposalRatio << "\n";
    std::cout << "\t\tlnPriorRatio\t" << lnPriorRatio << "\n";
#endif

    double r = 0.0;
    if (lnR > 0.0)
        r = 1.0;
    else if (lnR < -300.0)
        r = 0.0;
    else
        r = exp(lnR);

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << std::setw(5) << mcmcCycle << " -- " << std::fixed << std::setprecision(8) << lnLike << " -> " << newLnL;

    // accept or reject the state
    bool accept = false;
    if (ranPtr->uniformRv() < r)
        accept = true;

    std::string result;
    if (accept == true)
    {
        //pathRates[m] = newRate;
        for (int n=0; n<treePtr->getNumNodes(); n++)
        {
            Node* p = treePtr->getDownPassNode(n);
            if (p->getAncestor() != NULL)
            {
                TransitionProbability* tp = p->getTransitionProbability();
                tp->setProbs( (p->getAncestor()->getTime()-p->getTime()), pathRates[1], pathRates[0], numAreas );
            }
        }

        lnLike = newLnL;
        result = "Accept new carryingCapacity";
    }
    else
    {
        carryingCapacity = oldCapacity;
        result = "Reject new carryingCapacity";
    }

    if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
        std::cout << " " << lnR << " " << result << " (" << carryingCapacity << ")" << std::endl;
}


// this is unneeded if the root branch is long enough,
// but may speed up computation by resampling the consensus below the root
void Model::updateConsensus(int mcmcCycle, int na) {

	// pick areas to update
	std::set<int> areaSet;
	while ( areaSet.size() < na )
		areaSet.insert( (int)(ranPtr->uniformRv()*areasPtr->getNumAreas()) );

	// propose new consensus configuration below root
	drawConsensus(areaSet);

	// MH
	double lnLikelihoodRatio = 0.0;
	double lnPriorRatio = 0.0;
	double oldLnP = lnProposedConsensusProb(areaSet, 0);
	double newLnP = lnProposedConsensusProb(areaSet, 1);
	double lnProposalRatio = 0.0;//oldLnP - newLnP;
	double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

	double r = 0.0;
	if (lnR > 0.0)
		r = 1.0;
	else if (lnR < -300.0)
		r = 0.0;
	else
		r = exp(lnR);

	if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
		std::cout << std::setw(5) << mcmcCycle << " ~c " << std::fixed << std::setprecision(8) << oldLnP << " -> " << newLnP;

	// accept or reject the state
	bool accept = false;
	if (ranPtr->uniformRv() < r)
		accept = true;

	std::string result;
	Node* p = treePtr->getRoot();
	if (accept == true)
	{
		(*p->getBranchHistory(0)) = (*p->getBranchHistory(1));
		result = "Accept new consensus";
	}
	else
	{
		(*p->getBranchHistory(1)) = (*p->getBranchHistory(0));
		result = "Reject new consensus";
	}

	if ( mcmcCycle % settingsPtr->getPrintFrequency() == 0 )
		std::cout << " " << lnR << " " << result << std::endl;

}
