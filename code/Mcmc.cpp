#include "Mcmc.h"


#define INDEPENDENCE 1
#define CARRYING_CAPACITY 2
#define DISTANCE_NORM 3
#define CARRYING_CAPACITY_AND_DISTANCE_NORM 4

Mcmc::Mcmc(Areas* a, GeoCoords* g, MbRandom* r, Model* m, Tree* t, Settings* s) {

    // remember the address of important objects
    ranPtr = r;
    modelPtr = m;
    treePtr = t;
    settingsPtr = s;
    areasPtr = a;
    geoPtr = g;
    numSamples = 0;
    
    // set the proposal probabilities
    proposalProbs.push_back(20.0);	// pathHistory
    proposalProbs.push_back(1.0);	// areaRates
    // pathRates
    if (settingsPtr->getUseAuxillarySampling())
        proposalProbs.push_back(1.0);
    else
        proposalProbs.push_back(0.0);
    proposalProbs.push_back(settingsPtr->getModelType() == DISTANCE_NORM ||
                            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM
                            ? 1.0 : 0.0);   // distancePower
    proposalProbs.push_back(settingsPtr->getModelType() == CARRYING_CAPACITY ||
                            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM
                            ? 1.0 : 0.0);   // carryingCapacity
    proposalProbs.push_back(settingsPtr->getModelType() == CARRYING_CAPACITY ||
                            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM
                            ? 1.0 : 0.0);   // carryingPower
    proposalProbs.push_back(0.0);	// consensus
    double sum = 0.0;
    for (int i=0; i<proposalProbs.size(); i++)
        sum += proposalProbs[i];
    for (int i=0; i<proposalProbs.size(); i++)
        proposalProbs[i] /= sum;
    for (int i=1; i<proposalProbs.size(); i++)
        proposalProbs[i] += proposalProbs[i-1];

    outStrm.open( settingsPtr->getOutputFileName().c_str(), std::ios::out );
    if (!outStrm) 
        Msg::error("Cannot open file \"" + settingsPtr->getOutputFileName() + "\"");

    histStrm.open( settingsPtr->getHistoryFileName().c_str(), std::ios::out );
	if (!histStrm)
		Msg::error("Cannot open file \"" + settingsPtr->getHistoryFileName() + "\"");

	areaPosteriorStrm.open( settingsPtr->getAreaPosteriorFileName().c_str(), std::ios::out );
	if (!areaPosteriorStrm)
	    Msg::error("Cannot open file \"" + settingsPtr->getAreaPosteriorFileName() + "\"");

	nhxStrm.open( settingsPtr->getNhxFileName().c_str(), std::ios::out );
    if (!nhxStrm)
        Msg::error("Cannot open file \"" + settingsPtr->getNhxFileName() + "\"");

    std::cout << settingsPtr->getOutputFileName() << "\n";
	// start MCMC
    runChain();
    
    outStrm.close();
    histStrm.close();
    areaPosteriorStrm.close();
    nhxStrm.close();

}

Mcmc::~Mcmc(void) {

}

void Mcmc::runChain(void) {

    std::cout << "MCMC initated\n\n";
    
    for (int n=1; n<=settingsPtr->getChainLength(); n++)
    {
        
        // propose a new state
        double u = ranPtr->uniformRv();

        if (u <= proposalProbs[0])
        {
            //std::cout << "MCMC: updateNode()\n";
            double numAreas = areasPtr->getNumAreas();
            int m = ranPtr->poissonRv(numAreas * settingsPtr->getAreaProposalTuner()) + 1;
            if (m > numAreas)
                m = numAreas;
            modelPtr->updateNode(n, m);
        }
        else if (u > proposalProbs[0] && u <= proposalProbs[1])
        {
            //std::cout << "MCMC: updateAreaRate()\n";
            modelPtr->updateAreaRate(n, settingsPtr->getRateProposalTuner());
        }
        else if (u > proposalProbs[1] && u <= proposalProbs[2])
        {
            //std::cout << "MCMC: updatePathRate()\n";
            modelPtr->updatePathRate(n, settingsPtr->getRateProposalTuner());
        }
        else if (u > proposalProbs[2] && u <= proposalProbs[3])
        {
            //std::cout << "MCMC: updateDistancePower()\n";
            modelPtr->updateDistancePower(n, settingsPtr->getDistanceProposalTuner());
        }
        else if (u > proposalProbs[3] && u <= proposalProbs[4])
        {
            //std::cout << "MCMC: updateCarryingCapacity()\n";
            modelPtr->updateCarryingCapacity(n, 1000.0);
        }
        else if (u > proposalProbs[4] && u <= proposalProbs[5])
        {
            //std::cout << "MCMC: updateCarryingPower()\n";
            modelPtr->updateCarryingPower(n, 0.25);
        }
        else if (u > proposalProbs[5] && u <= proposalProbs[6])
        {
            //std::cout << "MCMC: updateConsensus()\n";
            double t = settingsPtr->getAreaProposalTuner();
            int n = 0;
            do
            {
                n = ranPtr->poissonRv(areasPtr->getNumAreas() * t) + 1;
            } while(n > areasPtr->getNumAreas());

        	//modelPtr->updateConsensus(n, areasPtr->getNumAreas() * .01 + 1);
            modelPtr->updateConsensus(n, areasPtr->getNumAreas() * .05 + 1);
        }
        
        // update numSamples
        if (n >= settingsPtr->getProbBurnIn() &&
            n % settingsPtr->getHistorySampleFrequency() == 0)
            numSamples++;

        // sample areas on tree
        modelPtr->getTreePtr()->incrementAreaCounts();
        
        // sample parameter values
        if ( n % settingsPtr->getParameterSampleFrequency() == 0)
            sampleParms(n);

        // sample path values
        if ( n % settingsPtr->getHistorySampleFrequency() == 0)
            samplePaths(n);

        // update path counts (used to compute posteriors of area occupancy)
        if ( n % settingsPtr->getHistorySampleFrequency() == 0)
            updateAncestralStateCounts(n);
        
    }

    std::cout << "MCMC complete\n";
    std::cout << "Writing output files\n";
    writeAncestralStateFreqs();
    writeNhxFile();
    std::cout << "Analysis complete\n";

}

void Mcmc::sampleParms(int n) {

    std::cout << std::fixed << std::setprecision(5);
    if (n == settingsPtr->getParameterSampleFrequency())
    {
    	outStrm << "n\tlnL\tgain\tloss\tgain-p\tloss-p";
    	if (settingsPtr->getModelType() == DISTANCE_NORM ||
    	    settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM)
    	{
    	    outStrm << "\tdistP";
    	}
        if (settingsPtr->getModelType() == CARRYING_CAPACITY ||
            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM)
        {
            outStrm << "\tcarryC";
            outStrm << "\tcarryP";
        }
    	outStrm << "\tnumGain\tnumLoss" << std::endl;
    }

    if (n >= settingsPtr->getChainBurnIn())
    {
        outStrm << n << '\t' << modelPtr->getLnLike() << '\t' << modelPtr->getGainRate() << '\t' << modelPtr->getLossRate();
        outStrm << '\t' << modelPtr->getGainPathRate() << '\t' << modelPtr->getLossPathRate();
        if (settingsPtr->getModelType() == DISTANCE_NORM ||
            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM)
        {
            outStrm << '\t' << modelPtr->getDistancePower();
        }
        if (settingsPtr->getModelType() == CARRYING_CAPACITY ||
            settingsPtr->getModelType() == CARRYING_CAPACITY_AND_DISTANCE_NORM)
        {
            outStrm << '\t' << modelPtr->getCarryingCapacity();
            outStrm << '\t' << modelPtr->getCarryingPower();
        }
        outStrm << '\t' << modelPtr->getNumGains(0) << '\t' << modelPtr->getNumLosses(0) << std::endl;
    }
}

void Mcmc::samplePaths(int n) {

	if (n == settingsPtr->getHistorySampleFrequency())
	{
        histStrm << "n\tlnL\tnode_index\tstates" << std::endl;
		// tips only recorded once
		for (int i = 0; i < treePtr->getNumNodes(); i++)
		{
			Node* p = treePtr->getNode(i);
			if (p->getNumDescendants() == 0)
			{
				std::vector<bool> s;
				histStrm << n << "\t" << modelPtr->getLnLike() << "\t" << p->getIndex() << "\t";
				for (int j = 0; j < areasPtr->getNumAreas(); j++)
				{
					histStrm << areasPtr->areaState(p->getIndex(), j);
				}
				histStrm << "\n";
			}
		}
	}

	if (n >= settingsPtr->getChainBurnIn())
	{
        // internal nodes recorded every time
        for (int i = 0; i < treePtr->getNumNodes(); i++)
        {
            Node* p = treePtr->getDownPassNode(i);
            if (p->getNumDescendants() != 0)
            {
                std::vector<bool> s;
                s = p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralState();

                histStrm << n << "\t" << modelPtr->getLnLike() << "\t" << p->getIndex() << "\t";
                for (int j = 0; j < s.size(); j++)
                {
                    histStrm << s[j];
                }
                histStrm << "\n";

            }
        }
	}
}


void Mcmc::updateAncestralStateCounts(int n)
{
    if (n >= settingsPtr->getProbBurnIn())
    {
        for (int i = 0; i < treePtr->getNumNodes(); i++)
        {
            Node* p = treePtr->getDownPassNode(i);
            if (p->getNumDescendants() != 0)
            {
                p->getDescendantIndexed(0)->getBranchHistory(0)->updateAncestralStateCounts();
            }
        }
    }
}

void Mcmc::writeAncestralStateFreqs(void)
{

    //int numSamples = (settingsPtr->getChainLength() - settingsPtr->getProbBurnIn()) / settingsPtr->getHistorySampleFrequency();
    for (int i = 0; i < treePtr->getNumNodes(); i++)
    {
        Node* p = treePtr->getDownPassNode(i);
        areaPosteriorStrm << p->getName() << "\t";

        if (p->getNumDescendants() != 0)
        {
            for (int j = 0; j < modelPtr->getNumAreas(); j++)
            {
                // get frequency of occupancy
                std::stringstream ss;
                double x = (double)p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralStateCountForArea(j) / numSamples;
                ss << x;
                areaPosteriorStrm << ss.str() << "\t";
            }
        }
        else
        {
            for (int j = 0; j < modelPtr->getNumAreas(); j++)
            {
                std::string s = ( areasPtr->isInArea(p->getIndex(), j) ? "1" : "0");
                areaPosteriorStrm << s << "\t";
            }
        }
        areaPosteriorStrm << "\n";
    }

}

void Mcmc::writeNhxFile(void)
{
    
    // begin nexus file
    nhxStrm << "#NEXUS" << "\n\n";

    // phylowood settings block
    nhxStrm << "Begin phylowood;\n";
    nhxStrm << "\tdrawtype pie\n";
    nhxStrm << "\tmodeltype biogeography\n";
    nhxStrm << "\tareatype discrete\n";
    nhxStrm << "\tmaptype clean\n";
    nhxStrm << "\tpieslicestyle even\n";
    nhxStrm << "\tpiefillstyle outwards\n";
    nhxStrm << "\ttimestart -" << treePtr->getRoot()->getTime() << "\n";
    nhxStrm << "\tmarkerradius " << 300 << "\n";
    nhxStrm << "\tminareaval " << 0.1 << "\n";
    nhxStrm << "End;\n\n";
    
    // bayarea-fig block
    nhxStrm << "Begin bayarea-fig;\n";
    nhxStrm << "\tmapheight\t100\n";
    nhxStrm << "\tmapwidth\t150\n";
    nhxStrm << "\tcanvasheight\t2000\n";
    nhxStrm << "\tcanvaswidth\t1000\n";
    nhxStrm << "\tminareaval\t0.15\n";
    nhxStrm << "\tareacolors black\n";
    nhxStrm << "\tareatypes";
    for (int i = 0; i < areasPtr->getNumAreas(); i++)
        nhxStrm << " 1";
    nhxStrm << "\n";
    nhxStrm << "\tareanames Default\n";
    nhxStrm << "End;\n\n";
    
    // taxa block
    nhxStrm << "Begin taxa;\n";
    nhxStrm << "\tDimensions ntax=" << areasPtr->getNumTaxa() << ";\n";
    nhxStrm << "\tTaxlabels\n";
    for (int i = 0; i < treePtr->getNumNodes(); i++)
        if (treePtr->getNode(i)->getNumDescendants() == 0)
            nhxStrm << "\t\t" << treePtr->getNode(i)->getName() << "\n";
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";

    // geo block
    nhxStrm << "Begin geo;\n";
    nhxStrm << "\tDimensions ngeo=" << areasPtr->getNumAreas() << ";\n";
    nhxStrm << "\tCoords\n";

    for (int i = 0; i < geoPtr->getNumAreas(); i++)
    {
        nhxStrm << "\t\t" << i << "\t" << geoPtr->getData(0.0, i, 0) << "\t" << geoPtr->getData(0.0, i, 1);
        if (i < (modelPtr->getNumAreas() - 1))
            nhxStrm << ",";
        nhxStrm << "\n";
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";

    // tree block
    nhxStrm << "Begin trees;\n";
    nhxStrm << "\tTranslate\n";
    for (int i = 0; i < treePtr->getNumNodes(); i++)
    {
        if (treePtr->getNode(i)->getNumDescendants() == 0)
        {
            nhxStrm << "\t\t" << treePtr->getNode(i)->getIndex() << "\t" << treePtr->getNode(i)->getName();
            if (i < (treePtr->getNumNodes() - 1))
                nhxStrm << ",";
            nhxStrm << "\n";
        }
    }
    nhxStrm << "\t\t;\n";

    // write tree string
    std::string treeStr = "";
    treeStr = addNodeNhxToString(treePtr->getRoot(), treeStr);
    //std::cout << "nhxStr\n" << treeStr << "\n";
    nhxStrm << "tree TREE1 = " << treeStr << "\n";
    nhxStrm << "End;\n";


}

std::string Mcmc::addNodeNhxToString(Node* p, std::string s)
{

    if (p != NULL)
    {

        // define divergence events
        std::vector<Node*> ndeDescendants = p->getDescendants();
        if (ndeDescendants.size() != 0)
        {
            s += "(";
            //std::cout << s << "\n";
            std::vector<Node*>::iterator it_last = ndeDescendants.end() - 1;
            for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            {
                s = addNodeNhxToString(*it, s);
                //std::cout << s << "\n";
                if (it != it_last)//ndeDescendants.end())
                    s += ",";
                else
                    s+= ")";
                //std::cout << s << "\n";
            }

        }

        // define node & branch values
        std::stringstream ss;
        //int numSamples = (settingsPtr->getChainLength() - settingsPtr->getProbBurnIn()) / settingsPtr->getHistorySampleFrequency();

        // only label tips
        if (p->getNumDescendants() == 0)
            ss << p->getIndex();

        ss << "[&area_pp={";
        for (int i = 0; i < modelPtr->getNumAreas(); i++)
        {
            if (i != 0)
                ss << ",";

            double x = 0.0;
            if (p->getNumDescendants() == 0)
            {
                x = areasPtr->areaState(p->getIndex(), i);
            }
            else
            {
                x = (double)p->getDescendantIndexed(0)->getBranchHistory(0)->getAncestralStateCountForArea(i) / numSamples;
            }
            ss << x;
        }
        ss << "}]:" << p->getLen();
        s += ss.str();
    }

    // string complete
    if (p == treePtr->getRoot())
        s += ";";

    //std::cout << s << "\n";
    return s;
}
