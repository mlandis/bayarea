#include "TransitionProbability.h"
#include <cmath>
#include <iomanip>
#include <iostream>





TransitionProbability::TransitionProbability(Node* n) {
    

    nodePtr = n;
    numStates = 2;
    scaleFactor = 1.0;
    
    tiProb = new double*[numStates];
    tiProb[0] = new double[numStates*numStates];
    for (int i=1; i<numStates; i++)
        tiProb[i] = tiProb[i-1] + numStates;
    for (int i=0; i<numStates; i++)
        for (int j=0; j<numStates; j++)
            tiProb[i][j] = 0.0;
}

TransitionProbability::~TransitionProbability(void) {
    
    delete [] tiProb[0];
    delete [] tiProb;
}

void TransitionProbability::print(void) {
    
    for (int i=0; i<numStates; i++)
        {
        for (int j=0; j<numStates; j++)
            std::cout << std::fixed << std::setprecision(6) << tiProb[i][j] << " ";
        std::cout << std::endl;
        }
}

void TransitionProbability::setProbs(double v, double& lambda01, double& lambda10, int numAreas) {


	double scaleFactor = 1.0; // (lambda01 + lambda10) / (2.0 * lambda01 * lambda10);
	//scaleFactor /= (double)numAreas;       // 09/30/12 MJL: don't think this is needed under current model parameterization...
	double expPart = exp( -(lambda01 + lambda10) * scaleFactor * v );
	double pi0 = lambda10 / (lambda01 + lambda10);
	double pi1 = 1.0 - pi0;
	tiProb[0][0] = pi0 + pi1 * expPart;
	tiProb[0][1] = pi1 - pi1 * expPart;
	tiProb[1][0] = pi0 - pi0 * expPart;
	tiProb[1][1] = pi1 + pi0 * expPart;

#if 0
	double scaleFactor = 1.0;
	if (scaleFactor == true)
	{
		double scaleFactor = (lambda01 + lambda10) / (2.0 * lambda01 * lambda10);
		scaleFactor /= (double)numAreas;
	}
    double expPart = exp( -(lambda01 + lambda10) * scaleFactor * v );
    double pi0 = lambda10 / (lambda01 + lambda10);
    double pi1 = 1.0 - pi0;
    tiProb[0][0] = pi0 + pi1 * expPart;
    tiProb[0][1] = pi1 - pi1 * expPart;
    tiProb[1][0] = pi0 - pi0 * expPart;
    tiProb[1][1] = pi1 + pi0 * expPart;
#endif
    
}
