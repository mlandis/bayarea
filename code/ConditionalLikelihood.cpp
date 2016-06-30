#include "ConditionalLikelihood.h"
#include <iomanip>
#include <iostream>



ConditionalLikelihood::ConditionalLikelihood(Node* n, int na, int ns) {
    
    nodePtr   = n;
    numAreas  = na;
    numStates = ns;
    
    cls = new double[numAreas * numStates];
    for (int i=0; i<numAreas*numStates; i++)
        cls[i] = 0.0;
}

ConditionalLikelihood::~ConditionalLikelihood(void) {
    
    delete [] cls;
}

void ConditionalLikelihood::print(void) {

    int k = 0;
    for (int aIdx=0; aIdx<numAreas; aIdx++)
        {
        std::cout << "(";
        for (int s=0; s<numStates; s++)
            {
            std::cout << std::fixed << std::setprecision(4) << cls[k++];
            if (s < numStates-1)
                std::cout << " ";
            }
        std::cout << ") ";
        }
}