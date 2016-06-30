#ifndef TransitionProbability_H
#define TransitionProbability_H

class Node;



class TransitionProbability {
    
    public:
                                TransitionProbability(Node* n);
                                //TransitionProbability(Node* nf);
                               ~TransitionProbability(void);
        double**                getTiProb(void) { return tiProb; }
        void                    print(void);
        void                    setProbs(double v, double& lambda01, double& lambda10, int numAreas);
        
    private:
        Node*                   nodePtr;
        bool                    scaleFactor;
        int                     numStates;
        double**                tiProb;
};


#endif
