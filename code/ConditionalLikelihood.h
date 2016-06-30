#ifndef ConditionalLikelihood_H
#define ConditionalLikelihood_H

class Node;



class ConditionalLikelihood {
    
    public:
                                ConditionalLikelihood(Node* n, int na, int ns);  
                               ~ConditionalLikelihood(void);
        double*                 getCls(void) { return cls; }
        void                    print(void);
        
    private:
        Node*                   nodePtr;
        int                     numStates;
        int                     numAreas;
        double*                 cls;
};


#endif
