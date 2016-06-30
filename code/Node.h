#ifndef Node_H
#define Node_H

#include <string>
#include <vector>
class BranchHistory;
class ConditionalLikelihood;
class TransitionProbability;



class Node {
    
    public:
                                    Node(void);
        void                        addDescendant(Node* p) { descendants.push_back(p); }
        void                        attachBranchHistory(BranchHistory* h, int i) { historyPtr[i] = h; }
        void                        attachConditionalLikelihoods(ConditionalLikelihood* c) { condLikePtr = c; }
        void                        attachTransitionProbability(TransitionProbability* p) { tiProbPtr = p; }
        void                        clearDescendants(void);
        BranchHistory*              getBranchHistory(int i) { return historyPtr[i]; }
        ConditionalLikelihood*      getConditionalLikelihood(void) { return condLikePtr; }
        TransitionProbability*      getTransitionProbability(void) { return tiProbPtr; }
        std::string                 getName(void) { return name; }
        int                         getNumAreas(void) { return (int)areaCounts.size(); }
        double                      getTime(void) { return time; }
        int                         getIndex(void) { return index; }
        Node*                       getAncestor(void) { return ancestor; }
        int                         getCountForArea(int i) { return areaCounts[i]; }
        Node*                       getDescendantIndexed(int idx) { return descendants[idx]; }
        std::vector<Node*>&         getDescendants(void) { return descendants; }
        double                      getLen(void) { return len; }
        int                         getNumDescendants(void) { return (int)descendants.size(); }
        void                        incrementArea(int i) { areaCounts[i]++; }
        bool                        isTip(void);
        void                        resizeAreaCountsTo(int n);
        void                        setName(std::string s) { name = s; }
        void                        setTime(double x) { time = x; }
        void                        setIndex(int x) { index = x; }
        void                        setAncestor(Node* p) { ancestor = p; }
        void                        setLen(double x) { len = x; }
        
    private:
        std::string                 name;
        double                      time;
        double                      len;
        int                         index;
        int                         curHistory;
        BranchHistory*              historyPtr[2];
        ConditionalLikelihood*      condLikePtr;
        TransitionProbability*      tiProbPtr;
        Node*                       ancestor;
        std::vector<Node*>          descendants;
        std::vector<int>            areaCounts;
};

#endif
