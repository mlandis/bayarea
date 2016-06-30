#include <iostream>
#include "BranchHistory.h"
#include "ConditionalLikelihood.h"
#include "Node.h"
#include "TransitionProbability.h"



Node::Node(void) {

    name = "";
    time = 0.0;
    index = 0;
    condLikePtr = NULL;
    tiProbPtr = NULL;
    historyPtr[0] = NULL;
    historyPtr[1] = NULL;
    ancestor = NULL;
    descendants.empty();
}

void Node::clearDescendants(void) {

    for (std::vector<Node*>::iterator p = descendants.begin(); p != descendants.end(); p++)
        delete (*p);
    descendants.clear();
    if (condLikePtr != NULL)
        delete condLikePtr;
    if (tiProbPtr != NULL)
        delete tiProbPtr;
    if (historyPtr[0] != NULL)
        delete historyPtr[0];
    if (historyPtr[1] != NULL)
        delete historyPtr[1];
}

bool Node::isTip(void) {

    if ( descendants.size() == 0 )
        return true;
    return false;
}

void Node::resizeAreaCountsTo(int n) {

    areaCounts.resize(n);
    for (int i=0; i<areaCounts.size(); i++)
        areaCounts[i] = 0;
}




