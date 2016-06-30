#ifndef Tree_H
#define Tree_H

#include <sstream>
#include <string>
#include <vector>
class MbRandom;
class Node;



class Tree {

    public:
                                Tree(std::string fileName, std::vector<std::string>& taxonNames);  
                                Tree(Tree &t);
                               ~Tree(void);
        int                     getNumNodes(void) { return (int)nodes.size(); }
        int                     getNumTaxa(void) { return numTaxa; }
        double                  getTreeLength(void) { return treeLength; }
        Node*                   getDownPassNode(int i) { return downPassSequence[i]; }
        Node*                   getNode(int i) { return nodes[i]; }
        Node*                   getRoot(void) { return root; }
        std::string             getNewick(void);
        void                    initializeDownPassSequence(void);
        void                    incrementAreaCounts(void);
        Node*                   randomInteriorNode(MbRandom* ranPtr);
        Node*                   randomInteriorNodeByLength(MbRandom* ranPtr);
        void                    showAreaProbs(Node* p, int indent, int numSamples);
        void                    writeTree(Node *p, std::stringstream &ss);
        void                    writeTreeState(std::stringstream &ss);
        void                    print(void);

    private:
        void                    passDown(Node* p);
        void                    readNewickFormattedTree(std::string ts, std::vector<std::string>& taxonNames);
        void                    showNodes(Node* p, int indent);
        int                     numTaxa;
        Node*                   root;
        std::vector<Node*>      nodes;
        std::vector<Node*>      downPassSequence;
        double                  treeLength;
};

#endif
