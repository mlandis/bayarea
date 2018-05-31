#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "BranchHistory.h"
#include "MbRandom.h"
//#include "Msg.h"
#include "Node.h"
#include "Tree.h"
#include "Util.h"



Tree::Tree(std::string fileName, std::vector<std::string>& taxonNames) {
    
	treeLength = 0.0;
    std::string newickStr = Util::getLineFromFile( fileName, 1 );
    readNewickFormattedTree(newickStr, taxonNames);
}

Tree::Tree(Tree& t) {
    
}

Tree::~Tree(void) {
    
    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        delete (*it);
}

std::string Tree::getNewick(void) {
    
    return "";
}

void Tree::incrementAreaCounts(void) {
    
    for (int n=0; n<getNumNodes(); n++)
        {
        Node* p = getDownPassNode(n);
        if (p->getNumDescendants() != 0)
            {
            BranchHistory* h = p->getDescendantIndexed(0)->getBranchHistory(0);
            std::vector<bool> a = h->getAncestralState();
            if (p->getNumAreas() == 0)
                p->resizeAreaCountsTo( (int)a.size() );
            for (int i=0; i<a.size(); i++)
                {
                if (a[i] == true)
                    p->incrementArea(i);
                }
            }
        }
}

void Tree::initializeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root);
}

void Tree::passDown(Node* p) {
    
    if (p != NULL)
        {
        std::vector<Node*> ndeDescendants = p->getDescendants();
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            passDown(*it);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {
    
    showNodes(root, 0);
}

Node* Tree::randomInteriorNode(MbRandom* ranPtr) {
    
    Node* p = NULL;
    while (p == NULL)
        {
        p = nodes[(int)(ranPtr->uniformRv()*nodes.size())];
        if (p->getNumDescendants() == 0)
            p = NULL;
        }
    return p;
}

Node* Tree::randomInteriorNodeByLength(MbRandom* ranPtr) {

	Node* p = NULL;
	double u = ranPtr->uniformRv() * treeLength;
	for (int i = 0; i < nodes.size(); i++)
	{
		p = downPassSequence[i];
		u -= p->getLen();
		if (u < 0.0)
			return p;
	}
	return NULL;
}

void Tree::readNewickFormattedTree(std::string ts, std::vector<std::string>& taxonNames) {
    
	/* parse the tree string, and put each token into a vector of strings */
	std::vector<std::string> parsedNewick;
	std::string temp = "";
	bool readingBrlen = false;
	int nt = 0;
	for (int i=0; i<ts.size(); i++)
		{
		char c = ts[i];
		if ( c == ' ' )
			continue;
		if ( c == '(' || c == ')' || c == ',' || c == ':' || c == ';' )
			{
			temp = c;
			parsedNewick.push_back( temp );
			if ( c == ':' )
				readingBrlen = true;
			else
				readingBrlen = false;
			}
		else
			{
			/* the character is part of a taxon name */
			int j = i;
			std::string taxonName = "";
			while ( ts[j] != '(' && ts[j] != ')' && ts[j] != ',' && ts[j] != ':' && ts[j] != ';' )
				{
				taxonName += ts[j];
				j++;
				}
			parsedNewick.push_back( taxonName );
			i = j - 1;
			if ( readingBrlen == false )
				nt++;
			readingBrlen = false;
			}
		if ( c == ';' )
			break;
		}
#   if 0
    for (std::vector<std::string>::iterator p = parsedNewick.begin(); p != parsedNewick.end(); p++)
        std::cout << (*p) << std::endl;
#   endif
    
	/* check that the number of taxa in the tree description is the same as the
     number of taxa in the alignment */
	if ( nt != taxonNames.size() )
	{
		std::cout << nt << "\t" << taxonNames.size() << "\n";
        std::cout << "The tree file is not of the right size\n";
	}

    /* build-up the tree */
    Node* p = NULL;
    readingBrlen = false;
    int intNodeIdx = (int)taxonNames.size();
    for (std::vector<std::string>::iterator t=parsedNewick.begin(); t != parsedNewick.end(); t++)
		{
		//std::cout << (*t) << std::endl;
		if ( (*t) == "(" )
			{
			/* add a new interior node */
			if (p == NULL)
				{
                p = new Node;
                p->setIndex(intNodeIdx++);
                std::stringstream ss;
                ss << intNodeIdx;
                p->setName("Taxon_" + ss.str());
                nodes.push_back(p);
                root = p;
				}
			else
				{
				Node* q = new Node;
                q->setIndex(intNodeIdx++);
                q->setAncestor(p);
                std::stringstream ss;
                ss << intNodeIdx;
                q->setName("Taxon_" + ss.str());
                nodes.push_back(q);
                p->addDescendant(q);
                p = q;
				}
			readingBrlen = false;
			}
		else if ( (*t) == ")" )
			{
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            else
                std::cout << "Problem going down tree\n";
			readingBrlen = false;
			}
		else if ( (*t) == "," )
			{
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            else
                std::cout << "Problem going down tree\n";
			readingBrlen = false;
			}
		else if ( (*t) == ":" )
			{
			readingBrlen = true;
			}
		else if ( (*t) == ";" )
			{
			/* We are at the end of the tree description. I guess we don't have
             to do anything. */
			}
		else
			{
			if (readingBrlen == false)
				{
                std::string tipName = (*t);
                int tipIdx = -1;
                for (int i=0; i<taxonNames.size(); i++)
                    {
                    if ( tipName == taxonNames[i] )
                        {
                        tipIdx = i;
                        break;
                        }
                    }
                if (tipIdx == -1)
                    std::cout << "Could not find taxon " + tipName + " in the list of taxon names\n";
                
				Node* q = new Node;
                nodes.push_back(q);
                q->setIndex(tipIdx);
                q->setAncestor(p);
                q->setName(taxonNames[tipIdx]);
                p->addDescendant(q);
                p = q;
				}
			else
				{
				/* reading a branch length */
				double x = 0.0;
				std::istringstream buf(*t);
				buf >> x;
				if (x < 0.00001)
					x = 0.0001;
                //p->setTime(x);
                p->setLen(x);
				readingBrlen = false;
				treeLength += x;
				}
			}
		}
    
    // initialize the down pass sequence
    initializeDownPassSequence();

    // set the node times
    for (int n=0; n<getNumNodes(); n++)
        {
        Node* p = getDownPassNode(n);
        if (p != NULL)
            {
            if (p->getNumDescendants() == 0)
                {
                p->setTime(0.0);
                }
            else
                {
                double x = 0.0;
                for (int i=0; i<p->getNumDescendants(); i++)
                    {
                    if (p->getDescendantIndexed(i)->getTime() + p->getDescendantIndexed(i)->getLen() > x)
                        x = p->getDescendantIndexed(i)->getTime() + p->getDescendantIndexed(i)->getLen();
                    }
                p->setTime(x);
                }
            }
        }
    
    // initialize root length to be twice the tree height when left unassigned
    if (root->getLen() == 0.0)
        root->setLen( 2 * root->getTime());
    
    numTaxa = (int)taxonNames.size();
}

void Tree::showAreaProbs(Node* p, int indent, int numSamples) {
    
	if (p != NULL)
		{
        std::vector<Node*> ndeDescendants = p->getDescendants();
		for (int i=0; i<indent; i++)
			std::cout << " ";
        std::cout << p->getIndex();
        std::cout << " (";
        std::cout << p->getTime() << ") ";
        if (p->getNumDescendants() == 0)
            std::cout << " (" << p->getName() << ") ";

        for (int i=0; i<p->getNumAreas(); i++)
            std::cout << std::fixed << std::setprecision(3) << (double)(p->getCountForArea(i))/numSamples << " ";

        std::cout << std::endl;
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            showAreaProbs(*it, indent+2, numSamples);
		}
    
}

void Tree::showNodes(Node* p, int indent) {
    
	if (p != NULL)
		{
        std::vector<Node*> ndeDescendants = p->getDescendants();
		for (int i=0; i<indent; i++)
			std::cout << " ";
        std::cout << p->getIndex();
        std::cout << " ( ";
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            std::cout << (*it)->getIndex() << " ";
        std::cout << ") " << p->getLen() << "/" << p->getTime();
        if (p->getNumDescendants() == 0)
            std::cout << " (" << p->getName() << ") ";
		if (p == root)
			std::cout << " <- Root";
        std::cout << std::endl;
        for (std::vector<Node*>::iterator it = ndeDescendants.begin(); it != ndeDescendants.end(); it++)
            showNodes(*it, indent+2);
		}
    
}

void Tree::writeTreeState(std::stringstream& ss)
{
	for (int n = 0; n < getNumNodes(); n++)
	{
		Node* p = getDownPassNode(n);
		if (p != NULL)
		{
			ss << p->getName() << "\t";
			std::vector<bool> s = p->getBranchHistory(0)->getAncestralState();
			for (int i = 0; i < s.size(); i++)
			{
				ss << s[i];
			}
			ss << "\n";
		}
	}
}

void Tree::writeTree(Node* p, std::stringstream& ss) {
}


