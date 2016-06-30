#include "Areas.h"
#include "GeoCoords.h"
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Settings.h"
#include "Tree.h"
#include <iostream>
#include <string>

class Areas;
class GeoCoords;
class Mcmc;
class Model;
class MbRandom;
class Settings;
class Tree;

int main(int argc, char* argv[]) {

    std::cout << "\nBayArea v1.0.3\n";
    
    // read the user settings
    Settings mySettings(argc, argv);

    mySettings.print();
   
    // instantiate a random number object
    MbRandom myRandom;
    if (mySettings.getSeed() > 0) {
  //  	std::cout << "seed:\t" << mySettings.getSeed() << "\n";
    	myRandom.setSeed(mySettings.getSeed());
    }
    
    // read the file containing the areas
    Areas myAreas( mySettings.getInputFilePath() + mySettings.getAreasFileName() );
    std::cout << "\nReading areas file:\t" << mySettings.getAreasFileName() << "\n";
    std::cout << "\tNumber of taxa  = " << myAreas.getNumTaxa() << "\n";
    std::cout << "\tNumber of areas = " << myAreas.getNumAreas() << "\n";
    myAreas.print();
    
    // read the file containing the fixed tree
    Tree myTree( mySettings.getInputFilePath() + mySettings.getTreeFileName(), myAreas.getTaxonList() );
    std::cout << "\nReading tree file:\t" << mySettings.getTreeFileName() << "\n";
    myTree.print();
    
    // read the file containing the geographical coordinates
    GeoCoords myGeoCoords( &mySettings );
    std::cout << "\nReading geography file:\t" << mySettings.getGeoFileName() << "\n";
    std::cout << myGeoCoords.getDataStr();

    // set up a model
    Model myModel( &myAreas, &myGeoCoords, &myTree, &myRandom, &mySettings );
    
    // run the MCMC analysis
    Mcmc myMarkovChain( &myAreas, &myGeoCoords, &myRandom, &myModel, &myTree, &mySettings );
     
    return 0;
}

