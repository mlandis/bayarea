//#include "Msg.h"
#include "Settings.h"
#include "Util.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>

// modelType aliases
#define INDEPENDENCE 1
#define CARRYING_CAPACITY 2
#define DISTANCE_NORM 3
#define CARRYING_CAPACITY_AND_DISTANCE_NORM 4

Settings::Settings(int argc, char *argv[]) {

    if (argc == 1)
    {
        printUsage();
        return;
    }
  
	/* set default values for parameters */
	seed             = -1;
	areaFileName     = "";
	outputFileName   = "";
	historyFileName  = "";
	areaPosteriorFileName = "";
	nhxFileName = "";
	treeFileName     = "";
	chainLength      = 10000000;
	printFrequency   = 1000;
	parameterSampleFrequency  = 1000;
	historySampleFrequency = 1000;
	chainBurnIn = chainLength * 0.0;
    probBurnIn = chainLength * 0.25;
	geoDistanceType  = 1;
	geoDistancePower = 0.0;
	carryingCapacity = 0.5;
	carryingPower = 1.0;
	areaGainRate = 0.0;
	areaLossRate = 0.0;
    areaGainPrior = 1.0;
    areaLossPrior = 1.0;
    areaProposalTuner = 0.2;
    rateProposalTuner = 0.5;
    distanceProposalTuner = 0.5;
    useAuxillarySampling = false;
    geoDistancePowerPrior = 1.0;
    geoDistancePositive = false;
    geoDistanceTruncate = false;
	pathGainRate = areaGainRate;
	pathLossRate = areaLossRate;

//	modelType = INDEPENDENCE;
//	modelType = CARRYING_CAPACITY;
	modelType = DISTANCE_NORM;
//	modelType = CARRYING_CAPACITY_AND_DISTANCE_NORM;

    guessInitialRates = true;
	useSimulate = false;
	useSteppingStone = false;
	betaSteppingStone = 1.0;

	treeScale = 1.0;

	inputFilePath = "./";
	outputFilePath = "./";
    timestampPrefix = "";
    useTimestampPrefix = true;
	outputPrefix = "";
	bgdFileName = "";
    
    //inputFilePath = "./examples/";
    //areaFileName = "vireya.areas.txt";
    //treeFileName = "vireya.tree.txt";
    //geoFileName = "vireya.geo.txt";

	// overwrite default arguments with command line arguments
	setArguments(argc, argv);
    
    if (seed == -1)
    {
        srand((unsigned int)time(NULL));
        seed = rand() % 100000 + 1;
    }

    // construct output file name
    std::stringstream ss;
    ss << "." << areaFileName;
    if (useSteppingStone == true)
    	ss << ".bss" << betaSteppingStone;
	ss << ".model" << modelType;
	ss << ".seed" << seed;
	if (modelType == DISTANCE_NORM || modelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
		ss << ".dp" << geoDistancePower;
	if (modelType == CARRYING_CAPACITY || modelType == CARRYING_CAPACITY_AND_DISTANCE_NORM)
	{
		ss << ".cc" << carryingCapacity;
		ss << ".cp" << carryingPower;
	}

    if (useTimestampPrefix)
    {
        timestampPrefix = Util::getTime();
        if (outputPrefix != "") outputPrefix += "." + timestampPrefix;
        else outputPrefix += timestampPrefix;
    }
	outputPrefix += ss.str();


    outputFileName = outputFilePath + outputPrefix + ".parameters.txt";
	historyFileName = outputFilePath + outputPrefix + ".area_states.txt";
	areaPosteriorFileName = outputFilePath + outputPrefix + ".area_probs.txt";
	
    nhxFileName = outputFilePath + outputPrefix + ".nhx";


}

void Settings::setArguments(int argc, char** argv)
{
	std::map<std::string, std::string> argMap;
	std::map<std::string, std::string>::iterator argMapIt;
	std::pair<std::map<std::string, std::string>::iterator, bool> ret;

	std::vector<std::string> argVector;

	for (int i = 0; i < argc; i++)
	{
		std::stringstream ss(argv[i]);
		argVector = Util::split(ss.str(), '=');
		if (argVector.size() == 2)
		{
			//std::cout << argVector[0] << "\t" << argVector[1] << "\n";
			ret = argMap.insert(std::pair<std::string, std::string>(argVector[0], argVector[1]));
		}
	}

	std::string argName = "";
	std::string argVal = "";

	for (argMapIt = argMap.begin(); argMapIt != argMap.end(); argMapIt++)
	{

		argName = (*argMapIt).first;
		argVal = (*argMapIt).second;

		//std::cout << "Setting " << argName << " = " << argVal << "\n";

		if (argName == "" || argVal == "")
			; // do nothing

		// INPUT SETTINGS
		else if (argName == "-inputFilePath")
			inputFilePath = argVal;
		else if (argName == "-areaFileName" || argName == "-bgdFileName")
			areaFileName = argVal;
		else if (argName == "-treeFileName")
			treeFileName = argVal;
		else if (argName == "-geoFileName")
		    geoFileName = argVal;


		// OUTPUT SETTINGS
		else if (argName == "-outputFilePath")
			outputFilePath = argVal;
		else if (argName == "-outputFileName")
			outputFileName = argVal;
		else if (argName == "-outputPrefix" || argName == "-simName")
			outputPrefix = argVal;
        else if (argName == "-outputTimestamp")
            useTimestampPrefix = Util::stringToBool(argVal);
		else if (argName == "-printFrequency")
			printFrequency = Util::stringToInt(argVal);

		// MODEL SETTINGS
		else if (argName == "-seed")
			seed = Util::stringToInt(argVal);
		else if (argName == "-modelType")
			modelType = Util::stringToInt(argVal);
		else if (argName == "-geoDistanceType")
			geoDistanceType = Util::stringToInt(argVal);
		else if (argName == "-geoDistancePower")
			geoDistancePower = Util::stringToDouble(argVal);
        else if (argName == "-geoDistancePowerPositive")
            geoDistancePositive = Util::stringToBool(argVal);
        else if (argName == "-geoDistanceTruncate")
            geoDistanceTruncate = Util::stringToBool(argVal);
		else if (argName == "-carryingCapacity")
			carryingCapacity = Util::stringToDouble(argVal);
		else if (argName == "-carryingPower")
			carryingPower = Util::stringToDouble(argVal);
		else if (argName == "-useSteppingStone")
			useSteppingStone = Util::stringToBool(argVal);
		else if (argName == "-betaSteppingStone")
			betaSteppingStone = Util::stringToDouble(argVal);
        else if (argName == "-distancePowerPrior")
            geoDistancePowerPrior = Util::stringToDouble(argVal);
        else if (argName == "-gainPrior")
            areaGainPrior = Util::stringToDouble(argVal);
        else if (argName == "-lossPrior")
            areaLossPrior = Util::stringToDouble(argVal);
        else if (argName == "-guessInitialRates")
            guessInitialRates = Util::stringToBool(argVal);
        else if (argName == "-areaProposalTuner")
        {
            areaProposalTuner = Util::stringToDouble(argVal);
            if (areaProposalTuner > 1.0 || areaProposalTuner < 0.0)
                std::cout << "WARNING: areaProposalTuner must be between 0 and 1. Setting to 0.2.\n";
        }
        else if (argName == "-rateProposalTuner")
            rateProposalTuner = Util::stringToDouble(argVal);
        else if (argName == "-distanceProposalTuner")
            distanceProposalTuner = Util::stringToDouble(argVal);
        else if (argName == "-useAuxillarySampling")
            useAuxillarySampling = Util::stringToBool(argVal);

		// MCMC SETTINGS
		else if (argName == "-chainLength")
			chainLength = Util::stringToInt(argVal);
		else if (argName == "-chainBurnIn")
			chainBurnIn = Util::stringToDouble(argVal);
        else if (argName == "-probBurnIn")
			probBurnIn = Util::stringToDouble(argVal);
		else if (argName == "-parameterSampleFrequency" || argName == "-parmSampleFrequency")
			parameterSampleFrequency = Util::stringToInt(argVal);
		else if (argName == "-historySampleFrequency" || argName == "-pathSampleFrequency")
			historySampleFrequency = Util::stringToInt(argVal);

	}
}


void Settings::print(void) {

    std::cout << "Input/Output settings:" << std::endl;
	std::cout << "   * Area file name                 = " << areaFileName    << std::endl;
    std::cout << "   * Geography file name            = " << geoFileName    << std::endl;
	std::cout << "   * Tree file name                 = " << treeFileName    << std::endl;
    std::cout << "   * Input file path                = " << inputFilePath << std::endl;
    std::cout << "   * Output prefix                  = " << outputPrefix << std::endl;
    std::cout << "   * Output timestamp               = " << timestampPrefix << std::endl;
    std::cout << "   * Output file path               = " << outputFilePath    << std::endl;
	std::cout << "   * Parameter sample frequency     = " << parameterSampleFrequency << std::endl;
    std::cout << "   * History sample frequency       = " << historySampleFrequency << std::endl;
    std::cout << "   * Print frequency                = " << printFrequency  << std::endl;
	std::cout << std::endl;
    
	std::cout << "Analysis settings:" << std::endl;
    std::cout << "   * Random number seed             = " << seed << std::endl;
    std::cout << "   * Chain length                   = " << chainLength << std::endl;
	std::cout << "   * Chain burn-in                  = " << chainBurnIn << std::endl;
    std::cout << "   * Probability file burn-in       = " << probBurnIn << std::endl;
    std::cout << "   * Model type                     = " << modelType << std::endl;
    std::cout << "   * Area gain prior                = " << areaGainPrior << std::endl;
    std::cout << "   * Area loss prior                = " << areaLossPrior << std::endl;
    std::cout << "   * Distance power prior           = " << geoDistancePowerPrior << std::endl;
    std::cout << "   * Distance power positive        = " << Util::boolToString(geoDistancePositive) << std::endl;
    std::cout << "   * Geo. distance truncate         = " << Util::boolToString(geoDistanceTruncate) << std::endl;
    std::cout << "   * Guess initial rates            = " << Util::boolToString(guessInitialRates) << std::endl;
    std::cout << "   * Area proposal tuning value     = " << areaProposalTuner << std::endl;
    std::cout << "   * Rate proposal tuning value     = " << rateProposalTuner << std::endl;
    std::cout << "   * Distance proposal tuning value = " << distanceProposalTuner << std::endl;
    std::cout << "   * Use auxillary sampling         = " << useAuxillarySampling << std::endl;
    std::cout << std::endl;

}

void Settings::printUsage(void) {

	std::cout << std::endl;
    std::cout << "Usage:" << std::endl;
	std::cout << "   -areaFileName <FILE NAME>          : Input area file name" << std::endl;
    std::cout << "   -geoFileName <FILE NAME>           : Input geography file name" << std::endl;
    std::cout << "   -treeFileName <FILE NAME>          : Input tree file name" << std::endl;
	std::cout << "   -inputFilePath <FILE NAME>         : Input file directory" << std::endl;
    std::cout << "   -outputPrefix <STRING>             : Identifying output prefix" << std::endl;
    std::cout << "   -outputTimestamp <BOOL>            : Append timestamp string to outputPrefx (T,F)" << std::endl;
    std::cout << "   -outputFilePath <FILE NAME>        : Output file directory" << std::endl;
    std::cout << "   -parameterSampleFrequency <NUMBER> : MCMC sample frequency of model parameters" << std::endl;
    std::cout << "   -historySampleFrequency <NUMBER>   : MCMC sample frequency of biogeographic histories" << std::endl;
    std::cout << "   -printFrequency <NUMBER>           : Stdout frequency of MCMC state" << std::endl;
    std::cout << "   -seed <NUMBER>                     : Random number seed" << std::endl;
	std::cout << "   -chainLength <NUMBER>              : Number of MCMC cycles" << std::endl;
    std::cout << "   -chainBurnIn <NUMBER>              : First MCMC cycle sampling point for .mcmc and .asr files" << std::endl;
    std::cout << "   -probBurnIn <NUMBER>               : First MCMC cycle sampling point for .asp and .nhx files" << std::endl;
    std::cout << "   -modelType <NUMBER>                : Model type (1=INDEPENDENCE, 3=DISTANCE)" << std::endl;
    std::cout << "   -gainPrior <NUMBER>                : Scale parameter for half-Cauchy prior on area gain rate" << std::endl;
    std::cout << "   -lossPrior <NUMBER>                : Scale parameter for half-Cauchy prior on area loss rate" << std::endl;
    std::cout << "   -distancePowerPrior <NUMBER>       : Scale parameter for Cauchy prior on distance power" << std::endl;
    std::cout << "   -geoDistancePowerPositive <BOOL>   : Forbid negative distance power parameters (T=enable, F=disable)" << std::endl;
    std::cout << "   -geoDistanceTruncate <BOOL>        : Faster distance computations by approximation (T=enable, F=disable)" << std::endl;
    std::cout << "   -guessInitialRates <BOOL>          : Initialize gain/loss rate with heuristic (T=enable, F=draw from prior)" << std::endl;
    std::cout << "   -areaProposalTuner <NUMBER>        : Tunes avg prob of resampling history for an area (must be from 0 to 1)" << std::endl;
    std::cout << "   -rateProposalTuner <NUMBER>        : Tunes scaler proposal for gain/loss rates" << std::endl;
    std::cout << "   -distanceProposalTuner <NUMBER>    : Tunes scaler proposal for distance power" << std::endl;
    std::cout << "   -useAuxillarySampling <BOOL>       : Use auxillary \"*\" path sampling, disabled by default (T=enable, F=disable)" << std::endl;
    
	std::cout << std::endl;
	std::cout << "Example:" << std::endl;
	std::cout << "   ./bayarea -areaFileName <area file> -treeFileName <tree file> -geoFileName <geo file>" << std::endl << std::endl;
	exit(1);
}
