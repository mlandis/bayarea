/*
 * Atlas.cpp
 *
 *  Created on: Jan 28, 2012
 *      Author: mlandis
 */


#include "Atlas.h"
#include "Util.h"


Atlas::Atlas(Settings* sp)
{
	settingsPtr = sp;
	timeStart = 0.0;
	timeEnd = 0.0;

	numLayers = 0;
	numAreas = 0;
	numMetrics = 0;

	//fileNameVector = settingsPtr->getAtlasFileNameVector();
	fileNameVector.push_back(settingsPtr->getInputFilePath() + settingsPtr->getGeoFileName());

	// user specifies maps in fileNameVector
	for (unsigned int i = 0; i < fileNameVector.size(); i++)
	{
		addLayer(fileNameVector[i]);
	}
}

void Atlas::addLayer(std::string fn)
{

	/*
	 e.g. layer input file
	 #	time
	 1.0	2.0	-1.0
	 1.0	1.0	0.0
	 ...

	*/

	// initialize data
	std::cout << "ATLAS: reading \"" << fn << "\"\n";
	std::string layerLine = "";
	std::ifstream layerStream;
	FileMgr* layerFile = new FileMgr(fn);
	if (layerFile->openFile(layerStream) == false)
	{
		std::cerr << "Cannot open file \"" + layerFile->getFileName() + "\"\n";
		exit(1);
	}

	// read layer into data
	bool header = true;
	double t = 1.0;
	//while (layerStream.good())
    while (!layerStream.eof())
	{
		Util::safeGetline(layerStream, layerLine);
		std::stringstream ss(layerLine);
		std::string token;
		std::vector<double> tokens;
		while ( ss >> token )
		{
			if (header == true)
			{
				// time
				ss >> token;
				t = atof(token.c_str());
				if (t < timeStart)
					timeStart = t;
				header = false;
				break;
			}
			else
			{
				double v = atof(token.c_str());
				tokens.push_back(v);
			}
		}
		if (tokens.size() != 0)
		{
			data[t].push_back(tokens);
		}
	}
	numLayers = (int)data.size();
	numAreas = (int)(*data.begin()).second.size();
	numMetrics = (int)(*data.begin()).second.front().size();

	std::cout << "\ttime\t(" << timeStart << "," << timeEnd << ")\n";
	std::cout << "\tnumLayers\t" << numLayers << "\n";
	std::cout << "\tnumAreas\t" << numAreas << "\n";
	std::cout << "\tnumMetrics\t" << numMetrics << "\n";
}

std::string Atlas::getDataStr(void)
{
	std::stringstream ss;
	std::map<double, std::vector<std::vector<double> > >::iterator it_t = data.begin();
	for (; it_t != data.end(); it_t++)
	{
		ss << getDataStr(it_t->first) << "\n";
	}
	return ss.str();
}

std::string Atlas::getDataStr(double t)
{
	std::stringstream ss;
	std::cout << "t:\t" << t << "\tin\t(" << timeStart << ", " << timeEnd << ")\n";
	for (int i = 0; i < numAreas; i++)
	{
		for (int j = 0; j < numMetrics; j++)
		{
			ss << data[t][i][j] << "\t";
		}
		ss << "\n";
	}
	ss << "\n";
	return ss.str();
}

std::string Atlas::getDistanceStr(std::string n, double t)
{
	std::stringstream ss;
	std::cout << "name:\t" << n << "\tt:\t" << t << "\tin\t(" << timeStart << ", " << timeEnd << ")\n";
	for (int i = 0; i < numAreas; i++)
	{
		for (int j = 0; j < numAreas; j++)
		{
			ss << distance[n][t][i][j] << "\t";
		}
		ss << "\n";
	}
	ss << "\n";
	return ss.str();
}



void Atlas::initializeDistance(std::string d)
{
	std::map<double,std::vector<std::vector<double> > >::iterator it_m = data.begin();
	for (; it_m != data.end(); it_m++)
	{
		double t = it_m->first;
		distance[d][t].resize(numAreas);
		for (int i = 0; i < numAreas; i++)
			distance[d][t][i].resize(numAreas);
		calculateDistance(d,t);
	}
}
