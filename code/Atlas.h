/*
 * Atlas.h
 *
 *  Created on: Jan 28, 2012
 *      Author: mlandis
 */

#ifndef ATLAS_H_
#define ATLAS_H_

#include <map>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>

#include "FileMgr.h"
#include "Settings.h"

class FileMgr;
class Settings;

class Atlas {
public:

	Atlas(Settings* sp);

	void addLayer(std::string fp);
	void initializeData(void);

	const std::vector<std::vector<double> >&					getData(double t)				{ return data[t]; }
	const std::vector<double>&									getData(double t, int i)		{ return data[t][i]; }
	double														getData(double t, int i, int j)	{ return data[t][i][j]; }

	const std::map<double, std::vector<std::vector<double> > >& getDistance(std::string d)							{ return distance[d]; }
	const std::vector<std::vector<double> >&					getDistance(std::string d, double t)				{ return distance[d][t]; }
	const std::vector<double>&									getDistance(std::string d, double t, int i)			{ return distance[d][t][i]; }
	double														getDistance(std::string d, double t, int i, int j)	{ return distance[d][t][i][j]; }

	//void setDistance(std::map<double, std::vector<std::vector<double> > >& d) { distance = d; };


	int getNumLayers(void)	{ return numLayers; }
	int getNumAreas(void)	{ return numAreas; }
	int getNumMetrics(void)	{ return numMetrics; }

	std::string getDataStr(void);
	std::string getDataStr(double t);

	std::string getDistanceStr(void);
	std::string getDistanceStr(std::string d);
	std::string getDistanceStr(std::string d, double t);


protected:

	void			initializeDistance(std::string d);
	virtual void	calculateDistance(std::string d, double t) = 0;

	Settings*		settingsPtr;

	// [layer][area][coord]
	std::map<double, std::vector<std::vector<double> > > data;

	// [distance][layer][area][area]
	std::map<std::string, std::map<double, std::vector<std::vector<double> > > > distance;

	std::vector<std::string>				fileNameVector;

	int										numLayers;
	int										numAreas;
	int										numMetrics;

	double									timeStart;
	double									timeEnd;
};

#endif /* ATLAS_H_ */
