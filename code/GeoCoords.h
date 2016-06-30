/*
 * GeoCoords.h
 *
 *  Created on: Jan 29, 2012
 *      Author: mlandis
 */

#ifndef GEOCOORDS_H_
#define GEOCOORDS_H_

#include "Atlas.h"
#include "Settings.h"

#include <math.h>

class Settings;

class GeoCoords: public Atlas {
public:
	GeoCoords(Settings*);
	double getDistancePower(void) { return distancePower; }
	void changeDistancePower(std::string d, double t, double x);
    void setTruncateThreshhold(double th);

private:
    void findOrderedDistanceIdx(double t);
	void calculateDistance(std::string d, double t);
	double distancePower;
    double truncateThreshhold;
    std::map<double,std::vector<std::vector<int> > > orderedDistanceIdx;
};

#endif /* GEOCOORDS_H_ */
