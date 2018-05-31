/*
 * GeoCoords.cpp
 *
 *  Created on: Jan 29, 2012
 *      Author: mlandis
 */

#include "GeoCoords.h"

// constants
#define EARTHRADIUSKM 6371
#define PI 3.14159265


GeoCoords::GeoCoords(Settings* sp) : Atlas(sp) {
	
	distancePower = settingsPtr->getGeoDistancePower();
	initializeDistance("haversineRaw");
    initializeDistance("haversinePower");
    initializeDistance("haversinePowerInverse");
    initializeDistance("haversinePowerInverseTruncated");
    
    truncateThreshhold = 1e-40;
    //

}

void GeoCoords::calculateDistance(std::string d, double t)
{
	//std::cout << "GEOCOORDS: calculating distances\n";
	for (int i = 0; i < numAreas; i++)
	{
        bool truncate = false;
		for (int j = 0; j < numAreas; j++)
		{
            
			if (d == "haversineRaw")
			{
                double lat0 = data[t][i][0] * PI/180.0;
                double lon0 = data[t][i][1] * PI/180.0;
                double lat1 = data[t][j][0] * PI/180.0;
                double lon1 = data[t][j][1] * PI/180.0;
                
				double r = 0.0;
                double dlat = lat0 - lat1;
				double dlon = lon0 - lon1;
				double sindlat = sin(dlat/2);
				double sindlon = sin(dlon/2);
				double a = sqrt(sindlat * sindlat + cos(lat0) * cos(lat1) * sindlon * sindlon);
				if (a < 0.9995)
					r = 2 * asin(a) * EARTHRADIUSKM;
				else
					r = 2 * asin(1.0) * EARTHRADIUSKM;
				//std::cout << i << " " << j << " " << r << "\n";
				distance[d][t][i][j] = r;
			}
			else if (d == "haversinePower")
			{
				distance[d][t][i][j] = pow(distance["haversineRaw"][t][i][j], distancePower);
				//std::cout << distance[d][t][i][j] << "\n";
			}
			else if (d == "haversinePowerInverse")
			{
				distance[d][t][i][j] = pow(distance["haversineRaw"][t][i][j], -distancePower);
				//std::cout << distance[d][t][i][j] << "\n";
			}
            else if (d == "haversinePowerInverseTruncated")
            {
                int k = orderedDistanceIdx[t][i][j];
                
                if (i == k)
                    distance[d][t][i][k] = 0.0;
                
                else if (truncate)
                    distance[d][t][i][k] = truncateThreshhold;
                
                else
                {
                    double v = pow(distance["haversineRaw"][t][i][k], -distancePower);
                    if (v < truncateThreshhold)
                    {
                        distance[d][t][i][k] = truncateThreshhold;
                        truncate = true;
                    }
                    else
                    {
                        distance[d][t][i][k] = v;
                    }
                }
            }
			else
			{
				std::cout << "GEOCOORDS: Error: distance measure \"" << d << "\" undefined.\n";
			}
		}
	}
    if (d == "haversineRaw")
        findOrderedDistanceIdx(t);
}

void GeoCoords::setTruncateThreshhold(double th)
{
    truncateThreshhold = th;
}

void GeoCoords::findOrderedDistanceIdx(double t)
{
    
//    double minDist = 0;
    orderedDistanceIdx[t].resize(numAreas);
    for (size_t i = 0; i < numAreas; i++)
    {
        //orderedDistanceIdx[i].resize(numAreas,0);
        for (size_t j = 0; j < numAreas; j++)
        {
            double d = distance["haversineRaw"][t][i][j];
            
            size_t k = 0;
            for (k = 0; k < orderedDistanceIdx[t][i].size(); k++)
            {
                if (d < distance["haversineRaw"][t][i][orderedDistanceIdx[t][i][k]])
                    break;
            }
            std::vector<int>::iterator it = orderedDistanceIdx[t][i].begin();
            orderedDistanceIdx[t][i].insert(it + k, (int)j);
        }
    }
}

void GeoCoords::changeDistancePower(std::string d, double t, double x)
{
    distancePower = x;
    
    if (d == "haversinePower")
    {
        for (int i = 0; i < numAreas; i++)
            for (int j = 0; j < numAreas; j++)
                distance[d][t][i][j] = pow(distance["haversineRaw"][t][i][j], distancePower);
    }

    else if (d == "haversinePowerInverse")
    {
        for (int i = 0; i < numAreas; i++)
            for (int j = 0; j < numAreas; j++)
                distance[d][t][i][j] = pow(distance["haversineRaw"][t][i][j], -distancePower);
    }
    
    else if (d == "haversinePowerInverseTruncated")
    {
        double factor = (distancePower > 0.0 ? 1e100 : 1e-100);
        for (int i = 0; i < numAreas; i++)
        {
            for (int j = 0; j < numAreas; j++)
            {
                // which dispersal event has the highest rate?
                //std::cout << distance["haversineRaw"][t][i][j] << "\n";
                if (i == j)
                    continue;
                else if (distancePower > 0.0 && distance["haversineRaw"][t][i][j] < factor)
                    factor = distance["haversineRaw"][t][i][j];
                else if (distancePower < 0.0 && distance["haversineRaw"][t][i][j] > factor)
                    factor = distance["haversineRaw"][t][i][j];
            }
        }
        double factorPower = pow(factor,-distancePower);
        //std::cout << factor << " " << factorPower << "\n";
        
        //factorPower = 1.0;
        
        for (int i = 0; i < numAreas; i++)
        {
            bool truncate = false;
     
            
            for (int k = 0; k < numAreas; k++)
            {
                int j = orderedDistanceIdx[t][i][k];
                
                if (i == j)
                {
                    distance[d][t][i][j] = 0.0;
                }
                else if (truncate)
                {
                    //std::cout << "truncated!\n";
                    distance[d][t][i][j] = truncateThreshhold;
                }
                else
                {
                    double v = pow(distance["haversineRaw"][t][i][j], -distancePower);

                    // truncate values if the ratio of event i,j is very small with respect to the largest rate event
                    if (v/factorPower < truncateThreshhold)
                    {
                        distance[d][t][i][j] = truncateThreshhold;
                        truncate = true;
                    }
                    else
                    {
                        distance[d][t][i][j] = v;
                    }
                }
            }
        }
    }
}
