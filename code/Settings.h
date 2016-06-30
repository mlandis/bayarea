#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

	public:
                        Settings(int argc, char* argv[]);
        int             getChainLength(void) { return chainLength; }
        int             getChainBurnIn(void) { return chainBurnIn; }
        int             getProbBurnIn(void) { return probBurnIn; }
        std::string     getAreasFileName(void) { return areaFileName; }
        int             getPrintFrequency(void) { return printFrequency; }
        std::string		getInputFilePath(void) { return inputFilePath; }
        std::string     getOutputFileName(void) { return outputFileName; }
        std::string     getHistoryFileName(void) { return historyFileName; }
        std::string     getAreaPosteriorFileName(void) { return areaPosteriorFileName; }
        std::string     getNhxFileName(void) { return nhxFileName; }
        int             getParameterSampleFrequency(void) { return parameterSampleFrequency; }
        int             getHistorySampleFrequency(void) { return historySampleFrequency; }
        std::string     getTreeFileName(void) { return treeFileName; }
        std::string     getGeoFileName(void) { return geoFileName; }
        int				getSeed(void) { return seed; }
        double			getBetaSteppingStone(void) { return betaSteppingStone; }
        bool			getUseSteppingStone(void) { return useSteppingStone; }

        void            print(void);
        void            setAreaFileName(std::string s) { areaFileName = s; }
        void            setOutputFileName(std::string s) { outputFileName = s; }
        void            setHistoryFileName(std::string s) { historyFileName = s; }
        void            setTreeFileName(std::string s ) { treeFileName = s; }
        void            setGeoFileName(std::string s) { geoFileName = s; }
        void            setAreaPosteriorFileName(std::string s) { areaPosteriorFileName = s; }
        void            setNhxFileName(std::string s) { nhxFileName = s; }
        void			setSeed(int x) { seed = x; }

        double			getCarryingCapacity(void) { return carryingCapacity; }
        double			getCarryingPower(void) { return carryingPower; }
        double			getGeoDistancePower(void) { return geoDistancePower; }
        int             getGeoDistanceType(void) { return geoDistanceType; }
        double			getAreaGainRate(void) { return areaGainRate; }
        double			getAreaLossRate(void) { return areaLossRate; }
        double			getPathGainRate(void) { return pathGainRate; }
        double			getPathLossRate(void) { return pathLossRate; }
        double          getAreaGainPrior(void) { return areaGainPrior; }
        double          getAreaLossPrior(void) { return areaLossPrior; }
        double          getGeoDistancePowerPrior(void) { return geoDistancePowerPrior; }
        bool            isGeoDistancePositive(void) { return geoDistancePositive; }
        int             getModelType(void) { return modelType; }
        bool			getUseSimulate(void) { return useSimulate; }
        bool            getGuessInitialRates(void) { return guessInitialRates; }
        double          getAreaProposalTuner(void) { return areaProposalTuner; }
        double          getRateProposalTuner(void) { return rateProposalTuner; }
        double          getDistanceProposalTuner(void) { return distanceProposalTuner; }
        bool            getUseAuxillarySampling(void) { return useAuxillarySampling; }
        bool            getGeoDistanceTruncate(void) { return geoDistanceTruncate; }
    
    private:
        void			setArguments(int argc, char** argv);
        void            printUsage(void);
    
        int             chainLength;
        int             chainBurnIn;
        int             probBurnIn;

        std::string		outputPrefix;
        std::string     timestampPrefix;
        std::string		bgdFileName;
        std::string     areaFileName;
        std::string     geoFileName;
        std::string     outputFileName;
        std::string     historyFileName;
        std::string     areaPosteriorFileName;
        std::string     nhxFileName;
        std::string     treeFileName;
        std::string		inputFilePath;
        std::string		outputFilePath;

        int             printFrequency;
        int             parameterSampleFrequency;
        int             historySampleFrequency;
        bool            useTimestampPrefix;

        double			carryingPower;
        double			carryingCapacity;
        double          carryingPowerPrior;
        double          carryingCapacityPrior;
    
        double          geoDistancePower;
        int             geoDistanceType;
        double          geoDistancePowerPrior;
        bool            geoDistancePositive;
        bool            geoDistanceTruncate;
    
        double          areaProposalTuner;
        double          rateProposalTuner;
        double          distanceProposalTuner;
        double			areaGainRate;
        double			areaLossRate;
        double			pathGainRate;
        double			pathLossRate;
        double          areaLossPrior;
        double          areaGainPrior;
        bool            useAuxillarySampling;
    
        bool            guessInitialRates;
        int             modelType;
        int				seed;

        bool			useSteppingStone;
        double			betaSteppingStone;

        double          treeScale;

        bool			useSimulate;
};


#endif
