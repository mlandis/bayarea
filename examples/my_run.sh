./bayarea -areaFileName=vireya.areas.txt -geoFileName=vireya.geo.txt -treeFileName=vireya.tree.txt -inputFilePath=./examples/ -outputFilePath=./output/ -outputTimestamp=T -outputPrefix=my_run -parameterSampleFrequency=1000 -historySampleFrequency=1000 -printFrequency=1000 -chainLength=10000000 -chainBurnIn=0 -probBurnIn=1000000 -modelType=3 -gainPrior=0.2 -lossPrior=0.2 -useAuxiliarySampling=F -geoDistancePowerPositive=F -distancePowerPrior=0.2 -geoDistanceTruncate=F -guessInitialRates=T -areaProposalTuner=0.5 -rateProposalTuner=0.5 -distanceProposalTuner=0.5