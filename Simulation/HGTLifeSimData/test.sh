#!/bin/bash

python simulation.py --outputDirectory "C:\Users\penga\Downloads\Simulation\HGTLifeSimData" --runID testrun --censusTimes 1 500 --mutationRate 0.001 --hgtTranRate 0.1 --hgtConjRate 0.1 --motility 2 --hgtMaxDistance 1 --replDeathRate 0.3 --spaceSize 50 --pGenSize 1000 --foodSources 25,1000 --foodShelfLife 40,80 --fedReplBoost 0.05 --plotYLim 100 --mutationRateChangeMax 0 --replRateChangeMax 0 --deathRateChangeMax 0 --motilityChangeMax 0 --hgtTranRateChangeMax 0 --hgtConjRateChangeMax 0 --geneNumber 20 --genomeType "relatedRandom" --foodGaussMean 0 --foodGaussSigma 1
