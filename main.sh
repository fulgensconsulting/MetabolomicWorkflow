#!bin/bash
python studyQuery.py 
Rscript peakExtract.r
python featureExtract.py
Rscript peakAnalyze.r
