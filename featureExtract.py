""" Create mz:rt pairs based on overlapping mass and RT windows by looking at
    values across samples and performing a basic binning procedure.
       NOTE: This is a very rough "binning" AUC calculation procedure and more 
            elaborate methods would be better for peak alignment and feature 
            selection

    Assume large files and memory issues w self-joins, perform iteration instead
    
    For each m/z:rt, within sample, sum intensities over that mass tolerance and rt window
    Normalize AUCs by their sample intensity (total sample sum, median-centered)

    Fulgens Consulting
"""

import pandas as pd 
import numpy as np
import sys, os, re
sys.path.append('./src')
import setup_featureEx as util

studyName = 'ST000934'
refDir = 'data/processed/ST000934'
fileKeep = ['S-402', 'S-401', 'S-403', 'U-222', 'U-43', 'U-162', 
            'S-0041', 'U-42', 'U-461', 'U-463', 'U-1122', 'S-391']

def main(refDir, fileKeep, studyName):
    files = os.listdir(refDir)
    if fileKeep is not None:
        fileCheck = [x for x in files if re.sub('\.csv', '', x) in fileKeep]
        if len(fileCheck) == 0:
            print('File filter left no files, returning to fill original set')
        else:
            files = fileCheck
    
    studyCollect = []
    
    for f in files:
        
        df = pd.read_csv('{0}/{1}'.format(refDir, f))
        studyCollect.append(df)
    
    studyCollect = pd.concat(studyCollect)
    studyCollect = util.extendWindowsDuo(df=studyCollect, lowCol='massLower', 
        highCol='massUpper')
    
    collect = [] 
    for _, mini in studyCollect.groupby(['massLower']):
        mini = util.extendWindowsDuo(df=mini, lowCol='rtLower', highCol='rtUpper')
        collect.append(mini)

    df = pd.concat(collect)

    collect = []
    currSamp = None
    for _, mini in df.groupby(['sample', 'massLower', 'rtLower']):
        mz = np.nanmean([mini['massUpper'].mean(), mini['massLower'].mean()])
        rt = np.nanmean([mini['rtLower'].mean(), mini['rtUpper'].mean()])
        auc = np.nansum(mini['intensity'])
        samples = mini['sample'].unique().tolist()
        if len(samples) != 1:
            raise Exception('Only one sample name expected, got {}'.format(samples))
        
        collect.append([mz, rt, auc, samples[0]])
        if samples[0] != currSamp:
            print('{} processed'.format(samples[0]))
            currSamp = samples[0]
        
    collect = pd.DataFrame(collect, columns = ['mz', 'rt', 'auc', 'sample'])
    normCollect = []
    for _, sampDF in collect.groupby('sample'):
        sampDF['auc'] = sampDF['auc'].divide(sampDF['auc'].sum())
        sampDF['auc'] = sampDF['auc'].divide(sampDF['auc'].median())
        normCollect.append(sampDF)
    
    normCollect = pd.concat(normCollect)
    collect.to_csv('data/{}/featureListRaw.csv'.format(studyName), index=False)
    normCollect.to_csv('data/{}/featureListNorm.csv'.format(studyName), index=False)


if __name__ == '__main__':
    main(refDir, fileKeep, studyName)
