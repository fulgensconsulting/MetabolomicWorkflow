# Quick analysis of group variance from mz_rt Features through PCA and outlier detection
# ... Subsequently, trim features down to most "important" on variable subset selection
# ... and plot these loadings to the dimensionality reduction approach 
#
# Fulgens Consulting 

library(readr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(reshape2,quietly=TRUE)
library(cowplot, quietly=TRUE)

studyName = c('ST000934')
filePath = file.path('src')
source(file.path(filePath, 'setup_msprocess.R'))
df <- read_csv(file.path('data/processed', studyName, 'featureListNorm.csv'))
colnames(df) <- c('mz', 'rt', 'AUC', 'Sample')

df$mz <- round(df$mz, 5)
df$rt <- round(df$rt / 60, 2)
df$mz_rt <- paste0(as.character(df$mz), '_', as.character(df$rt))
mzTable <- as.data.frame(acast(df, Sample ~ mz_rt, value.var = 'AUC'))
mzTable$Sample <- rownames(mzTable)
mzPCA <- plotPCA(data=mzTable, confidence=95, nonDatCols=c('Sample'), minVals=0.25)

mzTable <- mzTable[mzTable$Sample != c('U-1122'), ]
samples <- read_csv(file.path('ref', studyName, 'Samples.csv')) %>% select(Sample, Status, Stage)
mzTable <- merge(mzTable, samples, by = c('Sample'))
mzTable$Sample = paste0(mzTable$Status, '_', mzTable$Stage)
mzSub <- mzTable %>% select(-Status, -Stage)

mzTrimPCA <- plotPCA(data=mzSub, confidence=95, nonDatCols=c('Sample'), 
    minVals=0.25, ggTitle=c('Group Sample PCA after Outlier Removal'))
mzTrimPCA

trimFeatures <- ElasticNetVariableSelection(data=mzSub, alpha=0.5, minVals=0.25)
mzTrim <- mzSub %>% select(trimFeatures)
mzTrim$Sample <- mzSub$Sample
pcaFit <- fitPCA(data=mzTrim, nonDatCols=c('Sample'), returnFit=TRUE)
pcaLoadings <- PCALoadingsPlot(origDF=df, pcaFitObject=pcaFit, varNames=TRUE, 
    adjustScale=TRUE, ggTitle=c('Discriminatory Feature Subset'))

outPlot = plot_grid(mzTrimPCA, pcaLoadings, nrow = 2, align = 'vh')
ggplot2::ggsave(filename=file.path('plots', paste0(studyName, '.pdf')), outPlot, width=7, height=8)
ggplot2::ggsave(filename=file.path('plots', paste0(studyName, '.png')), outPlot, width=7, height=8)
