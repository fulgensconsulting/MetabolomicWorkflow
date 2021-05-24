# Setup for feature processing, statistical analysis on high-res LC-MS data 
#
# Fulgens Consulting 

library(readr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(reshape2, quietly=TRUE)
library(rawrr, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(ggrepel, quietly=TRUE)
library(glmnet, quietly=TRUE)

#' Build sample, scan, mz, retention time, intensity table from .raw
#' @param rawFile: Str - Raw file
#' @param dropMS2: Bool - Drop MS2 scans
#' @param scan: List - Optionally filter explicit scan numbers
#' @param minRT: Float - Optionally filter starting retention time (seconds)
#' @param maxRT: Float - Optionally filter maximum retention time (seconds)
#' @param minMZ: Float - Optionally filter starting mZ
#' @param maxMZ: Float - Optionally filter maximum mZ
#' @return mapTable with added scan and sample information
buildScanMZTbl = function(rawFile, dropMS2=TRUE, scan=NULL, minRT=30, maxRT=NULL, 
    minMZ=NULL, maxMZ=NULL, ...){

    header <- rawrr::readFileHeader(rawFile)
    scanLen <- header$`Number of scans`
    if(is.null(scan)){
        scan <- rawrr::readIndex(rawFile)
        if(dropMS2 == TRUE){
            scan <- scan[which(scan$MSOrder == 'Ms'), ][[1]]
        }
        mzRaw=rawrr::readSpectrum(rawFile, scan=scan)
    } else {
        if(dropMS2 == TRUE){
            scanRead <- rawrr::readIndex(rawFile)
            scanRead <- scanRead[scan, ]
            scan <- scanRead[which(scanRead$MSOrder == 'Ms'), ][[1]]
            mzRaw <- rawrr::readSpectrum(rawFile, scan=scan)
        } else {mzRaw=rawrr::readSpectrum(rawFile, scan=scan)}
    }

    masterList = NULL    
    for(sc in 1:length(scan)){
        
        scanInfo = mzRaw[[sc]]
        if("noises" %in% names(scanInfo)){
            scList = as.data.frame(cbind(scanInfo$mZ, scanInfo$intensity, scanInfo$noises))
            scList$intensity = scList$intensity - scList$noises
            scList = scList %>% select(mZ, intensity) %>% as.data.frame(scList)
        } else {
            scList = as.data.frame(cbind(scanInfo$mZ, scanInfo$intensity))
        }
        colnames(scList) = c('mz', 'intensity')
        scList$scan = scanInfo$scan
        scList$rt = scanInfo$rtinseconds
        masterList = rbind(masterList, scList)
    }

    if(header$`Sample name` == ""){
        masterList$sample = header$`Sample id`
    } else {masterList$sample = header$`Sample name`}

    if(is.numeric(minRT)){
        masterList = masterList[masterList$rt >= minRT, ]
        if(nrow(masterList) == 0){
            stop(paste('Nothing left after minimum RT filter of ', minRT))
        }
    }
    if(is.numeric(maxRT)){
        masterList = masterList[masterList$rt <= maxRT, ]
        if(nrow(masterList) == 0){
            stop(paste('Nothing left after maximum RT filter of ', minRT))
        } 
    }
    if(is.numeric(minMZ)){
        masterList = masterList[masterList$mz >= minMZ, ]
        if(nrow(masterList) == 0){
            stop(paste('Nothing left after minimum MZ filter of ', minMZ))
        }
    }
    if(is.numeric(maxMZ)){
        masterList = masterList[masterList$rt <= maxMZ, ]
        if(nrow(masterList) == 0){
            stop(paste('Nothing left after maximum MZ filter of ', maxMZ))
        } 
    }

    return(masterList)
}

#Add ppm and RTwindows to the mz:RT pairs. This saves time by calculating outside loops
#' @param masterList: data table containing m/z's and RTs
#' @param ppm : ppm tolerance
#' @param RTWindow: retention time window
#' @return masterList with added windows of ppm and RT
addPpmRTWindows = function(masterList, ppm, RTWindow, dropMZRT=TRUE, ...){
    masterList$massLower = ppmWindowLower(masterList$mz, ppm = ppm)
    masterList$massUpper = ppmWindowUpper(masterList$mz, ppm = ppm)
    masterList$rtLower = masterList$rt - RTWindow/2
    masterList$rtUpper = masterList$rt + RTWindow/2

    if(dropMZRT == TRUE){
        masterList <- masterList %>% select(sample, intensity, massLower, massUpper, rtLower, rtUpper)
    }
    return(masterList)
}

#Find lower mass with given ppm and mz
ppmWindowLower = function(mz, ppm){
    return(abs((mz*ppm/1e06) - mz))
}
#Find upper mass with given ppm and mz
ppmWindowUpper = function(mz, ppm){
    return((mz*ppm/1e06) + mz)
}

#Build ion counts from the buildScanMZTbl() object and RT windows from addPpmRTWindows()
#' @param masterList: data table containing m/z's and RTs
#' @return data table containing AUCs for discrete mass and retention bins
massRTCompile = function(masterList){

    fillDF = as.data.frame(cbind(masterList, 0)) %>% rename(auc = "0")  
    
    for(i in 1:nrow(fillDF)){
        massLowMatch = fillDF$massLower[i]
        massUpMatch = fillDF$massUpper[i]
        RTLowMatch = fillDF$rtLower[i]
        RTUpMatch = fillDF$rtUpper[i]

        fillDF$auc[i] = sum(masterList$intensity[(masterList$mz >= massLowMatch) &
                                                (masterList$mz <= massUpMatch) &
                                                (masterList$rt >= RTLowMatch) &
                                                (masterList$rt <= RTUpMatch)])
    
    }
    return(fillDF)
}

#' Fit PCA on aggregated MZ data
#'  Note: columns are removed which have 0 variance (either the MZ was 0 in all samples,
#'      or the intensity was sample across samples for which there was a real ion count,
#'      which would be odd..)
#' @param data: Individual MZ data aggregated into NxM
#' @param nonDatCols: List - Non-data-containing columns
#' @param returnFit: List - Return only the fit object
#' @return PCA scores
fitPCA = function(data, nonDatCols=c('Sample'), returnFit=FALSE){
    dataOnly = data %>% select(-nonDatCols)
    dataOnly[is.na(dataOnly)] = 0
    dataOnly = dataOnly[ , apply(dataOnly, 2, var) != 0]
    pcaFit = prcomp(dataOnly, center = TRUE, scale = TRUE)

    if(returnFit == TRUE){
        return(pcaFit)
    } else {
        pcaDF = as.data.frame(pcaFit$x)
        projections = cbind(data$Sample, pcaDF[, c('PC1' ,'PC2')])
        colnames(projections)[1] = 'Sample'
        return(projections)
    }
}

#' Calculate the ellipse for a confidence interval on the PCA projections
#' @param pcaProjections: Scores from a prcomp() fit
#' @param confidence: Hotellings confidence (usually 95%)
#' @return 2D Ellipse boundaries for plotting the confidence interval
calcPCAEllipse = function(pcaProjections, confidence){
    theta = c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle = as.matrix(cbind(cos(theta), sin(theta)))
    sigma = var(cbind(pcaProjections$PC1, pcaProjections$PC2))
    mu = c(mean(pcaProjections$PC1), mean(pcaProjections$PC2)) #should be centered already
    ed = sqrt(qchisq(confidence/100, df = 2))
    ell = data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                    mu, FUN = "+"), groups = 1)
    names(ell)[1:2] <- c("xvar", "yvar")
    return(ell)
}


#' Plot PCA for the aggregated MZ data
#' @param data: Individual MZ data aggregated into NxM
#' @param confidence: Hotellings confidence (usually 95%)
#' @param minVals: Float - Keep features with at least minVals percentage of non-null values
#' @return ggplot object of PCA scores
plotPCA = function(data, confidence, nonDatCols=c('Sample'), minVals=0.5,
    ggTitle=NULL){

    if(is.numeric(minVals)){
        minSum = minVals*nrow(data)
        dataOnly = data %>% select(-nonDatCols)
        dataOnly = dataOnly[, sapply(dataOnly, function(x) sum(is.na(x))) < minSum] 
        if(nrow(dataOnly) == 0){
            stop(paste('Nothing left after minimum non null value filter of ', minVals))
        }
        data = cbind(data$Sample, dataOnly)
        colnames(data)[1] = 'Sample'
    }
    
    pcaProjections = fitPCA(data, nonDatCols=nonDatCols)
    pcaEllipse = calcPCAEllipse(pcaProjections, confidence)

    #Re-run part of fitPCA just to get the variance explained
    dataOnly = data %>% select(-nonDatCols)
    dataOnly[is.na(dataOnly)] = 0
    dataOnly = dataOnly[ , apply(dataOnly, 2, var) != 0]
    pcaFit = prcomp(dataOnly, center = TRUE, scale = TRUE)
    varExp1 = round(summary(pcaFit)$importance[2,1] * 100, 2) 
    varExp2 = round(summary(pcaFit)$importance[2,2] * 100, 2)

    pcaPlot = ggplot(pcaProjections, aes(x = PC1, y = PC2, colour = Sample))+
        geom_point(size = 5)+
        geom_path(data = pcaEllipse, aes(x = xvar, y = yvar), color = 'black')+
        theme_bw()+
        ylab(sprintf('PC2 (%s%s variance explained)', varExp2, '%')) + 
        xlab(sprintf('PC1 (%s%s variance explained)', varExp1, '%')) + 
        ggtitle(ggTitle)+
        geom_abline(aes(intercept = 0, slope = 0), color = 'black', size = 0.25)+
        geom_vline(aes(xintercept = 0), color = 'black', size = 0.25)+
        theme(axis.text = element_blank(), axis.title = element_text(size = 11, face = "bold"), 
                legend.title = element_text(size = 11), legend.text = element_text(size = 10), 
                title = element_text(size = 12, face = "bold"),
                plot.title = element_text(hjust = 0.5))

    return(pcaPlot)
}


#Create an elastic net regression procedure to take "most important" features
#... then see how the trimmed set does on PCA. Assume "Sample" is main Y feature,
#... which may represent group or individual sample labels.
#' @param data: AUCs by Sample.
#' @param alpha: Elastic-net mixing parameter, where alpha = 0 is Ridge regression
#'              and alpha = 1 is LASSO regression
#' @return Subset of features from data which exhibit strongest regressions in data
ElasticNetVariableSelection = function(data, alpha=0.5, minVals=0.5){

    if(is.numeric(minVals)){
        minSum = minVals*nrow(data)
        dataOnly <- data %>% select(-Sample)
        dataOnly = dataOnly[, sapply(dataOnly, function(x) sum(is.na(x))) < minSum] 
        if(nrow(dataOnly) == 0){
            stop(paste('Nothing left after minimum non null value filter of ', minVals))
        }
    }

    dataOnly[is.na(dataOnly)] <- 0
    dataOnly$Sample <- as.numeric(as.factor(data$Sample))
    elasticX = model.matrix(Sample ~., dataOnly)
    elasticY = dataOnly$Sample
    elasticCV = cv.glmnet(elasticX, elasticY, alpha = alpha, type.measure = 'mse', 
        standardize = TRUE)
    lambda_opt = elasticCV$lambda.1se
    coefDF = as.matrix(coef(elasticCV, s = lambda_opt))
    subFeatures = names(coefDF[(coefDF[, 1] > 0), ]) #Only non-zero coefficients
    subFeatures = subFeatures[-c(grep('Intercept', subFeatures))]
    subFeatures = gsub('\\`', '', subFeatures)
    return(subFeatures)
}


#' Plot PCA Loadings from prcomp() fit
#' @param origDF: Original data table, which contains protein names
#' @param pcaFitObject: prcomp() fit object
#' @param varNames: Add the varNames to the text
#' @param adjustScale: Adjust axis scale based on scores values?
#' @return ggplot object of PCA loadings
PCALoadingsPlot = function(origDF, pcaFitObject, varNames = TRUE, adjustScale = TRUE,
    ggTitle=NULL){

    loadingPlot = as.data.frame(cbind(rownames(pcaFitObject$rotation), 
                pcaFitObject$rotation[, c('PC1', 'PC2')])) %>% rename(Variable=V1)
    loadingPlot$PC1 = as.numeric(as.character(unlist(loadingPlot$PC1)))
    loadingPlot$PC2 = as.numeric(as.character(unlist(loadingPlot$PC2)))

    loadingPlot = merge(loadingPlot, origDF, by.x = c('Variable'), by.y = c('mz_rt')) %>% select(Variable, PC1, PC2)
    loadingPlot <- unique(loadingPlot)

    loadingGG = ggplot(loadingPlot, aes(x = PC1, y = PC2))+
            geom_point(size = 4)+
            theme_bw()+
            ggtitle(ggTitle)+
            geom_hline(yintercept = 0, alpha = 0.5)+
            geom_vline(xintercept = 0, alpha = 0.5)+
            theme(axis.text = element_blank(), axis.title = element_text(size = 11, face = "bold"), 
                legend.title = element_text(size = 11), legend.text = element_text(size = 10), 
                title = element_text(size = 10, face = "bold"),
                plot.title = element_text(hjust = 0.5))
    
    if(varNames == TRUE){
        loadingGG = loadingGG + 
            geom_text_repel(aes(label = Variable), size = 4, max.iter = 10000)
    }
    
    if(adjustScale == TRUE){
        loadingGG = loadingGG + 
            xlim(c(max(pcaFit$x[, c('PC1')]), min(pcaFit$x[, c('PC1')]))) +
            ylim(c(max(pcaFit$x[, c('PC2')]), min(pcaFit$x[, c('PC2')])))
    }

    return(loadingGG)
}
