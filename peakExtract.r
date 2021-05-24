# Extract MS peaks from .raw data
# Procedure:
#   Select m/z values and corresponding intensities from MS1 scans
#   Create m/z and rt range based on expected ppm mass tolerance and peak widths
#
# Fulgens Consulting 

library(readr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(reshape2,quietly=TRUE)

studyID <- c('ST000934')
exclFiles <- c('S-221.raw')
filePath = file.path('src')
source(file.path(filePath, 'setup_msprocess.R'))
rawFiles <- list.files(file.path('data/raw', studyID))[c(grep('.raw', list.files(file.path('data/raw', studyID))))]

for(rF in rawFiles){
    if(!rF %in% exclFiles & !gsub('.raw', '.csv', rF) %in% list.files(file.path('data/processed', studyID))){
        collectDF <- NULL
        header <- rawrr::readFileHeader(file.path('data/raw', studyID, rF))
        for(i in seq(1000,2000,200)){
            subList <- buildScanMZTbl(rawFile=file.path('data/raw', studyID, rF), 
                                        dropMS2=TRUE, scan=i:(i+199), minMZ=300, maxMZ=600)
            subList <- addPpmRTWindows(masterList=subList, ppm=5, RTWindow=10)
            subList$sample = gsub('.raw', '', rF)
            if(i==1000){
                write_csv(subList, paste0('data/processed/', studyID, '/', gsub('.raw', '.csv', rF)))
            } else {
                write_csv(subList, paste0('data/processed/', studyID, '/', gsub('.raw', '.csv', rF)), append=TRUE)
            }
            rm(subList)
            gc()
        }
    }
}
