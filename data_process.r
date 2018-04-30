library(openxlsx)
library(RColorBrewer)
library(NMF)
library(reshape2)
library(qgraph)
library(igraph)
library(graph)
#library(xlsx)
library("factoextra")
require(XLConnect)
require(data.table)
require(ggplot2)
require(gridExtra)
library(randomForest)
library(cluster)
library(PCAmixdata)
library(Boruta)
library(corrplot)
library(dichromat)
source("D:/project/lipidomics/data/lipidomicLib.r")


## process for DRG TRG MRG species
getNewData <- function(data){
  colN = colnames(data)
  colN = gsub("DRG","DG",colN)
  colN = gsub("TRG","TG",colN)
  colN = gsub("MRG","MG",colN)
  colnames(data) = colN
  data
}
tData = getPatientSpecieMatrix(TRUE)
nData = getPatientSpecieMatrix(FALSE)

tData = getNewData(tData)
nData = getNewData(nData)

tData = normalisePatientSpecieMatrix(tData)
nData = normalisePatientSpecieMatrix(nData)
filterClassPSM_1 <- function(data, lipidClass){
  
  lipidSpecies = vector()
  tmp = colnames(data)
  if(lipidClass == "SG"){
    data1 = matrix(data[,"C18 Sphingosine"])
    colnames(data1) = lipidClass
    rownames(data1) = rownames(data)
    data = data1
  }else if(lipidClass == "S1P" | lipidClass == "SPC" | lipidClass == "CH"){
    data1 = matrix(data[,which(grepl(lipidClass, tmp))])
    colnames(data1) = lipidClass
    rownames(data1) = rownames(data)
    data = data1
  }else if(lipidClass == "PIP" | lipidClass == "PIP2" | lipidClass == "PIP3"){
    k = 1
    for(i in 1:length(tmp)){
      # only take species which are formatted as carbon:bond
      if(grepl(":",tmp[i])){
        lipidSpecies[k] = tmp[i]
        k = k + 1
      }
    }
    tSpeciesVec = unlist(strsplit(lipidSpecies, split ="-"))
    l = length(tSpeciesVec)
    data = data[,lipidSpecies]  
    data = data[,which(tSpeciesVec[seq(2,l,2)] == lipidClass)]
  }else{
    data = data[,which(grepl(paste("-",lipidClass,sep=""), tmp))]
  }
  data  
}
xlsx.writeMultipleData <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
}

plotGraphs1 <- function(tData, nData, isAbsolute){
  lipid_class = vector()
  lipid_group = vector()
  concentration = vector()
  lipid_group = vector()
  type = vector()
  phenotTypeL =  c("adenoma", "dukes a", "dukes b", "dukes c", "dukes d")
  lipidClass = c(getLipidClass(),"PIP","PIP2","PIP3")
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  s = vector()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    id = which(names(tData[,1]) %in% rownames(patientData[poi,]))
    if(length(id)>0){
      stagePatients[[stages[i]]] = tData[id,]  
      patientList[[i]] = id
      s = c(s, stages[i])
    }else{
      next
    }
    
  }

  stages = s
  if(isAbsolute){
    yLabel = "Concentration"
  }else{
    yLabel = "Normalised Concentration"
  }
  # extract data for each lipid class
  for(i in 1:length(stages)){
    # extract data for each lipid class
    for(j in 1:length(lipidClass)){
      tPS = stagePatients[[stages[i]]]
      nPS = nData[patientList[[i]],]
      
      dt = as.matrix(filterClassPSM_1(as.matrix(tPS),lipidClass[j]))
      cdt = as.matrix(filterClassPSM_1(as.matrix(nPS),lipidClass[j]))
      # get time point 
      type = c(type, rep("Tumour",nrow(dt)))
      type = c(type, rep("Normal",nrow(dt)))
      # get data for each time point
      # sum up to get total lipid
      dtt = rowSums(dt) 
      cdtt = rowSums(cdt)
      lipid_class = c(lipid_class,rep(lipidClass[j],2*nrow(dt)))
      lipid_group = c(lipid_group, rep(stages[i], 2*nrow(dt)))
      concentration = c(concentration,dtt)  
      concentration = c(concentration,cdtt)  
    }
  }
  
  # create a data frame
  df = data.frame(lipid_group, type,lipid_class, concentration)
  # group lipid concentrations by group, time point and lipid class
  mData = aggregate(df$concentration, by = list(stage = df$lipid_group, type = df$type, lipid_class = df$lipid_class), FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  # need to do a little manipulation since the previous function returned a matrices instead of vectors
  mData <- do.call(data.frame, mData)
  # Now compute standard error for each group
  mData$se <- mData$x.sd / sqrt(mData$x.n)
  
  colnames(mData) <- c("stage", "type", "lipid_class", "mean", "sd", "n", "se")
  colnames(df) <- c("stage", "type", "lipid_class", "concentration")
  for(i in 1:length(stages)){
    wb <- openxlsx::createWorkbook()
    
    for(j in 1:length(lipidClass)){
      df1 = df[df$stage == stages[i] & df$lipid_class == lipidClass[j],]
      df2 = mData[mData$stage == stages[i] & mData$lipid_class == lipidClass[j],]
      
      n = nrow(df1)/2
      dt = data.frame(Tumour = df1[df1$type=="Tumour",]$concentration, Normal = df1[df1$type=="Normal",]$concentration, Tumour.mean = c(df2[df2$type=="Tumour",]$mean, rep("",n-1)), 
                      Tumour.sd = c(df2[df2$type=="Tumour",]$sd, rep("",n-1)), Normal.mean = c(df2[df2$type=="Normal",]$mean, rep("",n-1)), 
                      Normal.sd = c(df2[df2$type=="Normal",]$sd, rep("",n-1)))
      addWorksheet(wb, lipidClass[j])
      writeData(wb, sheet = lipidClass[j], x = dt)
      #write.xlsx(dt, file=paste("D:/project/lipidomics/data/lipid analysis/",stages[i],".xlsx", sep=""), sheetName = lipidClass[j], append=TRUE)
      
    }
    saveWorkbook(wb, file = paste("D:/project/lipidomics/data/lipid analysis/",stages[i],".xlsx", sep=""),overwrite=F)
  }
  
}