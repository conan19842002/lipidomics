library(RMySQL)
library(gplots)
library(RColorBrewer)
source("D:/project/lipidomics/data/lipidomicLib.r")
require(XLConnect)
setwd("D:/project/lipidomics/data/lipid analysis")
# global variables

getPatientVec <- function(start,end, df, patientVec){
  for(r in start:end){
    patientNum = unlist(strsplit(df[r,8], split = " "))[1]
    sampleNum = df[r,7]
    patientVec[sampleNum] = patientNum
  }
  patientVec
}
doPIPAnalysis <- function(){
  wb = loadWorkbook("20130716_Colorectal 2-26 PIPx.xlsx")
  df = readWorksheet(wb, sheet = 3, header = TRUE)
  patientVec = vector()
  patientVec = getPatientVec(10,23,df,patientVec)
  patientVec = getPatientVec(25,33,df,patientVec)
  patientVec = getPatientVec(35,40,df,patientVec)
  patientVec = getPatientVec(43,44,df,patientVec)
  patientVec = getPatientVec(47,60,df,patientVec)
  patientVec = getPatientVec(62,79,df,patientVec)
  patientVec = getPatientVec(82,87,df,patientVec)
  patientVec = getPatientVec(90,95,df,patientVec)
  patientVec = getPatientVec(98,117,df,patientVec)
  patientVec = getPatientVec(120,125,df,patientVec)
  patientVec = getPatientVec(128,135,df,patientVec)
  patientVec = getPatientVec(138,148,df,patientVec)
  patientNum = length(patientVec)/2
  tData = matrix(nrow = 15, ncol = patientNum)
  nData = matrix(nrow = 15, ncol = patientNum)
  colnames(tData) = patientVec[seq(1,length(patientVec),2)]
  colnames(nData) = patientVec[seq(2,length(patientVec),2)]
  
  df = readWorksheet(wb, sheet = 2, header = TRUE)
  rownames(tData) = df[29:43,1]
  rownames(nData) = df[29:43,1]
  i = 1
  for(colNum in seq(3,16,2)){
    tData[1:15,i] = df[29:43,colNum]
    nData[1:15,i] = df[29:43,colNum+1]
    i = i + 1
  }
  for(colNum in seq(19,26,2)){
    tData[1:15,i] = df[29:43,colNum]
    nData[1:15,i] = df[29:43,colNum+1]
    i = i + 1
  }
  
  # process second file
  wb = loadWorkbook("20130717_Colorectal 27-52 PIPx.xlsx")
  df = readWorksheet(wb, sheet = 2, header = TRUE)
  for(colNum in seq(1,22,2)){
    tData[1:15,i] = df[2:16,colNum]
    nData[1:15,i] = df[2:16,colNum+1]
    i = i + 1
  }
  # process third file
  wb = loadWorkbook("20130717_Colorectal 53-77 PIPx.xlsx")
  df = readWorksheet(wb, sheet = 2, header = TRUE)
  i = ncols(tData)+1
  for(colNum in seq(1,23,2)){
    tData[1:15,i] = df[2:16,colNum]
    nData[1:15,i] = df[2:16,colNum+1]
    i = i + 1
  }
  # process fourth file
  wb = loadWorkbook("20130718_Colorectal 78-103 PIPx.xlsx")
  df = readWorksheet(wb, sheet = 2, header = TRUE)
  for(colNum in seq(1,22,2)){
    tData[1:15,i] = df[2:16,colNum]
    nData[1:15,i] = df[2:16,colNum+1]
    i = i + 1
  }
  # process fifth file
  wb = loadWorkbook("20130719_Colorectal 104-143 PIPx.xlsx")
  df = readWorksheet(wb, sheet = 2, header = TRUE)
  for(colNum in seq(1,30,2)){
    tData[1:15,i] = df[2:16,colNum]
    nData[1:15,i] = df[2:16,colNum+1]
    i = i + 1
  }
  
  tData = t(tData)
  nData = t(nData)
}
