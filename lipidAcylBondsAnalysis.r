#this script investigates pd specifically
library(RMySQL)
setwd("D:/project/lipidomics/data/lipid analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
lipidArray = detectionLimits[,1]
tData = getPatientSpecieMatrix(TRUE)
nData = getPatientSpecieMatrix(FALSE)
tData = getMolarData(tData,TRUE)
nData = getMolarData(nData,TRUE)
for(k in 1:length(lipidArray)){
  #tPS = filterClassPSM(applyDetectionLimit(tData),lipidArray[k])
  #nPS = filterClassPSM(applyDetectionLimit(nData),lipidArray[k])
  tPS = filterClassPSM(tData,lipidArray[k])
  nPS = filterClassPSM(nData,lipidArray[k])
  carbonAcylList_N = list()
  doubleBondList_N = list()
  carbonAcylList_T = list()
  doubleBondList_T = list()
  listOfLipid = colnames(tPS)
  listOfCarbonLength = vector()
  listOfDoubleBond = vector()
  
  if(length(listOfLipid) > 0){
    for(i in 1:length(listOfLipid)){
      acyl = unlist(strsplit(unlist(strsplit(listOfLipid[i], split ="-"))[1], split = ":"))
      listOfCarbonLength[i] = acyl[1]
      listOfDoubleBond[i] = acyl[2]
    }
    listOfCarbonLength = unique(listOfCarbonLength)
    listOfDoubleBond = unique(listOfDoubleBond)
    listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
    listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
    # Collect all lipid in normal 
    carbonLengthToRemove = vector()
    doubleBondToRemove = vector()
    tTestCarbonLengthToRemove = vector()
    tTestDoubleBondToRemove = vector()
    for(i in 1:length(listOfCarbonLength)){
      carbonLengthToRemove = vector()
      r = which(grepl(paste(listOfCarbonLength[i],":", sep=""),listOfLipid))
      # sum all the lipid of the same axyl carbon length
      if(length(r) > 1){
        lc1 = rowSums(tPS[,r])
        lc2 = rowSums(nPS[,r])
      }
      else {
        lc1 = tPS[,r]
        lc2 = nPS[,r]
      }
      ind = which(lc1 == 0)
      if(length(ind) > 0)
      for(i1 in 1:length(ind))
        carbonLengthToRemove[length(carbonLengthToRemove)+1] = ind[i1]
      carbonAcylList_T[[listOfCarbonLength[i]]] = lc1
      ind = which(lc2 == 0)
      if(length(ind) > 0)
      for(i1 in 1:length(ind))
        carbonLengthToRemove[length(carbonLengthToRemove)+1] = ind[i1]
      carbonAcylList_N[[listOfCarbonLength[i]]] = lc2 
      
      if(length(carbonLengthToRemove) > 0){
        carbonAcylList_T[[listOfCarbonLength[i]]] = carbonAcylList_T[[listOfCarbonLength[i]]][-carbonLengthToRemove]
        carbonAcylList_N[[listOfCarbonLength[i]]] = carbonAcylList_N[[listOfCarbonLength[i]]][-carbonLengthToRemove]
      }
      if(length(carbonAcylList_T[[listOfCarbonLength[i]]]) < 2 | length(carbonAcylList_N[[listOfCarbonLength[i]]]) < 2){
        tTestCarbonLengthToRemove[length(tTestCarbonLengthToRemove)+1] = i
      }
    }
    if(length(tTestCarbonLengthToRemove) > 0)
      listOfCarbonLength = listOfCarbonLength[-tTestCarbonLengthToRemove]
    for(i in 1:length(listOfDoubleBond)){
      doubleBondToRemove = vector()
      r = which(grepl(paste(":",listOfDoubleBond[i], sep = ""),listOfLipid))
      # sum all the lipid of the same double bond number
      if(length(r) > 1){
        lb1 = rowSums(tPS[,r])  
        lb2 = rowSums(nPS[,r])
      }
      else{
        lb1 = tPS[,r]
        lb2 = nPS[,r]
      }
      
      ind = which(lb1 == 0)
      if(length(ind) > 0)
      for(i1 in 1:length(ind))
        doubleBondToRemove[length(doubleBondToRemove)+1] = ind[i1]
      doubleBondList_T[[listOfDoubleBond[i]]] = lb1
      ind = which(lb2 == 0)
      if(length(ind) > 0)
      for(i1 in 1:length(ind))
        doubleBondToRemove[length(doubleBondToRemove)+1] = ind[i1]
      doubleBondList_N[[listOfDoubleBond[i]]] = lb2 
      
      if(length(doubleBondToRemove) > 0){
        doubleBondList_T[[listOfDoubleBond[i]]] = doubleBondList_T[[listOfDoubleBond[i]]][-doubleBondToRemove]
        doubleBondList_N[[listOfDoubleBond[i]]] = doubleBondList_N[[listOfDoubleBond[i]]][-doubleBondToRemove]
      }
      if(length(doubleBondList_T[[listOfDoubleBond[i]]]) < 2 | length(doubleBondList_N[[listOfDoubleBond[i]]]) < 2){
        tTestDoubleBondToRemove[length(tTestDoubleBondToRemove)+1] = i
      }
    }
    if(length(tTestDoubleBondToRemove) > 0)
      listOfDoubleBond = listOfDoubleBond[-tTestDoubleBondToRemove]
    carbonAcylList = list()
    carbonAcylTestPValue = vector()
    for(i in 1:length(listOfCarbonLength)){
      t = t.test(carbonAcylList_T[[listOfCarbonLength[i]]],carbonAcylList_N[[listOfCarbonLength[i]]], paired = T)
      carbonAcylList[[i]] = log10(carbonAcylList_T[[listOfCarbonLength[i]]]/carbonAcylList_N[[listOfCarbonLength[i]]])
      carbonAcylTestPValue[i] = t$p.value
    }
    tCarbonAcylTestColVec = vector(length = length(carbonAcylTestPValue))
    for(i in 1:length(carbonAcylTestPValue)){
      if(!is.nan(carbonAcylTestPValue[i]) && carbonAcylTestPValue[i] < 0.05){
        tCarbonAcylTestColVec[i] = "red"
      }else{
        tCarbonAcylTestColVec[i] = "white"  
      }
      
    }
    doubleBondList = list()
    doubleBondTestPValue = vector()
    for(i in 1:length(listOfDoubleBond)){
      t = t.test(doubleBondList_T[[listOfDoubleBond[i]]],doubleBondList_N[[listOfDoubleBond[i]]], paired = T)
      doubleBondList[[i]] = log10(doubleBondList_T[[listOfDoubleBond[i]]]/doubleBondList_N[[listOfDoubleBond[i]]])
      doubleBondTestPValue[i] = t$p.value
    }
    tDoubleBondTestColVec = vector(length = length(doubleBondTestPValue))
    for(i in 1:length(doubleBondTestPValue)){
      if(!is.nan(doubleBondTestPValue[i]) && doubleBondTestPValue[i] < 0.05){
        tDoubleBondTestColVec[i] = "red"
      }else{
        tDoubleBondTestColVec[i] = "white"  
      }
    }
    #axis(1, labels = rownames(lipidP)[o], las = 2,cex = 0.5, at = 1:length(orVecL))
    
    windows(width=10, height=8)
    par(mfrow=c(1,2))
    
    pdf(file = paste("Lipid_Acyl_Bonds_Analysis_", paste(lipidArray[k], ".jpg")), width = 9, height = 8, onefile = T)
    boxplot(carbonAcylList,main = lipidArray[k], las = 2, names = listOfCarbonLength , ylab = "Lipid Ratio (T:N)", col = tCarbonAcylTestColVec)
    abline(h=0, col = "blue")
    color = unique(tCarbonAcylTestColVec)
    if(length(color) == 1){
      if(color == "white"){
        legend("bottomright", legend = c("non signif."), fill =c("white"))  
      }else{
        legend("bottomright", legend = c("signif."), fill =c("red"))  
      }
    }else
      legend("bottomright", legend = c("non signif.", "signif."), fill =c("white","red"))  
    #axis(1, labels = listOfCarbonLength, at = 1:length(listOfCarbonLength))
    boxplot(doubleBondList,main = lipidArray[k], las = 2, names = listOfDoubleBond , ylab = "Lipid Ratio (T:N)", col = tDoubleBondTestColVec)
    abline(h=0, col = "blue")
    color = unique(tDoubleBondTestColVec)
    if(length(color) == 1){
      if(color == "white"){
        legend("topleft", legend = c("non signif."), fill =c("white"))  
      }else{
        legend("topleft", legend = c("signif."), fill =c("red"))  
      }
    }else
      legend("topleft", legend = c("non signif.", "signif."), fill = c("white","red"))
    dev.off()
    #axis(1, labels = listOfDoubleBond, at = 1:length(listOfDoubleBond))
  }
  
}
