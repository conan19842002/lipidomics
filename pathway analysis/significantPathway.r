library(RMySQL)
setwd("D:/project/lipidomics/data/pathway analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
# all the nodes of the pathway are sub-species of lipid classes
# therefore notice TRG, DRG not TG, DG
pathway = c("DHCer", "Cer")
lipidClass = c("DHCer", "Cer")
# get set of subpathways of the different layers
subpathways = list()
# get patient-species matrix from normal and tumour
tData = normalisePatientSpecieMatrix(applyDetectionLimit(getPatientSpecieMatrix(TRUE)))
nData = normalisePatientSpecieMatrix(applyDetectionLimit(getPatientSpecieMatrix(FALSE)))
iSpeciesList = list()
# get set of lipid species which involve in the pathway
for(i in 1:length(lipidClass)){
  # get all sub-species from the tumour
  tSpecies = colnames(filterClassPSM(tData,lipidClass[i]))
  # extract carbon:bond from the sub-species, then store in lists
  iSpeciesList[[pathway[i]]] = getFAs(tSpecies)
}
# find the sub-species of the same layer
iSharedList = Reduce(intersect, iSpeciesList)
if(length(iSharedList) == 0) return

# generate subpathways of different layers
z_score = vector()
# K is the number of edges
K = length(pathway)-1
zColor = vector()
tSubPathways = list()
nSubPathways = list()
for(i in 1:length(iSharedList)){
  for(j in 1:length(pathway)){
    species = paste(paste(iSharedList[i], "-", sep = ""), pathway[j], sep = "")
    # first get data of the lipid class
    tPS = filterClassPSM(tData,lipidClass[j])
    nPS = filterClassPSM(nData,lipidClass[j])
    tSubPathways[[species]] = tPS[,which(colnames(tPS) == species)]
    nSubPathways[[species]] = nPS[,which(colnames(nPS) == species)]  
  }
}
tSubPathWayNode = list()
nSubPathWayNode = list()
size = length(pathway)-1
for(i in 1:length(iSharedList)){
  zc = vector()
  for (j in 1:size){
    reactant = paste(paste(iSharedList[i], "-", sep = ""), pathway[j], sep = "")
    product = paste(paste(iSharedList[i], "-", sep = ""), pathway[j+1], sep = "")
    tSubPathWayNode[[j]] = tSubPathways[[product]]/ tSubPathways[[reactant]]
    nSubPathWayNode[[j]] = nSubPathways[[product]]/ nSubPathways[[reactant]]
    # perform t-test and get the p-value, then convert it to z-score
    t = t.test(tSubPathWayNode[[j]],nSubPathWayNode[[j]], alternative= "greater", paired = T)
    if(is.finite(t$p.value)){
      zc[j] = qnorm(1 - t$p.value)
    }else{
      zc[j] = 0
    }  
  }
  # Compute z score for each subpathway
  z_score[iSharedList[i]] = sum(zc)/sqrt(K)
  # if z score is greater than the standard z score which is 1.645 then the subpathway is significant
  if(abs(z_score[iSharedList[i]]) > 1.645){
    zColor[iSharedList[i]] = "red" # significant
  }else{
    zColor[iSharedList[i]] = "black" # unsignificant
  }
}
# plot layer verse z score
title = pathway[1]
for(i in 2:length(pathway)){
  title = paste(title,pathway[i], sep="-")
}
plot(z_score, pch=20,  main = title, col = zColor, xaxt="n", xlab="Layer", ylab="Z-Score")
axis(1, at=1:length(iSharedList), labels=iSharedList)
color = unique(zColor)
if(length(color) == 1){
  if(color == "white"){
    legend("bottomleft", legend = c("non active"), fill =c("black"))  
  }else{
    legend("bottomleft", legend = c("active"), fill =c("red"))  
  }
}else
legend("topleft", legend = c("non active", "active"), fill = c("black","red"))

