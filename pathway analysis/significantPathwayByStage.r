library(RMySQL)
setwd("D:/project/lipidomics/data/pathway analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
# all the nodes of the pathway are sub-species of lipid classes
# therefore notice TRG, DRG not TG, DG
pathway = c("SM", "Cer")
lipidClass = c("SM", "Cer")
# get set of subpathways of the different layers
subpathways = list()
# get patient-species matrix from normal and tumour
tData = normalisePatientSpecieMatrix((getPatientSpecieMatrix(TRUE)))
nData = normalisePatientSpecieMatrix((getPatientSpecieMatrix(FALSE)))
# get stages of patients
patientData = getPatientData()
pml = stageMeta(patientData)
patientMeta = pml$patientMeta
#split into stages
stages = unique(patientMeta)

stagePatients =  list()
patientList = list();
for(i in 1:length(stages)){
  poi = which(patientMeta == stages[i])  
  patientList[[i]] = poi
  stagePatients[[stages[i]]] = tData[poi,]
}
pathwayN = paste(pathway,collapse = "->")
for(s in 1:length(stages)){
  tPS = stagePatients[[stages[s]]]
  nPS = nData[patientList[[s]],]
  iSpeciesList = list()
  # get set of lipid species which involve in the pathway
  for(i in 1:length(lipidClass)){
    # get all sub-species from the tumour
    tSpecies = colnames(filterClassPSM(tPS,lipidClass[i]))
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
      tPS1 = filterClassPSM(tPS,lipidClass[j])
      nPS1 = filterClassPSM(nPS,lipidClass[j])
      tSubPathways[[species]] = tPS1[,which(colnames(tPS1) == species)]
      nSubPathways[[species]] = nPS1[,which(colnames(nPS1) == species)]
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
      tSubPathWayTest = tSubPathways[[product]]/ tSubPathways[[reactant]]
      nSubPathWayTest = nSubPathways[[product]]/ nSubPathways[[reactant]]
      tSubPathWayTest[is.infinite(tSubPathWayTest)] = .Machine$double.xmax
      nSubPathWayTest[is.infinite(nSubPathWayTest)] = .Machine$double.xmax
      
      # perform t-test and get the p-value, then convert it to z-score
      t = t.test(tSubPathWayTest,nSubPathWayTest, alternative= "greater", paired = T)
      if(is.finite(t$p.value)){
        zc[j] = qnorm(1 - t$p.value)
      }else{
        zc[j] = 0
      }  
    }
    # Compute z score for each subpathway
    z_score[iSharedList[i]] = sum(zc)/sqrt(K)
  }
  dt = data.frame(pathwayN = factor(c(iSharedList), ordered = TRUE),
                  score = c(z_score))
  colnames(dt) = c(pathwayN, "score")
  dt = dt[order(dt$score, decreasing = TRUE),]
  write.xlsx(dt,paste(paste(pathway,collapse="_"),".xlsx",sep=""), sheetName= stages[s], append= TRUE)  
}

