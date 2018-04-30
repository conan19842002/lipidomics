library(RMySQL)
library(graph)
library(xlsx)

setwd("D:/project/lipidomics/data/pathway analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
# get all reactions
reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
# get set of reactants and products
reactant = vector()
product = vector()
weight = vector()
tData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(TRUE))
nData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(FALSE))
tData = getMolarData(tData, TRUE)
nData = getMolarData(nData, TRUE)
# This is for searching changed pathways in different stages
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
for(s in 1:length(stages)){
  tPS = stagePatients[[stages[s]]]
  nPS = nData[patientList[[s]],]
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    weight[i] = getEdgeZScore(reactant[i], product[i], tPS, nPS,"greater")
  }
  # generate a list of found subpathways with its z_score
  subPathwayL = vector()
  zScoreL = vector()
  # build a graph
  df = cbind(reactant,product)
  gr = ftM2graphNEL(df, W= weight, V=NULL, edgemode="directed")
  node = nodes(gr)
  edgeW = edgeWeights(gr)
  visitted = vector()
  num = 1
  for(i in 1:length(node)){
    
    for(k in 1:length(node)){
      visitted[node[k]] = FALSE 
    }
    startNode = node[i]
    visitted[startNode] = TRUE
    subPathway = c(startNode)
    edgeWL = edgeW[[startNode]]
    if(length(edgeWL) == 0){
      next  
    }
    # get list of its neighbors with decreasing order of weights
    neighborL = sort(edgeWL, decreasing = TRUE)
    # get the next node
    nodeL = names(neighborL)
    nextNode = nodeL[1]
    visitted[nextNode] = TRUE
    subPathway  = c(subPathway, nextNode)
    z_score = neighborL[1]#getSubPathwayZScore(subPathway, tData, nData)
    finish = FALSE  
    j = 0
    currentNode = nextNode
    while(!finish){
      nEdgeWL = edgeW[[currentNode]]
      if(length(nEdgeWL) == 0)
        break
      
      neighborL = sort(nEdgeWL, decreasing = TRUE)
      count = length(neighborL)
      notFound = FALSE
      for(j in 1:count){
        nodeL = names(neighborL)
        nextNode = nodeL[j]
        if(!visitted[nextNode]){
          subPathway = c(subPathway, nextNode)
          z_tmp = getSubPathwayZScore(subPathway, tPS, nPS,"greater")#getSubPathwayZScoreByStage(subPathway, tData, nData)
          if(z_tmp > 1.645){
            # choose this node
            currentNode = nextNode
            visitted[nextNode] = TRUE
            z_score = z_tmp
            break  
          }else{
            # remove this node and try to find another
            subPathway = subPathway[-which(subPathway == nextNode)]
            visitted[nextNode] = FALSE
            notFound = TRUE
          }  
        }else{
          finish = TRUE
          break
        }
        
      }
      # if there is no more neighbor to explore
      # then stop growing 
      if(j == count && notFound){
        finish = TRUE
      }
    }
    subPathwayL = c(subPathwayL ,paste(subPathway, collapse = "->"))
    zScoreL = c(zScoreL,z_score)
    
  }
  dt = data.frame(pathway = factor(c(subPathwayL), ordered = TRUE),
                  score = c(zScoreL))
  dt = dt[order(dt$score, decreasing = TRUE),]
  #dt = dt[which(dt$score > 1.645),]
  write.xlsx(dt, paste("most_active_pathway_stage_",stages[s],".xlsx",sep=""))  
}

