library(RMySQL)
library(graph)
library(xlsx)
require(gridExtra)
library(NMF)
library(qgraph)
setwd("D:/project/lipidomics/data/pathway analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
## process for DRG TRG MRG species
getNewData <- function(data){
  colN = colnames(data)
  colN = gsub("DRG","DG",colN)
  colN = gsub("TRG","TG",colN)
  colN = gsub("MRG","MG",colN)
  colN = gsub("C18 Sphingosine","C18-SG",colN)
  colnames(data) = colN
  data
}
getSigPathways <- function(){
  # this data does not include PIPx
  
  tData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(TRUE))
  nData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(FALSE))
  tData = getMolarData(tData, TRUE)
  nData = getMolarData(nData, TRUE)
  tData = getNewData(tData)
  nData = getNewData(nData)
  # find the most active pathways
  findMostSigPathways(tData,nData,"greater")
  # find the most inactive pathways
  findMostSigPathways(tData,nData,"less")
  # find all active pathways
  findAllSigPathways(tData, nData,"greater")
  # find all inactive pathways
  findAllSigPathways(tData, nData,"less")
  #find the most active pathways by stages
  findMostSigPathwaysByStage(tData, nData, "greater")
  #find the most inactive pathways by stage
  findMostSigPathwaysByStage(tData, nData, "less")
  #find all acive pathways by stage
  findAllSigPathwaysByStage(tData, nData, "greater")
  #find all inactive pathways by stage
  findAllSigPathwaysByStage(tData, nData, "less")
}

findAllSigPathwaysByStage <- function(tData, nData, alt){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  ZScore = vector()
  
  if(alt=="greater"){
    atv = "active"
  }else{
    atv = "inactive"
  }
  # This is for searching changed pathways in different stages
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  
  
  stages = unique(patientMeta)
  s = vector()
  stagePatients =  list()
  patientList = list()
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
      ZScore[i] = getEdgeZScore(reactant[i], product[i], tPS, nPS,alt)
    }
    id = which(ZScore==0)
    if(length(id)>0){
      ZScore = ZScore[-id]
      reactant = reactant[-id]
      product = product[-id]
    }
    # generate a list of found subpathways with its z_score
    subPathwayL = vector()
    zScoreL = vector()
    # build a graph
    df = cbind(reactant,product)
    gr = ftM2graphNEL(df, W= ZScore, V=NULL, edgemode="directed")
    node = nodes(gr)
    edgeW = edgeWeights(gr)
    visited = vector()
    for(i in 1:length(node)){
      visited[node[i]] = FALSE  
    }
    dataL = list()
    dataL[[1]] = vector()
    for(i in 1:length(node)){
      visited[which(visited==TRUE)] = FALSE
      dataL[[2]] = visited
      if(!visited[node[i]]){
        pathway = c(node[i])
        dataL = getSigPathway(dataL,pathway, edgeW, node[i], tPS, nPS,alt,TRUE)
      }
    }
    dt = data.frame(pathway = factor(names(dataL[[1]])), score = unname(dataL[[1]]))
    dt = dt[order(dt$score, decreasing = TRUE),]
    write.xlsx(dt, paste("D:/project/lipidomics/data/pathway analysis/most_",atv,"_pathway_stage_",stages[s],".xlsx",sep=""))  
    #plot pathway graph and save to file
    fileName = paste("D:/project/lipidomics/data/pathway analysis/most_",atv,"_pathway_stage_",stages[s],sep="")
    plotPathwayGraph(reactant, product,ZScore,fileName)
  }
}
findMostSigPathwaysByStage <- function(tData, nData, alt){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  zScore = vector()
  if(alt=="greater"){
    atv = "active"
  }else{
    atv = "inactive"
  }
  # This is for searching changed pathways in different stages
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  
  stagePatients =  list()
  patientList = list()
  
    
  stages = unique(patientMeta)
  s = vector()
  stagePatients =  list()
  patientList = list()
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
  for(s in 1:length(stages)){
    tPS = stagePatients[[stages[s]]]
    nPS = nData[patientList[[s]],]
    for(i in 1:nrow(reactions)){
      reactant[i] = reactions[i,1]
      product[i] = reactions[i,2]  
      zScore[i] = getEdgeZScore(reactant[i], product[i], tPS, nPS,alt)
    }
    id = which(zScore==0)
    if(length(id)>0){
      zScore = zScore[-id]
      reactant = reactant[-id]
      product = product[-id]
    }
    # generate a list of found subpathways with its z_score
    subPathwayL = vector()
    zScoreL = vector()
    # build a graph
    df = cbind(reactant,product)
    gr = ftM2graphNEL(df, W= zScore, V=NULL, edgemode="directed")
    node = nodes(gr)
    edgeW = edgeWeights(gr)
    visitted = vector()
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
      if(z_score < 1.645)
        next
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
            z_tmp = getSubPathwayZScore(subPathway, tPS, nPS,alt)#getSubPathwayZScoreByStage(subPathway, tData, nData)
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
                    score = c(unname(zScoreL)))
    dt = dt[order(dt$score, decreasing = TRUE),]
    write.xlsx(dt, paste("D:/project/lipidomics/data/pathway analysis/most_",atv,"_pathway_stage_",stages[s],".xlsx",sep=""))  
    #plot pathway graph and save to file
    fileName = paste("D:/project/lipidomics/data/pathway analysis/most_",atv,"_pathway_stage_",stages[s],sep="")
    
    
    plotPathwayGraph(reactant, product,zScore,fileName)
  }
}
findAllSigPathways <- function(tData, nData, alt){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  ZScore = vector()
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    ZScore[i] = getEdgeZScore(reactant[i], product[i], tData, nData,alt)
  }
  id = which(ZScore==0)
  if(length(id)>0){
    ZScore = ZScore[-id]
    reactant = reactant[-id]
    product = product[-id]
  }
  # generate a list of found subpathways with its z_score
  subPathwayL = vector()
  zScoreL = vector()
  # build a graph
  df = cbind(reactant,product)
  gr = ftM2graphNEL(df, W= ZScore, V=NULL, edgemode="directed")
  node = nodes(gr)
  edgeW = edgeWeights(gr)
  visited = vector()
  for(i in 1:length(node)){
    visited[node[i]] = FALSE  
  }
  dataL = list()
  dataL[[1]] = vector()
  for(i in 1:length(node)){
    visited[which(visited==TRUE)] = FALSE
    dataL[[2]] = visited
    if(!visited[node[i]]){
      pathway = c(node[i])
      dataL = getSigPathway(dataL,pathway, edgeW, node[i], tData, nData,alt,TRUE)
    }
  }
  dt = data.frame(pathway = factor(c(names(dataL[[1]])), ordered=TRUE), score = c(unname(dataL[[1]])))
  dt = dt[order(dt$score, decreasing = TRUE),]
  if(alt=="greater"){
    fileName = "D:/project/lipidomics/data/pathway analysis/all_active_pathway"
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/all_active_pathway.xlsx")
  }else{
    fileName="D:/project/lipidomics/data/pathway analysis/all_inactive_pathway"
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/all_inactive_pathway.xlsx")
  }
  #plot pathway graph and save to file
  plotPathwayGraph(reactant, product,ZScore,fileName)  
}
findMostSigPathways <- function(tData,nData,alt){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction2s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  ZScore = vector()
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    ZScore[i] = getEdgeZScore(reactant[i], product[i], tData, nData,alt)
  }
  id = which(ZScore==0)
  if(length(id)>0){
    ZScore = ZScore[-id]
    reactant = reactant[-id]
    product = product[-id]
  }
  # generate a list of found subpathways with its z_score
  subPathwayL = vector()
  zScoreL = vector()
  # build a graph
  df = cbind(reactant,product)
  gr = ftM2graphNEL(df, W= ZScore, V=NULL, edgemode="directed")
  node = nodes(gr)
  edgeW = edgeWeights(gr)
  visitted = vector()
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
    # get list of its neighbors with decreasing order of ZScores
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
    if(z_score < 1.645)
      next
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
          z_tmp = getSubPathwayZScore(subPathway, tData, nData,alt)#getSubPathwayZScoreByStage(subPathway, tData, nData)
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
    subPathwayL = c(subPathwayL,paste(subPathway, collapse = "->"))
    zScoreL = c(zScoreL,z_score)
  }
  dt = data.frame(pathway = factor(c(subPathwayL), ordered = TRUE),
                  score = c(unname(zScoreL)))
  dt = dt[order(dt$score, decreasing = TRUE),]
  if(alt=="greater"){
    fileName = "D:/project/lipidomics/data/pathway analysis/most_active_pathway"
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/most_active_pathway.xlsx")
  }else{
    fileName = "D:/project/lipidomics/data/pathway analysis/most_inactive_pathway"
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/most_inactive_pathway.xlsx")
  }
  plotPathwayGraph(reactant, product, ZScore, fileName)  
}
