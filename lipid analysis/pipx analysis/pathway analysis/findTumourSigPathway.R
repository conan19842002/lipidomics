library(RMySQL)
library(graph)
library(xlsx)

setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis")
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
getEdgeWeight_1 <- function(reactant, product, tData, nData,alt){
  trData = filterClassPSM_1(tData,reactant)
  tpData = filterClassPSM_1(tData,product)
  nrData = filterClassPSM_1(nData,reactant)
  npData = filterClassPSM_1(nData,product)
  tSumV = list()
  nSumV = list()
  if(is.vector(tpData)){
    tSumV[[1]] = tpData
  }else{
    tSumV[[1]] = rowSums(tpData)  
  }
  if(is.vector(trData)){
    tSumV[[2]] = trData
  }else{
    tSumV[[2]] = rowSums(trData)  
  }
  if(is.vector(npData)){
    nSumV[[1]] = npData
  }else{
    nSumV[[1]] = rowSums(npData)  
  }
  if(is.vector(nrData)){
    nSumV[[2]] = nrData
  }else{
    nSumV[[2]] = rowSums(nrData)  
  }
  ind = which(tSumV[[2]] > 0)
  tDataTest = vector()
  nDataTest = vector()
  tDataTest = tSumV[[1]][ind]/tSumV[[2]][ind]
  ind = which(nSumV[[2]] > 0)                       
  nDataTest = nSumV[[1]][ind]/nSumV[[2]][ind]
  z_score = 0
  out <- tryCatch({
    #if(flag == "active")
    t = t.test(tDataTest,nDataTest, alternative= alt)
    #else
    #  t = t.test(tDataTest,nDataTest)
  },
  error=function(cond){
    message(cond)
    return(0)
  },
  warning=function(cond){
    message(cond)
    return(0)
  }
  
  )
  if(is.finite(t$p.value)){
    z_score = qnorm(1 - t$p.value)
  }  
  z_score
}
savetoFildDiffPathway <- function(tData, nData){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/data/lipid analysis/pipx analysis/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  weight = vector()
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    weight[i] = getEdgeWeight_1(reactant[i], product[i], tData, nData,"different")
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
    subPathway = list(startNode)
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
    finish = FALSE  
    j = 0
    currentNode = nextNode
    z_score = neighborL[1]#getSubPathwayZScore(subPathway, tData, nData)
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
          z_tmp = getSubPathwayZScore(subPathway, tData, nData)#getSubPathwayZScoreByStage(subPathway, tData, nData)
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
    subPathwayL[num] = paste(subPathway, collapse = "->")
    zScoreL[num] = z_score
    num = num + 1
  }
  dt = data.frame(pathway = factor(c(subPathwayL), ordered = TRUE),
                  score = c(zScoreL))
  dt = dt[order(dt$score, decreasing = TRUE),]
  dt = dt[which(dt$score > 1.645),]
  if(length(dt) > 0)
    write.xlsx(dt, "pathway_sig_diff.xlsx")  
}
saveToFileSigPathway <- function(tData, nData, c_pip3,alt){
  # get all reactions
  reactions = as.matrix(read.csv("D:/project/lipidomics/data/lipid analysis/pipx analysis/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  weight = vector()
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    weight[i] = getEdgeWeight_1(reactant[i], product[i], tData, nData,alt)
  }
  if(alt=="greater"){
    atv = "active"
  }else{
    atv="inactive"
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
    subPathway = list(startNode)
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
    finish = FALSE  
    j = 0
    currentNode = nextNode
    z_score = neighborL[1]#getSubPathwayZScore(subPathway, tData, nData)
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
    subPathwayL[num] = paste(subPathway, collapse = "->")
    zScoreL[num] = z_score
    num = num + 1
  }
  dt = data.frame(pathway = factor(c(subPathwayL), ordered = TRUE),
                  score = c(zScoreL))
  dt = dt[order(dt$score, decreasing = TRUE),]
  dt = dt[which(dt$score > 1.645),]
  if(length(dt) > 0){
    #plot pathway graph and save to file
    fileName = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis/",c_pip3,"_tumour_most_",atv,"_pathway",sep="")
    plotPathwayGraph(reactant, product,weight,fileName)
    write.xlsx(dt, paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis/",c_pip3,"_tumour_most_",atv,"_pathway.xlsx",sep=""))  
  }
  
}
