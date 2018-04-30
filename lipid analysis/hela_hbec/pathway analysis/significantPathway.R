library(RMySQL)
library(graph)
library(xlsx)
require(XLConnect)
setwd("D:/project/lipidomics/data/lipid analysis/hela_hbec/pathway analysis")
source("D:/project/lipidomics/data/lipidomicLib.r")
reactions = as.matrix(read.csv("D:/project/lipidomics/data/lipid analysis/hela_hbec/reaction1s.csv", header = T))
findSigPathway <- function(){
  wb = loadWorkbook("D:/project/lipidomics/data/lipid analysis/hela_hbec/hela_hbec.xlsx")
  sheets = getSheets(wb)
  excepTL = c("Cer","SM")
  unInfectedHela = list()
  HRV1BHela = list()
  HRV16Hela = list()
  unInfectedHbec = list()
  HRV1BHbec = list()
  HRV16Hbec = list()
  lUnInfectedHela = matrix(nrow = 2)
  lHRV1BHela = matrix(nrow = 2)
  lHRV16Hela = matrix(nrow = 2)
  
  
  lUnInfectedHbec = matrix(nrow = 2)
  lHRV1BHbec = matrix(nrow = 2)
  lHRV16Hbec = matrix(nrow = 2)
  
  l = 0
  for(i in 1:25){  
    
    df = readWorksheet(wb, sheet = i, header = TRUE)
    if(i == 23 || i == 15)
    {
      l = 3 
    }else{
      lastItem = which(df[,1] == "Total")
      l = lastItem[1]    
    }
    
    if(sheets[i] %in% excepTL){
      df = df[3:l-1,] # ignore the constant species
      l1 = which(grepl(paste("-", sheets[i], sep=""),df[,1]))
      l2 = seq(1,l-2)
      aL = list()
      aL[[1]] = l1 # ignore the constant species
      aL[[2]] = l2[-l1]
      s = c(sheets[i], paste("dh", sheets[i], sep=""))
      for(k in 1:2){
        len = length(aL[[k]])
        lUnInfectedHela = matrix(ncol = len, nrow = 2)
        lHRV1BHela = matrix(ncol = len, nrow = 2)
        lHRV16Hela = matrix(ncol = len, nrow = 2)
        
        lUnInfectedHbec = matrix(ncol = len, nrow = 2)
        lHRV1BHbec = matrix(ncol = len, nrow = 2)
        lHRV16Hbec = matrix(ncol = len, nrow = 2)
        
        lUnInfectedHela[1,1:len] = df[aL[[k]],2]
        lUnInfectedHela[2,1:len] = df[aL[[k]],3]
        
        lHRV1BHela[1,1:len] = df[aL[[k]],4]
        lHRV1BHela[2,1:len] = df[aL[[k]],5]
        
        lHRV16Hela[1,1:len] = df[aL[[k]],6]
        lHRV16Hela[2,1:len] = df[aL[[k]],7]
        
        lUnInfectedHbec[1,1:len] = df[aL[[k]],8]
        lUnInfectedHbec[2,1:len] = df[aL[[k]],9]
        
        lHRV1BHbec[1,1:len] = df[aL[[k]],10]
        lHRV1BHbec[2,1:len] = df[aL[[k]],11]
        
        lHRV16Hbec[1,1:len] = df[aL[[k]],12]
        lHRV16Hbec[2,1:len] = df[aL[[k]],13]
        
        
        class(lUnInfectedHela) = "numeric"
        class(lHRV1BHela) = "numeric"
        class(lHRV16Hela) = "numeric"
        
        class(lUnInfectedHbec) = "numeric"
        class(lHRV1BHbec) = "numeric"
        class(lHRV16Hbec) = "numeric"
        
        colnames(lUnInfectedHela) = c(df[aL[[k]],1])
        colnames(lHRV1BHela) = c(df[aL[[k]],1])
        colnames(lHRV16Hela) = c(df[aL[[k]],1])
        
        colnames(lUnInfectedHbec) = c(df[aL[[k]],1])
        colnames(lHRV1BHbec) = c(df[aL[[k]],1])
        colnames(lHRV16Hbec) = c(df[aL[[k]],1])
        
        unInfectedHela[[s[k]]] = lUnInfectedHbec
        HRV1BHela[[s[k]]] = lHRV1BHbec
        HRV16Hela[[s[k]]] = lHRV16Hbec
        
        unInfectedHbec[[s[k]]] = lUnInfectedHbec
        HRV1BHbec[[s[k]]] = lHRV1BHbec
        HRV16Hbec[[s[k]]] = lHRV16Hbec
        
      }
      
    }else{
      lUnInfectedHela = matrix(ncol = l-2, nrow = 2)
      lHRV1BHela = matrix(ncol = l-2, nrow = 2)
      lHRV16Hela = matrix(ncol = l-2, nrow = 2)
      
      lUnInfectedHbec = matrix(ncol = l-2, nrow = 2)
      lHRV1BHbec = matrix(ncol = l-2, nrow = 2)
      lHRV16Hbec = matrix(ncol = l-2, nrow = 2)
      
      lUnInfectedHela[1,1:(l-2)] = df[3:l-1,2]
      lUnInfectedHela[2,1:(l-2)] = df[3:l-1,3]
      
      lHRV1BHela[1,1:(l-2)] = df[3:l-1,4]
      lHRV1BHela[2,1:(l-2)] = df[3:l-1,5]
      
      lHRV16Hela[1,1:(l-2)] = df[3:l-1,6]
      lHRV16Hela[2,1:(l-2)] = df[3:l-1,7]
      
      lUnInfectedHbec[1,1:(l-2)] = df[3:l-1,8]
      lUnInfectedHbec[2,1:(l-2)] = df[3:l-1,9]
      
      lHRV1BHbec[1,1:(l-2)] = df[3:l-1,10]
      lHRV1BHbec[2,1:(l-2)] = df[3:l-1,11]
      
      lHRV16Hbec[1,1:(l-2)] = df[3:l-1,12]
      lHRV16Hbec[2,1:(l-2)] = df[3:l-1,13]
      
      class(lUnInfectedHela) = "numeric"
      class(lHRV1BHela) = "numeric"
      class(lHRV16Hela) = "numeric"
      
      class(lUnInfectedHbec) = "numeric"
      class(lHRV1BHbec) = "numeric"
      class(lHRV16Hbec) = "numeric"
      
      
      colnames(lUnInfectedHela) = c(df[3:l-1,1])
      colnames(lHRV1BHela) = c(df[3:l-1,1])
      colnames(lHRV16Hela) = c(df[3:l-1,1])
      
      colnames(lUnInfectedHbec) = c(df[3:l-1,1])
      colnames(lHRV1BHbec) = c(df[3:l-1,1])
      colnames(lHRV16Hbec) = c(df[3:l-1,1])
      
      
      unInfectedHela[[sheets[i]]] = lUnInfectedHela
      HRV1BHela[[sheets[i]]] = lHRV1BHela
      HRV16Hela[[sheets[i]]] = lHRV16Hela
      
      unInfectedHbec[[sheets[i]]] = lUnInfectedHbec
      HRV1BHbec[[sheets[i]]] = lHRV1BHbec
      HRV16Hbec[[sheets[i]]] = lHRV16Hbec
    }
  }
  #lipidClass = c(sheets[1:25],"dhCer","dhSM")
  getSigPathways(unInfectedHela, HRV1BHela, 0)
  
}

getSubPathwayPvalue <- function(pathway, unInfData, infData){  
  
  size = length(pathway)-1
  unInfWeight = vector()
  infWeight = vector()
  p = 0
  for(i in 1:size){
    re = pathway[i]
    pro = pathway[i+1]
    unInfWeight[i] = sum(unInfData[[pro]][1,])/sum(unInfData[[re]][1,])
    infWeight[i] = sum(infData[[pro]][1,])/sum(infData[[re]][1,])
  }
  t = t.test(infWeight, unInfWeight, altervative = "greater")
  
  if(is.finite(t$p.value))
  p = t$p.value
  p
}
getSigPathways <- function(unInfData, infData, flag){
  reactant = vector()
  product = vector()
  weight = vector()
  
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    s1 = sum(unInfData[[product[i]]][1,])/sum(unInfData[[reactant[i]]][1,])
    s2 = sum(infData[[product[i]]][1,])/sum(infData[[reactant[i]]][1,])
    weight[i] = s1/s2
  }
  # generate a list of found subpathways with its z_score
  subPathwayL = vector()
  pScoreL = vector()
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
    subPathway = vector()
    subPathway[1] = startNode
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
    p_score = neighborL[1]
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
          p_tmp = getSubPathwayPvalue(subPathway, unInfData, infData)
          if(p_tmp < 0.05){
            # choose this node
            currentNode = nextNode
            visitted[nextNode] = TRUE
            p_score = p_tmp
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
    pScoreL[num] = p_score
    num = num + 1
  }
  dt = data.frame(pathway = factor(c(subPathwayL), ordered = TRUE),
                  score = c(pScoreL))
  dt = dt[order(dt$score, decreasing = TRUE),]
  dt = dt[which(dt$score < 0.05),]
  if(flag == 0){
    write.xlsx(dt, "most_active_pathway_hela_hrv1b.xlsx")  
  } else if(flag == 1){
    write.xlsx(dt, "most_active_pathway_hela_hrv16.xlsx")
  } else if(flag == 2){
    write.xlsx(dt, "most_active_pathway_hbec_hrv1b.xlsx")
  } else {
    write.xlsx(dt, "most_active_pathway_hbec_hrv16.xlsx")
  }
    
}
