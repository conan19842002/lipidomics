library(RMySQL)
library(gplots)
library(RColorBrewer)
library(NMF)
library(PCAmixdata)
library("factoextra")
source("D:/project/lipidomics/data/lipidomicLib.r")
source("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis/findTumourSigPathway.r")
require(XLConnect)
setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis")
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
  wb = XLConnect::loadWorkbook("20130716_Colorectal 2-26 PIPx.xlsx")
  df = XLConnect::readWorksheet(wb, sheet = 3, header = TRUE)
  patientVec = vector()
  patientVec = getPatientVec(10,23,df,patientVec)
  patientVec = getPatientVec(26,33,df,patientVec)
  patientVec = getPatientVec(35,40,df,patientVec)
  patientVec = getPatientVec(43,44,df,patientVec)
  patientVec = getPatientVec(47,60,df,patientVec)
  patientVec = getPatientVec(62,79,df,patientVec)
  patientVec = getPatientVec(82,87,df,patientVec)
  patientVec = getPatientVec(90,95,df,patientVec)
  patientVec = getPatientVec(98,117,df,patientVec)
  patientVec = getPatientVec(120,125,df,patientVec)
  patientVec = getPatientVec(128,135,df,patientVec)
  patientVec = getPatientVec(138,147,df,patientVec)
  patientNum = length(patientVec)/2
  tData = matrix(nrow = 15, ncol = patientNum)
  nData = matrix(nrow = 15, ncol = patientNum)
  colnames(tData) = patientVec[seq(1,length(patientVec),2)]
  colnames(nData) = patientVec[seq(2,length(patientVec),2)]
  
  df = XLConnect::readWorksheet(wb, sheet = 2, header = TRUE)
  rownames(tData) = df[29:43,1]
  rownames(nData) = df[29:43,1]
  colLabel = colnames(df)
  colLabel = colLabel[-1]
  for(colNum in seq(2,22,2)){
    sampleNum = substr(colLabel[colNum],2,nchar(colLabel[colNum]))
    tData[1:15,patientVec[sampleNum]] = as.numeric(df[29:43,colLabel[colNum]])
    nData[1:15,patientVec[sampleNum]] = as.numeric(df[29:43,colLabel[colNum+1]])
  }
  # process second file
  wb = XLConnect::loadWorkbook("20130717_Colorectal 27-52 PIPx.xlsx")
  df = XLConnect::readWorksheet(wb, sheet = 2, header = TRUE)
  colLabel = colnames(df)
  colLabel = colLabel[-1]
  for(colNum in seq(1,21,2)){
    sampleNum = substr(colLabel[colNum],2,nchar(colLabel[colNum]))
    tData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum]])
    nData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum+1]])
  }
  # process third file
  wb = XLConnect::loadWorkbook("20130717_Colorectal 53-77 PIPx.xlsx")
  df = XLConnect::readWorksheet(wb, sheet = 2, header = TRUE)
  colLabel = colnames(df)
  colLabel = colLabel[-1] # get rid of "X" character in the column name
  for(colNum in seq(1,22,2)){
    sampleNum = substr(colLabel[colNum],2,nchar(colLabel[colNum]))
    tData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum]])
    nData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum+1]])
  }
  # process fourth file
  wb = XLConnect::loadWorkbook("20130718_Colorectal 78-103 PIPx.xlsx")
  df = XLConnect::readWorksheet(wb, sheet = 2, header = TRUE)
  colLabel = colnames(df)
  colLabel = colLabel[-1]
  for(colNum in seq(1,21,2)){
    sampleNum = substr(colLabel[colNum],2,nchar(colLabel[colNum]))
    tData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum]])
    nData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum+1]])
  }
  # process fifth file
  wb = XLConnect::loadWorkbook("20130719_Colorectal 104-143 PIPx.xlsx")
  df = XLConnect::readWorksheet(wb, sheet = 2, header = TRUE)
  colLabel = colnames(df)
  colLabel = colLabel[-1]
  for(colNum in seq(1,29,2)){
    sampleNum = substr(colLabel[colNum],2,nchar(colLabel[colNum]))
    tData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum]])
    nData[1:15,patientVec[sampleNum]] = as.numeric(df[2:16,colLabel[colNum+1]])
  }
  
  tData = t(tData)
  nData = t(nData)
  # combine PIPx data with lipid data
  tData1 = getPatientSpecieMatrix(TRUE)
  nData1 = getPatientSpecieMatrix(FALSE)
  colN <- rownames(tData)
  pmtr1 <- tData1[colN,]
  td <- cbind(pmtr1,tData)
  pmtr1 <- nData1[colN,]
  nd <- cbind(pmtr1,nData)
  # this is to filter PI, PIPx
  #nd = nd[,which(grepl("PI",colnames(nd)))]
  #td = td[,which(grepl("PI",colnames(td)))]
  # convert to molar concentration
  td = getMolarData(td,tumour=TRUE)
  nd = getMolarData(nd,tumour=TRUE)
  td = normalisePatientSpecieMatrix(td)
  nd = normalisePatientSpecieMatrix(nd)
  
  td = getNewData(td)
  nd = getNewData(nd)
  
  
}

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
getSigSpeciesMt <- function(df,tData, nData, species, stage){
  if(stage==""){
    s = ""
  }else{
    s = paste(stage,"/",sep="")
  }
  sigSpecVec = vector()
  newSigSpecVec = vector()
  for(i in 1: length(species)){
    controlData = nData[,which(colnames(nData)==species[i])]
    refData = tData[,which(colnames(tData)==species[i])]
    t = t.test(controlData,refData,paired=T)
    specVec = c(specVec,species[i])
    if(is.finite(t$p.value))
    {
      if(t$p.value <.05){
        sigSpecVec = c(sigSpecVec, species[i])
        #massVec = c(massVec, mass[i])
      }
    }
    
  }
  #specVec = as.factor(specVec)
  #massVec = as.factor(massVec)
  #modeVec = rep('positive',length(specVec))
  #modeVec = as.factor(modeVec)
  #penVec = rep(1,length(specVec))
  #penVec = as.factor(penVec)
  #typeVec = rep('metabolite',length(specVec))
  #typeVec = as.factor(typeVec)
  
  sigSpecVec = gsub("DRG", "DG", sigSpecVec)
  sigSpecVec = gsub("TRG", "TG", sigSpecVec)
  sigSpecVec = gsub("MRG", "MG", sigSpecVec)
  
  for(i in 1:length(sigSpecVec)){
    ids = which(df[,4]==as.character(sigSpecVec[i]))
    newSigSpecVec[i] = as.character(df[ids[1],2])
  }
  id1 = which(!is.na(newSigSpecVec))
  d1 = data.frame(newSigSpecVec[id1])
  #write.table(d1, paste("D:/project/lipidomics/data/lipid analysis/hela_hbec/pathway analysis/mass_peak_time_",preTime,"_",currTime,".txt",sep=""), sep="\t",  col.names=FALSE,quote=FALSE,row.names=FALSE)
  write.table(d1, file=paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/sig_terminal_set.txt",sep=""), sep="\t",  col.names=FALSE,quote=FALSE,row.names=FALSE)  
  
}
getNodes <- function(df, na.file){
  tryCatch({
    naDf = read.table(na.file, sep="\t", header = T, fill = T)
    eaDf = read.table(ea.file, sep="\t", header = T, fill = T)
  }, 
  error = function(c){
    c$message <- paste(c$message, " ( in ", na.file, ")")
    stop(c)
    return(NULL)
  })
  if(nrow(naDf)){
    # convert factor column to character
    naDf[, 1] <- sapply(naDf[, 1], as.character)
    mappL = vector()
    #### convert HMDB IDs to real names #########
    for(k in 1:nrow(naDf)){
      re <- df[which(df$HMDB.ID == as.character(naDf[k,1])),]
      if(nrow(re)>0){
        mappL[naDf[k,1]] = as.character(re$abbreviation)
        naDf[k,1] = mappL[naDf[k,1]]
      }
      
    }
    
  }
  naDf$Protein
  
}
getPrizeVas <- function(tData, nData, species){
  pVec = vector()
  for(i in 1: length(species)){
    controlData = nData[,which(colnames(nData)==species[i])]
    refData = tData[,which(colnames(tData)==species[i])]
    t = t.test(controlData,refData,paired=T)
    if(is.finite(t$p.value))
    {
      pVec = c(pVec,-log(t$p.value))
    }else{
      pVec = c(pVec,0)
    }
  }
  
  names(pVec) = species
  pVec
}
computeNetRobustness <- function(df, stage){
  ###################### compute robustsness of solution #####################
  #fract = 0.3 # percentage of the most frequent nodes wil be retained in the optimal forest
  beta = seq(1,20,1)
  mu = seq(0.05,0.4,0.01)
  robustVec = vector()
  bVec = vector()  
  mVec = vector()
  if(stage==""){
    s = ""
  }else{
    s = paste(stage,"/",sep="")
  }
  for(b in 1:length(beta)){
    for(w in 1:length(mu)){
      bb = b - 1
      ww = w - 1
      nodeL = vector()
      edgeL = vector()
      
      noisy_na.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_noisy_nodeattributes_",bb,"_",ww,".tsv",sep="")
      noisy_ea.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_noisy_edgeattributes_",bb,"_",ww, ".tsv",sep="")
      
      na.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_nodeattributes_",bb,"_",ww,".tsv",sep="")
      ea.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_edgeattributes_",bb,"_",ww,".tsv",sep="")
      
      orgNodes =  getNodes(df, na.file)
      noisyNodes = getNodes(df, noisy_na.file)
      if(!is.null(orgNodes) & !is.null(noisyNodes)){
        bVec = c(bVec,beta[b])
        mVec = c(mVec,mu[w])
        overlapNode = intersect(orgNodes, noisyNodes)
        robust = length(overlapNode)/length(orgNodes)
        robustVec = c(robustVec, robust)  
      }
      
      
    }
  }
  dt = data.frame(beta = bVec, mu = mVec, robustness = robustVec)
  dt = dt[order(robustVec, decreasing = T),]
  config.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/config.txt",sep="")
  write.table(dt, file = config.file, sep = "\t", col.names = T,quote=FALSE,row.names = FALSE,
              qmethod = "double")
  
}
getMergedOptimalNet <- function(df, stage){
  if(stage==""){
    s = ""
  }else{
    s = paste(stage,"/",sep="")
  }
  naF = "result_nodeattributes.tsv"
  eaF = "result_edgeattributes.tsv"
  nodeRobustF = "node_robust.txt"
  cytonetF = "cyto_net.sif"
  robustNodeV = vector()
  beta = seq(1,20,1)
  mu = as.character(seq(0.05,0.4,0.01))
  config.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/config.txt",sep="")
  netTb = read.table(config.file, sep="", header = T, fill = T)
  # select parameter sets with highest robustness score for subnetwork defined by a threshold (0.6)
  netTb = netTb[which(netTb$robustness >= 0.6),]
  grVec = list(length = nrow(netTb))
  nodeL = vector()
  for(i in 1:nrow(netTb)){
    bb = which(beta == netTb$beta[i])
    mm = which(mu==as.character(netTb$mu[i]))
    b = bb-1
    m = mm-1
    na.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_nodeattributes_",b,"_",m,".tsv",sep="")
    ea.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/noise/result_edgeattributes_",b,"_",m,".tsv",sep="")
    nodeL = c(nodeL, getNodes(df, na.file))
    # extract subnetwork for each selected parameter sets
    grVec[[i]] = getSubNetwork(df, na.file, ea.file)
  }
  # merge selected subnetworks
  g = graph.union(grVec[[1]], grVec[[2]], byname=TRUE)
  for(i in 3:length(grVec)){
    g = graph.union(g, grVec[[i]], byname=TRUE)
  }
  # compute frequency of each node in the merged network as robustness scores for nodes
  g1 = table(nodeL)
  newNode = sort(g1 , decreasing = T)
  
  reactant = names(head_of(g, E(g)))
  product = names(tail_of(g,E(g)))
  
  
  
  #nodes = union(reactant, product)
  ## write a node prize file
  
  reDf = data.frame(protein = names(newNode), robust = newNode)
  ref = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",nodeRobustF,sep="")
  write.table(reDf, file = ref, sep = "\t", col.names = TRUE,quote=FALSE,row.names = FALSE, 
              qmethod = "double")
  sif = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",cytonetF,sep="")
  write.table(data.frame(reactant, "pp", product), file = sif, sep = "\t", col.names = FALSE,quote=FALSE,row.names = FALSE,
              qmethod = "double")
}
getSubNetwork <- function(df, naDf, eaDf){
  naDf = read.table(na.file, sep="\t", header = T, fill = T)
  eaDf = read.table(ea.file, sep="\t", header = T, fill = T)
  # convert factor column to character
  naDf[, 1] <- sapply(naDf[, 1], as.character)
  mappL = vector()
  #### convert HMDB IDs to real names #########
  for(j in 1:nrow(naDf)){
    re <- df[which(df$HMDB.ID == as.character(naDf[j,1])),]
    if(nrow(re)>0){
      mappL[naDf[j,1]] = as.character(re$abbreviation)
      naDf[j,1] = mappL[naDf[j,1]]
    }
    
  }
  ###### extract interactions from edge attributes ########
  eaDf[, 1] <- sapply(eaDf[, 1], as.character)
  
  reactant = vector(length = nrow(eaDf))
  product = vector(length = nrow(eaDf))
  weight = vector()
  for(j in 1:nrow(eaDf)){
    tokens = unlist(strsplit(eaDf[j,1], split = " "))
    weight = c(weight, eaDf[j,2])
    if(grepl("HMDB", tokens[1])){
      reactant[j] = mappL[tokens[1]]
    }else{
      reactant[j] = tokens[1]
    }
    if(grepl("HMDB", tokens[3])){
      product[j] = mappL[tokens[3]]
    }else{
      product[j] = tokens[3]  
    }
  }
  #nodes = union(reactant, product)
  #sizeV = vector()
  #for(i in 1:length(nodes)){
  #  sizeV[i] = mean(naDf[which(naDf$Protein == nodes[i]),]$Prize)
  #}
  edgeL = data.frame(from = reactant, to = product, weight = weight)
  #ig <- graphBAM(edgeL, edgemode="undirected")
  ig = graph.data.frame(edgeL,directed = F)
  #V(ig)$size<- sizeV
  ig <- simplify(ig, remove.multiple = T, remove.loops = T)
  ig
}
getOptimalNet <- function(stage){
  naF = "result_nodeattributes.tsv"
  eaF = "result_edgeattributes.tsv"
  if(stage==""){
    s = ""
  }else{
    s = paste(stage,"/",sep="")
  }
  #naF = "shuffle/result_randomTerminals_nodeattributes.tsv"
  #eaF = "shuffle/result_randomTerminals_edgeattributes.tsv"
  nodePriF = "node_prize.txt"
  cytonetF = "cyto_net.sif"
  
  na.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",naF,sep="")
  ea.file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",eaF,sep="")
  
  naDf = read.table(na.file, sep="\t", header = T, fill = T)
  eaDf = read.table(ea.file, sep="\t", header = T, fill = T)
  # convert factor column to character
  naDf[, 1] <- sapply(naDf[, 1], as.character)
  mappL = vector()
  #### convert HMDB IDs to real names #########
  for(j in 1:nrow(naDf)){
    re <- df[which(df$HMDB.ID == as.character(naDf[j,1])),]
    if(nrow(re)>0){
      mappL[naDf[j,1]] = as.character(re$abbreviation)
      naDf[j,1] = mappL[naDf[j,1]]
    }
    
  }
  ###### extract interactions from edge attributes ########
  eaDf[, 1] <- sapply(eaDf[, 1], as.character)
  
  reactant = vector(length = nrow(eaDf))
  product = vector(length = nrow(eaDf))
  weight = vector()
  for(j in 1:nrow(eaDf)){
    tokens = unlist(strsplit(eaDf[j,1], split = " "))
    weight = c(weight, eaDf[j,2])
    if(grepl("HMDB", tokens[1])){
      reactant[j] = mappL[tokens[1]]
    }else{
      reactant[j] = tokens[1]
    }
    if(grepl("HMDB", tokens[3])){
      product[j] = mappL[tokens[3]]
    }else{
      product[j] = tokens[3]  
    }
  }
  edgeL = data.frame(from = reactant, to = product, weight = weight)
  nodes = union(reactant, product)
  #ig = graph.data.frame(edgeL,directed = F)
  #ig <- simplify(ig, remove.multiple = T, remove.loops = T)
  #gr = ftM2graphNEL(as.matrix(edgeL[,1:2]), W= edgeL[,3], V=unique(naDf$Protein), edgemode="undirected")
  #ig <- graph_from_graphnel(gr)
  #nodes = V(ig)$name
  sizeV = vector()
  for(i in 1:length(nodes)){
    if(grepl("-", nodes[i])){
      sizeV[i] = naDf[which(naDf$Protein == nodes[i]),]$Prize
    }else{
      sizeV[i] = 0
    }
  }
  #V(ig)$label.cex <- 1         ## text size
  #V(ig)$size      <- sizeV
  #idx = which(grepl("-",V(ig)$name))
  #V(ig)[idx]$shape <- "circle"
  #V(ig)[-idx]$shape <- "square"
  #ig$layout = layout_with_fr
  #pdf(file = paste("D:/project/lipidomics/data/lipid analysis/hela_hbec/steiner forest/results/","time_", preTime,"_",currTime,"/optimal_wei_net.pdf", sep=""))
  #ig = simplify(ig, edge.attr.comb=list(weight="sum"))
  #plot(ig)
  #dev.off()
  # save infor to file
  reDf = data.frame(protein = nodes, prize = sizeV)
  ref = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",nodePriF,sep="")
  sif = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/",cytonetF,sep="")
  write.table(reDf, file = ref, sep = "\t", col.names = F,quote=FALSE,row.names = F, 
              qmethod = "double")
  write.table(data.frame(reactant, "pp", product), file = sif, sep = "\t", col.names = F,quote=FALSE,row.names = F,
              qmethod = "double")
  
}
prepareDataForSteiner <- function(df, tData, nData, stage){
  if(stage==""){
    s = ""
  }else{
    s = paste(stage,"/",sep="")
  }
  infor = read.table("D:/project/lipidomics/data/lipid analysis/hela_hbec/pathway analysis/iRef13_hmdb_recon_net_woSpace.txt", sep="", header = T, fill = T)
  colnames(infor) = c("r","p","s")
  speciesL = colnames(tData)
  
  reactant = infor[,1]
  product = infor[,2]
  metabo= union(reactant, product)
  
  # look for species in HMDB database
  mappedSpecies =  as.vector(df[which(df[,4] %in% speciesL),1])
  exprStr = paste(mappedSpecies,collapse="|")
  tmp1<-as.character(infor$r)
  tmp2<-as.character(infor$p)
  mappedSpecies = mappedSpecies[-which(mappedSpecies == "Sphingosine" | mappedSpecies == "Sphingosine 1-phosphate")]
  exptL = c("sphingosine(1+)","sphingosine_1-phosphate(1-)")
  mappL = c("Sphingosine","Sphingosine 1-phosphate")
  for(i in 1:length(mappedSpecies)){
    ids = which(grepl(as.character(mappedSpecies[i]),infor$r,fixed=T))
    tmp1[ids]<- as.character(df[which(df[,1]==as.character(mappedSpecies[i])),2])
    infor$r<-factor(tmp1)
    
    ids = which(grepl(as.character(mappedSpecies[i]),infor$p,fixed=T))
    tmp2[ids]<- as.character(df[which(df[,1]==as.character(mappedSpecies[i])),2])
    infor$p<-factor(tmp2)
  }
  for(i in 1:2){
    ids = which(grepl(as.character(exptL[i]),infor$r,fixed=T))
    tmp1[ids]<- as.character(df[which(df[,1]==as.character(mappL[i])),2])
    infor$r<-factor(tmp1)
    
    ids = which(grepl(as.character(exptL[i]),infor$p,fixed=T))
    tmp2[ids]<- as.character(df[which(df[,1]==as.character(mappL[i])),2])
    infor$p<-factor(tmp2)
  }
  if(stage == ""){
    write.table(infor, file = "D:/project/lipidomics/data/lipid analysis/iRef13_hmdb_recon_net_woSpace_a.txt", sep = "\t", col.names = T,quote=FALSE,row.names = FALSE,
                qmethod = "double")
    
    infor = read.table("D:/project/lipidomics/data/lipid analysis/iRef13_hmdb_recon_net_woSpace_a.txt", sep="", header = T, fill = T)
    #colnames(infor) = c("r","p","s")
    dt1 = infor[which(grepl("HMDB",infor$r)),]
    dt2 = infor[which(grepl("HMDB",infor$p)),]
    dt = rbind(dt1,dt2)
    hmdbIdL = union(dt1$r,dt2$p)
    
    ###### interactome file #########
    write.table(dt, file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/edgeFile.txt",sep=""), sep = "\t", col.names = FALSE,quote=FALSE,row.names = FALSE, 
                qmethod = "double")
  }else{
    write.table(infor, file = "D:/project/lipidomics/data/lipid analysis/iRef13_hmdb_recon_net_woSpace_colorectal_stage.txt", sep = "\t", col.names = T,quote=FALSE,row.names = FALSE,
                qmethod = "double")
    
    infor = read.table("D:/project/lipidomics/data/lipid analysis/iRef13_hmdb_recon_net_woSpace_colorectal_stage.txt", sep="", header = T, fill = T)
    colnames(infor) = c("r","p","s")
    dt1 = infor[which(grepl("HMDB",infor$r)),]
    dt2 = infor[which(grepl("HMDB",infor$p)),]
    dt = rbind(dt1,dt2)
    hmdbIdL = union(dt1$r,dt2$p)
    
    ###### interactome file #########
    write.table(dt, file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/edgeFile.txt",sep=""), sep = "\t", col.names = FALSE,quote=FALSE,row.names = FALSE, 
                qmethod = "double")
  }
  
  pVals = getPrizeVas(tData,nData, speciesL)
  WéWWW
  #### dt for significant species only ####
  #dt = rbind(dt1[which(dt1$r %in% as.character(prize[,1])),],dt2[which(dt2$p %in% as.character(prize[,1])),])
  #dtt1 = df[which(df$HMDB.ID %in% as.character(prize[,1])),]
  dtt1 = df[which(df$HMDB.ID %in% hmdbIdL),]
  dtt2 = data.frame(metabolite = dtt1$HMDB.ID, mapped = dtt1$abbreviation)
  for(i in 1:nrow(dtt2)){
    score = pVals[which(names(pVals)== as.character(dtt2$mapped[i]))]
    dtt2$score[i] = score
  }
  dtt = data.frame(dtt2$metabolite,dtt2$score)
  
  write.table(dtt, file = paste("D:/project/lipidomics/data/lipid analysis/steiner/",s,"results/terminal_set.txt",sep=""), sep = "\t", col.names = FALSE,quote=FALSE,row.names = FALSE,
              qmethod = "double")
  #prepare significant species list
  #getSigSpeciesMt(df,tData, nData, speciesL, stage)
  
}
runSteinerForest <- function(tData, nData){
  df = read.csv("D:/project/lipidomics/data/lipid analysis/hmdb_metabolites/hmdb_id.csv", header = TRUE)
  # prepare data for whole data
  prepareDataForSteiner(df, tData,nData,"")
  computeNetRobustness(df, "")  
  getOptimalNet("")
  # prepare data for stages
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
  for(s in 1:length(stages)){
    tPS = stagePatients[[stages[s]]]
    nPS = nData[patientList[[s]],]
    prepareDataForSteiner(df,tPS,nPS,stages[s])
    
  }
  for(s in 1:length(stages)){
    computeNetRobustness(df,stages[s])
  }
  
}
# this function is to find whether there are differences in pathway analyses between two tumour classes
# class 1 : tumours which have more C36-PIP3
# class 2 : tumours which have more C38-PIP3
doTumourSigPathway <- function(){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  # tumours which have more c36-PIP3
  c36Id = vector()
  # tumours which have more c38-PIP3
  c38Id = vector()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    td = filterClassPSM_1(stagePatients[[stages[s]]],"PIP3")
    nd = filterClassPSM_1(nData[patientList[[s]],],"PIP3")
    carbonAcylList_td = vector()
    carbonAcylList_nd = vector()
    listOfLipid = colnames(td)
    listOfCarbonLength = vector()
    
    if(length(listOfLipid) > 0){
      
      for(j in 1:length(listOfLipid)){
        acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
        listOfCarbonLength[j] = acyl[1]  
      }
      listOfCarbonLength = unique(listOfCarbonLength)
      mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
      rownames(mtAcylLen) = rownames(td)
      colnames(mtAcylLen) = listOfCarbonLength
      
      for(j in 1:length(listOfCarbonLength)){
        
        r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
        # sum all the lipid of the same acyl carbon length
        if(length(r) > 1){
          lc1 = rowSums(td[,r])
          lc2 = rowSums(nd[,r])  
        }else{
          lc1 = td[,r]
          lc2 = nd[,r]
        }
        
        mtAcylLen[,j] = lc1/lc2
      }
      
      
      ind = which(mtAcylLen == 0)
      mtAcylLen = log10(mtAcylLen)
      mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
      mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
      mtAcylLen[which(is.nan(mtAcylLen))] = 0
      if(all(mtAcylLen==0 | is.infinite(mtAcylLen)))
        next
      # find max values in each column
      #maxMt = apply(mtAcylLen,2,max)
      # subtract each column with the max values
      #mtAcylLen = sweep(-mtAcylLen,2,maxMt,"+")
      # divide by the max values
      #mtAcylLen = sweep(mtAcylLen,2,maxMt,"/")
      # Take all tumour of 70% upper quartile
      strId = rownames(mtAcylLen)
      c36Id = append(c36Id,strId[which(mtAcylLen[,"36"] > 0 & mtAcylLen[,"36"] > mtAcylLen[,"38"])])
      c38Id = append(c38Id,strId[which(mtAcylLen[,"38"] > 0 & mtAcylLen[,"38"] > mtAcylLen[,"36"])])
    }
    
  }
  
  #c36Id = c( "102",  "1245", "1192", "93" ,  "85"  , "1234" ,"1382", "1054", "73" ,  "84"  , "86"  , "90" , "97" ,  "98" ,  "1375", "68" ,  "69",   "94" ,  "1058", "37",   "1088")
  
  #c38Id = c("1361", "81",   "91",   "77",   "1286",   "1357",   "82")
  c36TData = tData[c36Id,]
  c38TData = tData[c38Id,]
  c36NData = nData[c36Id,]
  c38NData = nData[c38Id,]
  #setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis")
  saveToFileSigPathway(c36TData, c36NData,"C36_PIP3","greater")
  saveToFileSigPathway(c36TData, c36NData,"C36_PIP3","less")
  saveToFileSigPathway(c38TData, c38NData,"C38_PIP3","greater")
  saveToFileSigPathway(c38TData, c38NData,"C38_PIP3","less")
}

getEdgeWeight <- function(re, pro, tData, nData,alt){
  trData = filterClassPSM_1(tData,re)
  tpData = filterClassPSM_1(tData,pro)
  nrData = filterClassPSM_1(nData,re)
  npData = filterClassPSM_1(nData,pro)
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
  tDataTest = tSumV[[1]][ind]/tSumV[[2]][ind]
  ind = which(nSumV[[2]] > 0)                       
  nDataTest = nSumV[[1]][ind]/nSumV[[2]][ind]
  z_score = 0
  out <- tryCatch({
    if(length(tDataTest)==length(nDataTest)){
      t = t.test(tDataTest,nDataTest, alternative= alt, paired=T)  
    }else{
      t = t.test(tDataTest,nDataTest, alternative= alt)  
    }
    
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

findSigPathways <- function(tData,nData,alt){
  reactions = as.matrix(read.csv("reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  weight = vector()
  for(i in 1:nrow(reactions)){
    reactant[i] = reactions[i,1]
    product[i] = reactions[i,2]  
    weight[i] = getEdgeWeight(reactant[i], product[i], tData, nData,alt)
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
  if(alt=="greater"){
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/most_active_pathway.xlsx")  
  }else{
    write.xlsx(dt, "D:/project/lipidomics/data/pathway analysis/most_inactive_pathway.xlsx")  
  }
  
}

findSigPathwayByStage <- function(tData, nData){
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  weight = vector()
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
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
      if(z_score < 1.645)
        next
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
    write.xlsx(dt, paste("D:/project/lipidomics/data/pathway analysis/most_active_pathway_stage_",stages[s],".xlsx",sep=""))  
  }
}
plotSigPathwayLayerByStage <- function(tData, nData){
  reactions = as.matrix(read.csv("D:/project/lipidomics/reaction1s.csv", header = T))
  # get set of reactants and products
  reactant = vector()
  product = vector()
  # all the nodes of the pathway are sub-species of lipid classes
  # therefore notice TRG, DRG not TG, DG
  pathway = c("PIP3", "PIP2","PIP")
  lipidClass = c("PIP3", "PIP2","PIP")
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
  for(s in 1:length(stages)){
    # get set of subpathways of the different layers
    subpathways = list()
    # get patient-species matrix from normal and tumour
    
    iSpeciesList = list()
    # get set of lipid species which involve in the pathway
    
    for(i in 1:length(lipidClass)){
      # get all sub-species from the tumour
      tSpecies = colnames(filterClassPSM_1(stagePatients[[stages[s]]],lipidClass[i]))
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
        tPS = filterClassPSM_1(stagePatients[[stages[s]]],lipidClass[j])
        nPS = filterClassPSM_1(nData[patientList[[s]],],lipidClass[j])
        tSubPathways[[species]] = tPS[,which(colnames(tPS) == species)]
        nSubPathways[[species]] = nPS[,which(colnames(nPS) == species)]
        # now check for non-zero values
        id1 = which(tSubPathways[[species]] == 0)
        id2 = which(nSubPathways[[species]] == 0)
        id = c(id1,id2)
        if(length(id) > 0){
          tSubPathways[[species]] = tSubPathways[[species]][-id]
          nSubPathways[[species]] = nSubPathways[[species]][-id]
        }
      }
    }
    tSubPathWayNode = list()
    nSubPathWayNode = list()
    size = length(pathway)-1
    for(i in 1:length(iSharedList)){
      zc = vector()
      k = 1
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
    pdf(file = paste("active_pathway_by_layer_",stages[s],".pdf",sep=""))
    plot(z_score, pch=20,  main = title, col = zColor, xaxt="n", xlab="Layer", ylab="Z-Score")
    axis(1, at=1:length(iSharedList), labels=iSharedList)
    color = unique(zColor)
    if(length(color) == 1){
      if(color == "black"){
        legend("bottomleft", legend = c("non active"), fill =c("black"))  
      }else{
        legend("bottomleft", legend = c("active"), fill =c("red"))  
      }
    }else
      legend("topleft", legend = c("non active", "active"), fill = c("black","red"))
    dev.off()
  }
}
plotSigPathwayByLayer <- function(tData, nData){
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/pathway analysis")
  # all the nodes of the pathway are sub-species of lipid classes
  # therefore notice TRG, DRG not TG, DG
  pathway = c("PIP","PI")
  lipidClass = c("PIP","PI")
  # get set of subpathways of the different layers
  subpathways = list()
  # get patient-species matrix from normal and tumour
  
  iSpeciesList = list()
  # get set of lipid species which involve in the pathway
  
  for(i in 1:length(lipidClass)){
    # get all sub-species from the tumour
    tSpecies = colnames(filterClassPSM_1(tData,lipidClass[i]))
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
      tPS = filterClassPSM_1(tData,lipidClass[j])#tData[,which(tSpeciesVec[seq(2,l,2)] == lipidClass[j])]
      nPS = filterClassPSM_1(nData,lipidClass[j])#nData[,which(tSpeciesVec[seq(2,l,2)] == lipidClass[j])]
      tSubPathways[[species]] = tPS[,which(colnames(tPS) == species)]
      nSubPathways[[species]] = nPS[,which(colnames(nPS) == species)]
      # now check for non-zero values
      id1 = which(tSubPathways[[species]] == 0)
      id2 = which(nSubPathways[[species]] == 0)
      id = c(id1,id2)
      if(length(id) > 0){
        tSubPathways[[species]] = tSubPathways[[species]][-id]
        nSubPathways[[species]] = nSubPathways[[species]][-id]
      }
    }
  }
  tSubPathWayNode = list()
  nSubPathWayNode = list()
  size = length(pathway)-1
  for(i in 1:length(iSharedList)){
    zc = vector()
    k = 1
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
    if(z_score[iSharedList[i]] > 1.645){
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
  pdf(file = "active_pathway_by_layer.pdf")
  plot(z_score, pch=20,  main = title, col = zColor, xaxt="n", xlab="Layer", ylab="Z-Score")
  axis(1, at=1:length(iSharedList), labels=iSharedList)
  color = unique(zColor)
  if(length(color) == 1){
    if(color == "black"){
      legend("bottomleft", legend = c("non active"), fill =c("black"))  
    }else{
      legend("bottomleft", legend = c("active"), fill =c("red"))  
    }
  }else
    legend("topleft", legend = c("non active", "active"), fill = c("black","red"))
  dev.off()
}
plotPipxRatioByAcylBond <- function(tData, nData, speciesClass){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pipxL = c("PIP3","PIP2","PIP")
  
  for(k in 1:3){
    tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
    npipxVec = unlist(strsplit(colnames(nData), split ="-"))
    td = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == pipxL[k])]
    nd = nData[,which(npipxVec[seq(2,length(npipxVec),2)] == pipxL[k])]
    colN = colnames(td)
    carbonAcylList_N = list()
    doubleBondList_N = list()
    carbonAcylList_T = list()
    doubleBondList_T = list()
    listOfCarbonLength = vector()
    carbonLengthToRemove = vector()
    doubleBondToRemove = vector()
    tTestCarbonLengthToRemove = vector()
    tTestDoubleBondToRemove = vector()
    listOfDoubleBond = vector()
    for(i in 1:length(colN)){
      acyl = unlist(strsplit(unlist(strsplit(colN[i], split ="-"))[1], split = ":"))
      listOfCarbonLength[i] = acyl[1]
      listOfDoubleBond[i] = acyl[2]
    }
    listOfCarbonLength = unique(listOfCarbonLength)
    listOfDoubleBond = unique(listOfDoubleBond)
    listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
    listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
    
    for(i in 1:length(listOfCarbonLength)){
      carbonLengthToRemove = vector()
      
      r = which(grepl(paste(listOfCarbonLength[i],":", sep=""),colN))
      # sum all the lipid of the same acyl carbon length
      if(length(r) > 1){
        lc1 = rowSums(td[,r])  
        lc2 = rowSums(nd[,r])
      }
      else {
        lc1 = td[,r]  
        lc2 = nd[,r]
      }
      ind1 = which(lc1 == 0)
      ind2 = which(lc2 == 0)
      carbonLengthToRemove = c(id1,id2)
      
      if(length(carbonLengthToRemove) > 0){
        carbonAcylList_T[[listOfCarbonLength[i]]] = lc1[-carbonLengthToRemove]
        carbonAcylList_N[[listOfCarbonLength[i]]] = lc2[-carbonLengthToRemove]
      }
      
      if(length(carbonAcylList_T[[listOfCarbonLength[i]]]) < 2 | length(carbonAcylList_N[[listOfCarbonLength[i]]]) < 2){
        tTestCarbonLengthToRemove[length(tTestCarbonLengthToRemove)+1] = i
      }
    }
    # check whether there are enough data for acyl chain length test
    if(length(tTestCarbonLengthToRemove) > 0)
      listOfCarbonLength = listOfCarbonLength[-tTestCarbonLengthToRemove]
    # do the same for double bonds
    for(i in 1:length(listOfDoubleBond)){
      doubleBondToRemove = vector()
      r = which(grepl(paste(":",listOfDoubleBond[i], sep = ""),colN))
      # sum all the lipid of the same double bond number
      if(length(r) > 1){
        lb1 = rowSums(td[,r])  
        lb2 = rowSums(nd[,r])
      }
      else{
        lb1 = td[,r]
        lb2 = nd[,r]
      }
      
      ind1 = which(lb1 == 0)
      ind2 = which(lb2 == 0)
      doubleBondToRemove = c(id1,id2)
      
      if(length(doubleBondToRemove) > 0){
        doubleBondList_T[[listOfDoubleBond[i]]] = lb1[-doubleBondToRemove]
        doubleBondList_N[[listOfDoubleBond[i]]] = lb2[-doubleBondToRemove]
      }
      if(length(doubleBondList_T[[listOfDoubleBond[i]]]) < 2 | length(doubleBondList_N[[listOfDoubleBond[i]]]) < 2){
        tTestDoubleBondToRemove[length(tTestDoubleBondToRemove)+1] = i
      }
    }
    # check whether there are enough data for double bonds test
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
    # box plots
    pdf(file = paste("Acyl_Bonds_",speciesClass[k], ".pdf", sep = ""))
    #lmts <- range(carbonAcylList,doubleBondList)
    par(mfrow=c(2,1))
    boxplot(carbonAcylList,main = speciesClass[k], las = 1, names = listOfCarbonLength , ylab = paste(speciesClass[k]," Ratio (T:N)"), xlab = "Acyl carbon chain length", col = tCarbonAcylTestColVec)
    abline(h=0, col = "blue")
    boxplot(doubleBondList, las = 1, names = listOfDoubleBond , ylab = paste(speciesClass[k]," Ratio (T:N)"), xlab = "Number of double bonds" , col = tDoubleBondTestColVec)
    abline(h=0, col = "blue")
    dev.off()
  }
  
}
plotPipxRatio <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pipxL = c("PIP3","PIP2","PIP","PI")
  for(i in 1:length(pipxL)){
    tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
    npipxVec = unlist(strsplit(colnames(nData), split ="-"))
    td = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == pipxL[i])]
    nd = nData[,which(npipxVec[seq(2,length(npipxVec),2)] == pipxL[i])]
    pipxRatio = list()
    pipxRatioPvalue = vector()
    testColVec = vector()
    idToRemove = vector()
    td1 = vector()
    nd1 = vector()
    colN = colnames(td)
    for(j in 1:ncol(td)){
      td1 = td[,j]
      nd1 = nd[,j]
      id1 = which(nd1 == 0)
      id2 = which(td1 == 0)
      idToRemove = c(id1,id2)
      if(length(idToRemove) > 0){
        td1 = td1[-idToRemove]
        nd1 = nd1[-idToRemove]  
      }
      pipxRatio[[j]] = log10(td1/nd1)
      t = t.test(td1,nd1,paired=T)
      pipxRatioPvalue[j] = t$p.value  
    }
    
    for(j in 1:length(pipxRatioPvalue)){
      if(!is.nan(pipxRatioPvalue[j]) && pipxRatioPvalue[j] < 0.05){
        testColVec[j] = "red"
      }else{
        testColVec[j] = "white"  
      }
      
    }
    pdf(file = paste(pipxL[i],"_ratio.pdf", sep = ""))
    boxplot(pipxRatio,main = paste(pipxL[i], " ratio (T:N)", sep= ""), las = 2, names = colN , ylab = paste(pipxL[i], " ratio (T:N)",sep=""), col = testColVec)
    abline(h=0, col = "blue")
    dev.off()
  }
}
plotPipxRatioBox <- function(tData, nData){
  tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
  npipxVec = unlist(strsplit(colnames(nData), split ="-"))
  td1 = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == "PIP3")]
  td2 = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == "PIP2")]
  td =  td1/td2
  nd1 = nData[,which(tpipxVec[seq(2,length(npipxVec),2)] == "PIP3")]
  nd2 = nData[,which(tpipxVec[seq(2,length(npipxVec),2)] == "PIP2")]
  nd =  nd1/nd2
  colN = colnames(td1)
  acylVec = vector()
  listOfPipxRatio = list()
  pipxRatioTestPValue = vector()
  for(i in 1:length(colN)){
    acylVec[i] = unlist(strsplit(colN[i], split ="-"))[1]
  }
  for(i in 1:length(acylVec)){
    testTd = td[,i]
    testNd = nd[,i]
    id1 = which(testTd == 0 | is.infinite(testTd))
    id2 = which(testNd == 0 | is.infinite(testNd))
    id = c(id1,id2)
    if(length(id) > 0){
      testTd = testTd[-id]
      testNd = testNd[-id]
    }
    t = t.test(testTd,testNd, paired = T)
    listOfPipxRatio[[i]] = log10(testTd/testNd)
    pipxRatioTestPValue[i] = t$p.value
  }
  testColVec = vector(length = length(pipxRatioTestPValue))
  for(i in 1:length(pipxRatioTestPValue)){
    if(!is.nan(pipxRatioTestPValue[i]) && pipxRatioTestPValue[i] < 0.05){
      testColVec[i] = "red"
    }else{
      testColVec[i] = "white"  
    }
  }
  pdf(file = "pip3_pip2_ratio_box_plot.pdf")
  boxplot(listOfPipxRatio,main = "PIP3:PIP2 ratio (T:N)", las = 1, names = acylVec , ylab = "PIP3:PIP2 ratio (T:N)" , xlab = "Acyl chain", col = testColVec)
  abline(h=0, col = "blue")
  dev.off()
  
}
plotPipxSpeciesHeatMap <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pr = matrix(nrow = nrow(tData), ncol = ncol(tData))
  pr = tData/nData
  ind = which(pr == 0)
  prMat = log10(pr)
  prMat[ind] = min(prMat[which(is.finite(prMat))])
  prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
  pdf(file = "pipx_species_heat_map.pdf", width=15,height=8, onefile = T)
  Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
            rep("brown",323))
  heatmap.2(t(prMat), main = paste("PIPx species Heat Map (T:N)"),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
            scale="column", key=T, keysize=0.5, density.info="none",
            trace="none",cexCol=1.2, RowSideColors=Label,
            lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
            lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
  dev.off()
}
plotPipxSpeciesHeatMapByStage <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
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
    pr = matrix(nrow = length(patientList[[s]]), ncol = ncol(tData))
    rownames(pr) = rownames(stagePatients[[stages[s]]])
    td = stagePatients[[stages[s]]]
    nd = nData[patientList[[s]],]
    pr = td/nd
    ind = which(pr == 0)
    prMat = log10(pr)
    prMat[ind] = min(prMat[which(is.finite(prMat))])
    prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
    prMat[which(is.nan(prMat))] = 0
    pdf(file = paste("pipx_species_heat_map_stage_",stages[s],".pdf",sep=""), width=15,height=8, onefile = T)
    Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
              rep("brown",323))
    heatmap.2(t(prMat), main = paste("PIPx Species Heat Map (T:N) (", stages[s],")", sep=""),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
              scale="column", key=T, keysize=0.5, density.info="none",
              trace="none",cexCol=1.2, RowSideColors=Label,
              lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
              lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
    dev.off()
  }
  
}
plotPipxHeatMapByStage <- function(tData,nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  lipidArray = c("PIP3","PIP2","PIP","PI")
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
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
    pr = matrix(nrow = length(patientList[[s]]), ncol = length(lipidArray))
    colnames(pr) = lipidArray
    rownames(pr) = rownames(stagePatients[[stages[s]]])
    tSumV = vector()
    nSumV = vector()
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      if(is.vector(td)){
        tSumV = td
      }else{
        tSumV = rowSums(td)  
      }
      if(is.vector(nd)){
        nSumV = nd
      }else{
        nSumV = rowSums(nd)  
      }
      pr[,i] = tSumV/nSumV
    }
    ind = which(pr == 0)
    prMat = log10(pr)
    prMat[ind] = min(prMat[which(is.finite(prMat))])
    prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
    prMat[which(is.nan(prMat))] = 0
    pdf(file = paste("pipx_heat_map_stage_",stages[s],".pdf",sep=""), width=15,height=8, onefile = T)
    Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
              rep("brown",323))
    heatmap.2(t(prMat), main = paste("Pipx Heat Map (T:N) (", stages[s],")", sep=""),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
              scale="column", key=T, keysize=0.5, density.info="none",
              trace="none",cexCol=1.2, RowSideColors=Label,
              lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
              lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
    dev.off()
  }
  
}
plotPipxHeatMap <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pr = matrix(nrow = nrow(tData), ncol = 4)
  pipxL = c("PIP3","PIP2","PIP", "PI")
  colnames(pr) = pipxL
  rownames(pr) = rownames(tData)
  tSumV = vector()
  nSumV = vector()
  
  for(i in 1:4){
    tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
    npipxVec = unlist(strsplit(colnames(nData), split ="-"))
    td = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == pipxL[i])]
    nd = nData[,which(npipxVec[seq(2,length(npipxVec),2)] == pipxL[i])]
    tSumV = rowSums(td)
    nSumV = rowSums(nd)  
    pr[,i] = tSumV/nSumV
  }
  
  ind = which(pr == 0)
  prMat = log10(pr)
  prMat[ind] = min(prMat[which(is.finite(prMat))])
  prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
  #prMat[which(is.nan(prMat))] = min(prMat[which(is.finite(prMat))])-4
  pdf(file = "D:/project/lipidomics/data/lipid analysis/pipx analysis/pipx_heat_map.pdf", width=15,height=8, onefile = T)
  Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
            rep("brown",323))
  heatmap.2(t(prMat), main = paste("PIPx Heat Map (T:N)"),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
            scale="column", key=T, keysize=0.5, density.info="none",
            trace="none",cexCol=1.2, RowSideColors=Label,
            lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
            lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
  dev.off()
}
plotPipxAcylSpeciesPatientHeatMap <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pr = matrix(nrow = nrow(tData), ncol = 4)
  pipxL = c("PI","PIP","PIP2", "PIP3")
  colnames(pr) = pipxL
  rownames(pr) = rownames(tData)
  tSumV = vector()
  nSumV = vector()
  mtAcylLengthList = list()
  mtDoubleBondList = list()
  for(i in 1:4){
    tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
    npipxVec = unlist(strsplit(colnames(nData), split ="-"))
    td = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == pipxL[i])]
    nd = nData[,which(npipxVec[seq(2,length(npipxVec),2)] == pipxL[i])]
    carbonAcylList_td = vector()
    doubleBondList_td = vector()
    carbonAcylList_nd = vector()
    doubleBondList_nd = vector()
    listOfLipid = colnames(td)
    listOfCarbonLength = vector()
    carbonLengthToRemove = vector()
    doubleBondToRemove = vector()
    
    if(length(listOfLipid) > 0){
      listOfDoubleBond = vector()
      listOfCarbonLength = c("36","38")
      listOfDoubleBond = c("1","2","3","4")
      listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
      listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
      mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
      mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
      rownames(mtAcylLen) = rownames(td)
      colnames(mtAcylLen) = listOfCarbonLength
      rownames(mtDoubleBond) = rownames(td)
      colnames(mtDoubleBond) = listOfDoubleBond
      for(j in 1:length(listOfCarbonLength)){
        
        r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
        # sum all the lipid of the same acyl carbon length
        if(length(r) > 1){
          lc1 = rowSums(td[,r])
          lc2 = rowSums(nd[,r])  
        }else{
          lc1 = td[,r]
          lc2 = nd[,r]  
        }
        
        mtAcylLen[,j] = lc1/lc2
      }
      for(j in 1:length(listOfDoubleBond)){
        r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
        # sum all the lipid of the same double bond number
        if(length(r) > 1){
          lb1 = rowSums(td[,r])
          lb2 = rowSums(nd[,r]) 
        }else{
          lb1 = td[,r]
          lb2 = nd[,r]  
        }
        mtDoubleBond[,j] = lb1/lb2
      }
      ind = which(mtAcylLen == 0)
      mtAcylLen = log10(mtAcylLen)
      mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
      mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
      
      ind = which(mtDoubleBond == 0)
      mtDoubleBond = log10(mtDoubleBond)
      mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
      mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
      mtAcylLengthList[[pipxL[i]]] = mtAcylLen
      mtDoubleBondList[[pipxL[i]]] = mtDoubleBond
    }
  }
  patientNum = rownames(td)
  for(i in 1:nrow(td)){
    prAcylLength = matrix(nrow = 4, ncol = 2)
    prDoubleBond = matrix(nrow = 4, ncol = 4)
    prAcylLength = rbind(rbind(rbind(mtAcylLengthList[["PI"]][i,],mtAcylLengthList[["PIP"]][i,]),mtAcylLengthList[["PIP2"]][i,]),mtAcylLengthList[["PIP3"]][i,])
    
    pdf(file = paste("acyl_species_patient_", patientNum[i], "_heatmap.pdf",sep = ""), width=15,height=8, onefile = T)
    
    par(mfrow=c(2,1))
    prDoubleBond = rbind(rbind(rbind(mtDoubleBondList[["PI"]][i,],mtDoubleBondList[["PIP"]][i,]),mtDoubleBondList[["PIP2"]][i,]),mtDoubleBondList[["PIP3"]][i,])
    rownames(prAcylLength) = pipxL
    rownames(prDoubleBond) = pipxL
    Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
              rep("brown",323))
    heatmap.2(t(prAcylLength), main = paste("acyl chain length (patient ",patientNum[i],") (T:N)",sep=""),col=redblue(256), dendrogram = "none", na.rm = TRUE, na.color=par("bg"),
              scale="column", key=T, keysize=0.5, density.info="none",
              trace="none",cexCol=1.2, RowSideColors=Label,
              lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
              lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE,Colv=FALSE)   
    heatmap.2(t(prDoubleBond), main = paste("Number of double bonds (patient ",patientNum[i],") (T:N)",sep=""),col=redblue(256), dendrogram = "none", na.rm = TRUE, na.color=par("bg"),
              scale="column", key=T, keysize=0.5, density.info="none",
              trace="none",cexCol=1.2, RowSideColors=Label,
              lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
              lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE,Colv=FALSE) 
    dev.off()  
  }
}
plotPipxAcylSpeciesHeatMap <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  pr = matrix(nrow = nrow(tData), ncol = 4)
  pipxL = c("PI","PIP","PIP2", "PIP3")
  colnames(pr) = pipxL
  rownames(pr) = rownames(tData)
  tSumV = vector()
  nSumV = vector()
  mtAcylLengthList = list()
  mtDoubleBondList = list()
  for(i in 1:4){
    tpipxVec = unlist(strsplit(colnames(tData), split ="-"))
    npipxVec = unlist(strsplit(colnames(nData), split ="-"))
    td = tData[,which(tpipxVec[seq(2,length(tpipxVec),2)] == pipxL[i])]
    nd = nData[,which(npipxVec[seq(2,length(npipxVec),2)] == pipxL[i])]
    carbonAcylList_td = vector()
    doubleBondList_td = vector()
    carbonAcylList_nd = vector()
    doubleBondList_nd = vector()
    listOfLipid = colnames(td)
    listOfCarbonLength = vector()
    carbonLengthToRemove = vector()
    doubleBondToRemove = vector()
    
    if(length(listOfLipid) > 0){
      listOfDoubleBond = vector()
      for(j in 1:length(listOfLipid)){
        acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
        listOfCarbonLength[j] = acyl[1]
        listOfDoubleBond[j] = acyl[2]
      }
      listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
      listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
      mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
      mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
      rownames(mtAcylLen) = rownames(td)
      colnames(mtAcylLen) = listOfCarbonLength
      rownames(mtDoubleBond) = rownames(td)
      colnames(mtDoubleBond) = listOfDoubleBond
      for(j in 1:length(listOfCarbonLength)){
        
        r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
        # sum all the lipid of the same acyl carbon length
        if(length(r) > 1){
          lc1 = rowSums(td[,r])/rowSums(td)
          lc2 = rowSums(nd[,r])/rowSums(nd)  
        }else{
          lc1 = td[,r]/rowSums(td)
          lc2 = nd[,r]/rowSums(nd)  
        }
        
        mtAcylLen[,j] = lc1/lc2
      }
      for(j in 1:length(listOfDoubleBond)){
        r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
        # sum all the lipid of the same double bond number
        if(length(r) > 1){
          lb1 = rowSums(td[,r])/rowSums(td)
          lb2 = rowSums(nd[,r])/rowSums(nd)  
        }else{
          lb1 = td[,r]/rowSums(td)
          lb2 = nd[,r]/rowSums(nd)  
        }
        mtDoubleBond[,j] = lb1/lb2
      }
      ind = which(mtAcylLen == 0)
      mtAcylLen = log10(mtAcylLen)
      mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
      mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
      
      ind = which(mtDoubleBond == 0)
      mtDoubleBond = log10(mtDoubleBond)
      mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
      mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
      acylTitle = "Acyl chain length"
      bondTitle = "Number of double bonds"
      
      pdf(file = paste("acyl_species_", pipxL[i],".pdf", sep = ""), width=15,height=8, onefile = T)
      
      par(mfrow=c(2,1))
      Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
                rep("brown",323))
      heatmap.2(t(mtAcylLen), main = paste("acyl chain length (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
                scale="column", key=T, keysize=0.5, density.info="none",
                trace="none",cexCol=1.2, RowSideColors=Label,
                lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
                lwid=c(1.5,0.2,2.5,2.5), srtCol=90)   
      heatmap.2(t(mtDoubleBond), main = paste("Number of double bonds (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
                scale="column", key=T, keysize=0.5, density.info="none",
                trace="none",cexCol=1.2, RowSideColors=Label,
                lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
                lwid=c(1.5,0.2,2.5,2.5), srtCol=90) 
      dev.off()
      
      
    }
  }
  
}
# This function is to identify what changes in some specific molecular species
# between 36 and 38 PIP3 groups
plotLipidSpeciesHeatMapByGroup <- function(c36TData,c36NData,c38TData,c38NData){
  lipidArray = c("Cer","TRG","DRG","MRG")
  for(i in 1:length(lipidArray)){
    td36 = filterClassPSM_1(c36TData,lipidArray[i])
    nd36 = filterClassPSM_1(c36NData,lipidArray[i])
    td38 = filterClassPSM_1(c38TData,lipidArray[i])
    nd38 = filterClassPSM_1(c38NData,lipidArray[i])
    
    pr36 = matrix(nrow = nrow(td36), ncol = ncol(td36))
    pr38 = matrix(nrow = nrow(td38), ncol = ncol(td38))
    pr36 = td36/nd36
    pr38 = td38/nd38
    colnames(pr36) = colnames(td36)
    colnames(pr38) = colnames(td38)
    ind = which(pr36 == 0)
    prMat36 = log10(pr36)
    prMat36[ind] = min(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.infinite(prMat36))] = max(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.nan(prMat36))] = 0
    nf <- layout(matrix(c(2,1),2,1, byrow=T), widths=c(1,1), heights=c(2,2), respect=T)
    ind = which(pr38 == 0)
    prMat38 = log10(pr38)
    prMat38[ind] = min(prMat36[which(is.finite(prMat38))])
    prMat38[which(is.infinite(prMat38))] = max(prMat38[which(is.finite(prMat38))])
    prMat38[which(is.nan(prMat38))] = 0
    
    mi = min(min(prMat36),min(prMat38))
    ma = max(max(prMat36),max(prMat38))
    pdf(file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/lipid_species_heat_map_",lipidArray[i],".pdf",sep=""))
    breaks = seq(mi,ma,by=(ma-mi)/100)
    aheatmap(t(prMat36), main = paste("Lipid Species Heat Map Group 36 (T:N) (", lipidArray[i],")", sep=""), Rowv = NA, breaks = breaks,Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    aheatmap(t(prMat38), main = paste("Lipid Species Heat Map Group 38 (T:N) (", lipidArray[i],")", sep=""), Rowv = NA, breaks = breaks, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    
    dev.off() 
    
    
  }
}
#This function is to plot boxplots of species for groups 36, 38 together with 
# significance levels
plotLipidSpeciesBoxPlotByGroup <- function(c36TData,c36NData,c38TData,c38NData){
  lipidArray = c("Cer","TRG","DRG","MRG")
  g <- lapply(lipidArray, function(lipidArray){
    td36 = filterClassPSM_1(c36TData,lipidArray)
    nd36 = filterClassPSM_1(c36NData,lipidArray)
    td38 = filterClassPSM_1(c38TData,lipidArray)
    nd38 = filterClassPSM_1(c38NData,lipidArray)
    
    pr36 = matrix(nrow = nrow(td36), ncol = ncol(td36))
    pr38 = matrix(nrow = nrow(td38), ncol = ncol(td38))
    pr36 = td36/nd36
    pr38 = td38/nd38
    colnames(pr36) = colnames(td36)
    colnames(pr38) = colnames(td38)
    ind = which(pr36 == 0)
    prMat36 = log10(pr36)
    prMat36[ind] = min(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.infinite(prMat36))] = max(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.nan(prMat36))] = 0
    #nf <- layout(matrix(c(2,1),2,1, byrow=T), widths=c(1,1), heights=c(2,2), respect=T)
    ind = which(pr38 == 0)
    prMat38 = log10(pr38)
    prMat38[ind] = min(prMat36[which(is.finite(prMat38))])
    prMat38[which(is.infinite(prMat38))] = max(prMat38[which(is.finite(prMat38))])
    prMat38[which(is.nan(prMat38))] = 0
    
    # create dataset
    vals = vector()
    sigL = vector()
    groups = vector()
    species = vector()
    labels = colnames(prMat36)
    n1 = nrow(prMat36)
    n2 = nrow(prMat38)
    #asterisk = vector()
    for(j in 1:ncol(prMat36)){
      vals = c(vals,c(prMat36[,j],prMat38[,j]))
      #asterisk[j] = 0.8*max(max(prMat36[,j]),max(prMat38[,j]))
      groups = c(groups, rep('36',n1), rep('38',n2))
      species = c(species, rep(labels[j], n1+n2))
      t = t.test(prMat36[,j],prMat38[,j])
      if(t$p.value < 0.05){
        sigL[j] = 'Sig.'
      }else{
        sigL[j] = ''
      }
    }
    
    dt <- data.frame(
      Amounts = vals,
      Group = groups,
      Species = species)
    label.df = data.frame(Species = colnames(td36), Amount =  as.numeric(quantile(vals)[5]) , Group = rep("36",ncol(td36)))
    
    ggplot(dt, aes(x = Species, y = Amounts, fill = Group)) +
      geom_boxplot() + geom_text(data = label.df, aes(x = colnames(td36), y = Amount, label = sigL))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))   
  })
  for(i in 1:length(lipidArray)){
    jpeg(file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/lipid_species_box_plot_",lipidArray[i],".jpeg",sep=""))
    print(g[[i]])
    dev.off()
  }
}
# This function is to identify what changes in some specific molecular species
# between group 36 and overall
plotLipidSpeciesHeatMapGroup36 <- function(c36TData,c36NData,tData,nData){
  lipidArray = c("TG","MG","DG")
  for(i in 1:length(lipidArray)){
    td36 = filterClassPSM(c36TData,lipidArray[i])
    nd36 = filterClassPSM(c36NData,lipidArray[i])
    td = filterClassPSM(tData,lipidArray[i])
    nd = filterClassPSM(nData,lipidArray[i])
    
    pr36 = matrix(nrow = nrow(td36), ncol = ncol(td36))
    pr = matrix(nrow = nrow(td), ncol = ncol(td))
    pr36 = td36/nd36
    pr = td/nd
    colnames(pr36) = colnames(td36)
    colnames(pr) = colnames(td)
    ind = which(pr36 == 0)
    prMat36 = log10(pr36)
    prMat36[ind] = min(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.infinite(prMat36))] = max(prMat36[which(is.finite(prMat36))])
    prMat36[which(is.nan(prMat36)  | is.na(prMat36) )] = 0
    nf <- layout(matrix(c(2,1),2,1, byrow=T), widths=c(1,1), heights=c(2,2), respect=T)
    ind = which(pr == 0)
    prMat = log10(pr)
    prMat[ind] = min(prMat36[which(is.finite(prMat))])
    prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
    prMat[which(is.nan(prMat) | is.na(prMat)) ] = 0
    
    mi = min(min(prMat36),min(prMat))
    ma = max(max(prMat36),max(prMat))
    pdf(file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/lipid_species_heat_map_",lipidArray[i],".pdf",sep=""), width=15,height=8, onefile = T)
    breaks = seq(mi,ma,by=(ma-mi)/100)
    aheatmap(t(prMat36), main = paste("Lipid Species Heat Map Group 36 (T:N) (", lipidArray[i],")", sep=""), Rowv = NA, breaks = breaks,Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    aheatmap(t(prMat), main = paste("Lipid Species Heat Map Overall (T:N) (", lipidArray[i],")", sep=""), Rowv = NA, breaks = breaks, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    dev.off() 
  }
}

plotAcylSpeciesHeatMapByStage <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  lipidArray = c("PIP3","PIP2","PIP","PI")
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
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
    for(i in 1:4){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      doubleBondList_td = vector()
      carbonAcylList_nd = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      carbonLengthToRemove = vector()
      doubleBondToRemove = vector()
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        for(j in 1:length(listOfLipid)){
          acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
          listOfCarbonLength[j] = acyl[1]
          listOfDoubleBond[j] = acyl[2]
        }
        listOfCarbonLength = unique(listOfCarbonLength)
        listOfDoubleBond = unique(listOfDoubleBond)
        listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
        listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
        mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
        mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
        rownames(mtAcylLen) = rownames(td)
        colnames(mtAcylLen) = listOfCarbonLength
        rownames(mtDoubleBond) = rownames(td)
        colnames(mtDoubleBond) = listOfDoubleBond
        for(j in 1:length(listOfCarbonLength)){
          
          r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
          # sum all the lipid of the same acyl carbon length
          if(length(r) > 1){
            lc1 = rowSums(td[,r])/rowSums(td)
            lc2 = rowSums(nd[,r])/rowSums(nd)  
          }else{
            lc1 = td[,r]/rowSums(td)
            lc2 = nd[,r]/rowSums(nd)  
          }
          
          mtAcylLen[,j] = lc1/lc2
        }
        for(j in 1:length(listOfDoubleBond)){
          r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
          # sum all the lipid of the same double bond number
          if(length(r) > 1){
            lb1 = rowSums(td[,r])/rowSums(td)
            lb2 = rowSums(nd[,r])/rowSums(nd)  
          }else{
            lb1 = td[,r]/rowSums(td)
            lb2 = nd[,r]/rowSums(nd)  
          }
          mtDoubleBond[,j] = lb1/lb2
        }
        
        # do t-test
        acylTitle = "Acyl chain length"
        bondTitle = "Number of double bonds"
        
        jpeg(file = paste("acyl_species", "_", lipidArray[i],"_",stages[s],".jpeg", sep = ""),width = 748, height = 779, units = "px")
        
        nf <- layout(matrix(c(1,2),1,2, byrow=T), widths=c(1,1), heights=c(2,2), respect=T)
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        #Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
        #          rep("brown",323))
        #heatmap.2(t(mtAcylLen), main = paste("acyl chain length (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="col",na.rm = TRUE, na.color=par("bg"),
        #          scale="column", key=T, keysize=0.5, density.info="none",
        #          trace="none",cexCol=1.2, RowSideColors=Label,
        #          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
        #          lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE)   
        #heatmap.2(t(mtDoubleBond), main = paste("Number of double bonds (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="col",na.rm = TRUE, na.color=par("bg"),
        #          scale="column", key=T, keysize=0.5, density.info="none",
        #          trace="none",cexCol=1.2, RowSideColors=Label,
        #          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
        #          lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE) 
        aheatmap(t(mtAcylLen), main = paste("acyl chain length (",lipidArray[i],") (T:N) ",stages[s],sep=""), Rowv = NA, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
        aheatmap(t(mtDoubleBond), main = paste("Number of double bonds (",lipidArray[i],") (T:N) ",stages[s],sep=""),Rowv = NA, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL) 
        dev.off()
        
      }
    }
  }
}
# This function is to produce acyl chain length heatmaps for phospholipid and DG
PlotAcylchainPhosphoLipidHeatMapByStage <- function(tData, nData){
  
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/Phospholipid heatmaps/acyl chain")
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[which(grepl("P",lipidArray))] # extract only phospholipid
  lipidArray = c(lipidArray,"DRG","PIP3","PIP2","PIP")
  
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste(stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:20),5,4, byrow=T), respect=T)
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      carbonAcylList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfCarbonLength[1] = "<34"
          mtAcylLen = td/nd  
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = listOfCarbonLength
        }
        else{
          for(j in 1:length(listOfLipid)){
            acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
            listOfCarbonLength[j] = acyl[1]  
          }
          listOfCarbonLength = as.numeric(listOfCarbonLength)
          rr = list()
          rr[[1]] = which(listOfCarbonLength <= 32)
          rr[[2]] = which(listOfCarbonLength >= 34 & listOfCarbonLength <= 36)
          rr[[3]] = which(listOfCarbonLength > 36)
          r = vector()
          if(length(rr[[1]]) > 0){
            r["<34"] = 1
          }else{
            r["<34"] = 0
          }
          if(length(rr[[2]]) > 0){
            r["34-36"] = 1
          }else{
            r["34-36"] = 0
          }
          if(length(rr[[3]]) > 0){
            r[">36"] = 1
          }else{
            r[">36"] = 0
          }
          
          name = names(r)
          id = which(r > 0)
          name = name[id]
          rr = rr[id]
          mtAcylLen = matrix(nrow = nrow(td), ncol = length(name))
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = name
          
          for(j in 1:length(name)){
            if(length(rr[[j]]) > 1){
              lb1 =  rowSums(td[,rr[[j]]])/rowSums(td)
              lb2 = rowSums(nd[,rr[[j]]])/rowSums(nd)  
            }else{
              lb1 = td[,rr[[j]]]/rowSums(td)
              lb2 = nd[,rr[[j]]]/rowSums(nd)
            }
            
            mtAcylLen[,j] = lb1/lb2
          }
        }
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        if(all(mtAcylLen==0 | is.infinite(mtAcylLen)))
          next
        aheatmap(t(mtAcylLen), main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = NA, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)       
      }
    }
    dev.off()
  }
  graphics.off()
  
}
# This function is to produce double bond heatmaps for phospholipid and DG
PlotDoubleBondPhosphoLipidHeatMapByStage <- function(tData, nData){
  
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/Phospholipid heatmaps/double bond")
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[which(grepl("P",lipidArray))] # extract only phospholipid
  lipidArray = c(lipidArray,"DRG","PIP3","PIP2","PIP")
  
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste(stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:20),5,4, byrow=T), respect=T)
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      doubleBondList_td = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfDoubleBond = vector()
      
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfDoubleBond[1] = "1,2,3"
          mtDoubleBond = td/nd  
          rownames(mtDoubleBond) = rownames(td)
          colnames(mtDoubleBond) = listOfDoubleBond
        }
        else{
          for(j in 1:length(listOfLipid)){
            acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
            listOfDoubleBond[j] = acyl[2]  
          }
          
          listOfDoubleBond = as.numeric(listOfDoubleBond)
          rr = list()
          rr[[1]] = which(listOfDoubleBond == 0)
          rr[[2]] = which(listOfDoubleBond > 0 & listOfDoubleBond < 4)
          rr[[3]] = which(listOfDoubleBond > 3)
          r = vector()
          if(length(rr[[1]]) > 0){
            r["0"] = 1
          }else{
            r["0"] = 0
          }
          if(length(rr[[2]]) > 0){
            r["1,2,3"] = 1
          }else{
            r["1,2,3"] = 0
          }
          if(length(rr[[3]]) > 0){
            r[">4"] = 1
          }else{
            r[">4"] = 0
          }
          
          name = names(r)
          id = which(r > 0)
          name = name[id]
          rr = rr[id]
          mtDoubleBond = matrix(nrow = nrow(td), ncol = length(name))
          rownames(mtDoubleBond) = rownames(td)
          colnames(mtDoubleBond) = name
          for(j in 1:length(name)){
            if(length(rr[[j]]) > 1){
              lb1 = rowSums(td[,rr[[j]]])/rowSums(td)
              lb2 = rowSums(nd[,rr[[j]]])/rowSums(nd)  
            }else{
              lb1 = td[,rr[[j]]]/rowSums(td)
              lb2 = nd[,rr[[j]]]/rowSums(nd)
            }
            
            mtDoubleBond[,j] = lb1/lb2
          }
        }
        
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        if(all(mtDoubleBond==0 | is.infinite(mtDoubleBond)))
          next
        aheatmap(t(mtDoubleBond), main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = NA, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
        
      }
    }
    dev.off()
  }
  graphics.off()
  
}
# This function is to find the correlation between PIP3 and other lipid classes
# Change in double bonds
plotDoubleBondLipidClassHeatMapByStage <- function(tData,nData){
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIP3 correlation/double bond")
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[-which(lipidArray == "CH")]
  lipidArray = c(lipidArray,"PIP3","PIP2","PIP")
  id = which(lipidArray == "TG")
  lipidArray[id] = "TRG"
  id = which(lipidArray == "DG")
  lipidArray[id] = "DRG"
  id = which(lipidArray == "MG")
  lipidArray[id] = "MRG"
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste(stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:28),7,4, byrow=T), respect=T)
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      doubleBondList_td = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfDoubleBond = vector()
      
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfDoubleBond[1] = "1"
          mtDoubleBond = td/nd  
          rownames(mtDoubleBond) = rownames(td)
          colnames(mtDoubleBond) = listOfDoubleBond
        }
        else{
          for(j in 1:length(listOfLipid)){
            acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
            listOfDoubleBond[j] = acyl[2]  
          }
          listOfDoubleBond = unique(listOfDoubleBond)
          listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
          mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
          rownames(mtDoubleBond) = rownames(td)
          colnames(mtDoubleBond) = listOfDoubleBond
          
          for(j in 1:length(listOfDoubleBond)){
            
            r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
            # sum all the lipid of the same double bond number
            if(length(r) > 1){
              lb1 = rowSums(td[,r])/rowSums(td)
              lb2 = rowSums(nd[,r])/rowSums(nd)  
            }else{
              lb1 = td[,r]/rowSums(td)
              lb2 = nd[,r]/rowSums(nd)  
            }
            mtDoubleBond[,j] = lb1/lb2
          }
        }
        
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        if(all(mtDoubleBond==0 | is.infinite(mtDoubleBond)))
          next
        aheatmap(t(mtDoubleBond), main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
        
      }
    }
    dev.off()
  }
  graphics.off()
}
# This function is to compute the correlation between PIPx and other lipid species
# double bond
plotDoubleBondCorr <- function(tData, nData){
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIP3 correlation/double bond")
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[-which(lipidArray == "CH")]
  lipidArray = c(lipidArray,"PIP3","PIP2","PIP")
  id = which(lipidArray == "TG")
  lipidArray[id] = "TRG"
  id = which(lipidArray == "DG")
  lipidArray[id] = "DRG"
  id = which(lipidArray == "MG")
  lipidArray[id] = "MRG"
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(file = paste("double_bond_species_cluster (",stages[s],").pdf", sep = "")) 
    dataList = list()
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      doubleBondList_td = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfDoubleBond = vector()
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfDoubleBond[1] = "1"
          mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
          rownames(mtDoubleBond) = rownames(td)
          colnames(mtDoubleBond) = paste(listOfDoubleBond,"-",lipidArray[i],sep="")
          mtDoubleBond = td/nd  
        }
        else{
          if(lipidArray[i] == "PIP3" | lipidArray[i] == "PIP2" | lipidArray[i] == "PIP"){  
            listOfDoubleBond = c("1,2","3,4")
            mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
            rownames(mtDoubleBond) = rownames(td)
            colnames(mtDoubleBond) = paste(listOfDoubleBond,"-",lipidArray[i],sep="")
            r1 = which(grepl(":1",listOfLipid))
            r2 = which(grepl(":2",listOfLipid))
            r3 = which(grepl(":3",listOfLipid))
            r4 = which(grepl(":4",listOfLipid))
            # sum all the lipid of the same double bond number
            r12 = c(r1,r2)
            r34 = c(r3,r4)
            if(length(r12) > 1){
              lb1 = rowSums(td[,r12])/rowSums(td)
              lb2 = rowSums(nd[,r12])/rowSums(nd)  
            }else{
              lb1 = td[,r12]/rowSums(td)
              lb2 = nd[,r12]/rowSums(nd)  
            }
            
            if(length(r34) > 1){
              lbb1 = rowSums(td[,r34])/rowSums(td)
              lbb2 = rowSums(nd[,r34])/rowSums(nd)  
            }else{
              lbb1 = td[,r34]/rowSums(td)
              lbb2 = nd[,r34]/rowSums(nd)  
            }
            
            mtDoubleBond[,1] = lb1/lb2
            mtDoubleBond[,2] = lbb1/lbb2
            
          }else{
            for(j in 1:length(listOfLipid)){
              acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
              listOfDoubleBond[j] = acyl[2]  
            }
            listOfDoubleBond = unique(listOfDoubleBond)
            listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
            mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
            rownames(mtDoubleBond) = rownames(td)
            colnames(mtDoubleBond) = paste(listOfDoubleBond,"-",lipidArray[i],sep="")
            for(j in 1:length(listOfDoubleBond)){
              
              r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
              # sum all the lipid of the same double bond number
              if(length(r) > 1){
                lb1 = rowSums(td[,r])/rowSums(td)
                lb2 = rowSums(nd[,r])/rowSums(nd)  
              }else{
                lb1 = td[,r]/rowSums(td)
                lb2 = nd[,r]/rowSums(nd)  
              }
              mtDoubleBond[,j] = lb1/lb2
            }
            
          }
          
          
        }
        
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        mtDoubleBond = mtDoubleBond[,which(!apply(mtDoubleBond==0 | is.infinite(mtDoubleBond) ,2,all)),drop=F]
        dataList[[lipidArray[i]]] = mtDoubleBond
      }
    }
    mtData = dataList[[1]]
    for(i in 2:length(dataList)){
      mtData = cbind(mtData,dataList[[i]])
    }
    #par(mfrow=c(2,1))
    # do PCA
    #my.active = mtData
    #matrix = as.numeric(as.matrix(my.active))
    #dim(matrix) = dim(my.active)
    #rownames(matrix) = rownames(my.active)
    mtData = mtData[,apply(mtData, 2, var, na.rm=TRUE) != 0] # remove zero variance columns
    res.pca = prcomp(mtData, center = TRUE, scale = TRUE)
    extractCorrCycleQuadrant(stages[s],res.pca,1)
    #the matrix of variable loadings (columns are eigenvectors)
    #head(unclass(res.pca$rotation)[, 1:4])
    # Eigenvalues
    #eig = (res.pca$sdev)^2
    
    # Variances in percentage
    #variance = eig*100/sum(eig)
    
    # Cumulative variances
    #cumvar = cumsum(variance)
    
    #eig.matrix = data.frame(eig = eig, variance = variance,
    #                        cumvariance = cumvar)
    #head(eig.matrix)
    #eig.val = get_eigenvalue(res.pca)
    #head(eig.val)
    
    # Variable correlation/coordinates
    # fviz_pca_var(res.pca)
    #loadings = res.pca$rotation
    #sdev = res.pca$sdev
    #var <- get_pca_var(res.pca)
    #var.coord = var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
    #head(var.coord[, 1:4])
    
    # Plot the correlation circle
    #a <- seq(0, 2*pi, length = 100)
    #plot( cos(a), sin(a), type = 'l', col="gray",
    #      xlab = "PC1",  ylab = "PC2")
    
    #abline(h = 0, v = 0, lty = 2)
    
    # Add active variables
    #arrows(0, 0, var.coord[, 1], var.coord[, 2], 
    #       length = 0.1, angle = 15, code = 2)
    
    # Add labels
    #text(var.coord, labels=colnames(my.active), cex = 0.5, adj=1)
    
    # clustering based on the absolute correlation between variables
    mtData = mtData[,which(!apply(mtData==0,2,all))]
    zeroSd = which(apply(mtData,2,sd)==0)
    if(length(zeroSd) > 0)
      mtData = mtData[,-zeroSd]
    library(psych)
    corr = cor(na.omit(mtData))
    # plot upper correlation matrix
    #corrplot(corr, type="upper", order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.5)
    M = hclust(dist(abs(corr)))
    M$call = NULL
    plot(M, xlab = "Double bond species", cex= 0.3)
    # cut the tree (k: number of clusters)
    #rect.hclust(M, k=5, border="red")
    
  }
  graphics.off()
}

plotPIPxDoubleBondLipidClassHeatMapStage <- function(tData, nData){
  lipidArray = c("PIP3","PIP2","PIP","PI")
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIPx heatmaps/double_bond_",stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:4),2,2, byrow=T), respect=T)
    mtList = list()
    mi = vector()
    ma = vector()
    for(i in 1:length(lipidArray)){
      td = filterClassPSM_1(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM_1(nData[patientList[[s]],],lipidArray[i])
      doubleBondList_td = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfDoubleBond = vector()
      
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        
        for(j in 1:length(listOfLipid)){
          acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
          listOfDoubleBond[j] = acyl[2]  
        }
        listOfDoubleBond = unique(listOfDoubleBond)
        listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
        mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
        rownames(mtDoubleBond) = rownames(td)
        colnames(mtDoubleBond) = listOfDoubleBond
        
        for(j in 1:length(listOfDoubleBond)){
          
          r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
          # sum all the lipid of the same double bond number
          if(length(r) > 1){
            lb1 = rowSums(td[,r])/rowSums(td)
            lb2 = rowSums(nd[,r])/rowSums(nd)  
          }else{
            lb1 = td[,r]/rowSums(td)
            lb2 = nd[,r]/rowSums(nd)  
          }
          mtDoubleBond[,j] = lb1/lb2
        }
        
        
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        if(all(mtDoubleBond==0 | is.infinite(mtDoubleBond)))
          next
        #ind = which(mtDoubleBond > 2, arr.ind = TRUE) # identify the element which messes up the scale
        #if(length(ind) > 0){
        #  mtDoubleBond[ind[1],ind[2]] = max(mtDoubleBond[-ind[1],ind[2]])
        #}
        mi[i] = min(mtDoubleBond)
        ma[i] = max(mtDoubleBond)
        mtList[[i]] = t(mtDoubleBond)
        
        
      }
      
    }
    mmi = min(mi)
    mma = max(ma)
    breaks = seq(mmi,mma,by=(mma-mmi)/100)
    for(i in 1:length(lipidArray)){
      aheatmap(mtList[[i]], main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = NA, breaks = breaks,revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    }
    dev.off()
  }
  graphics.off()
  
}
plotPIPxAcylChainLipidClassHeatMapStage <- function(tData, nData){
  lipidArray = c("PIP3","PIP2","PIP","PI")
  
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIPx heatmaps/acyl_chain_",stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:4),2,2, byrow=T), respect=T)
    mi = vector()
    ma = vector()
    mtList = list()
    for(i in 1:length(lipidArray)){
      td = filterClassPSM_1(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM_1(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      carbonAcylList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      
      if(length(listOfLipid) > 0){
        
        
        
        listOfDoubleBond = vector()
        
        
        for(j in 1:length(listOfLipid)){
          acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
          listOfCarbonLength[j] = acyl[1]  
        }
        listOfCarbonLength = unique(listOfCarbonLength)
        listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
        mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
        rownames(mtAcylLen) = rownames(td)
        colnames(mtAcylLen) = listOfCarbonLength
        
        for(j in 1:length(listOfCarbonLength)){
          
          r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
          # sum all the lipid of the same acyl carbon length
          if(length(r) > 1){
            lc1 = rowSums(td[,r])
            lc2 = rowSums(nd[,r]) 
          }else{
            lc1 = td[,r]
            lc2 = nd[,r]
          }
          
          mtAcylLen[,j] = lc1/lc2
        }
        
        
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        if(all(mtAcylLen==0 | is.infinite(mtAcylLen)))
          next
        mi[i] = min(mtAcylLen)
        ma[i] = max(mtAcylLen)
        mtList[[i]] = t(mtAcylLen)
      }
    }
    mmi = min(mi)
    mma = max(ma)
    breaks = seq(mmi,mma,by=(mma-mmi)/100)
    for(i in 1:length(lipidArray)){
      aheatmap(mtList[[i]], main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = NA, breaks = breaks, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
    }
    dev.off()
  }
  graphics.off()
  
}
plotAcylChainLipidClassHeatMapByStage <- function(tData, nData){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[-which(lipidArray == "CH")]
  lipidArray = c(lipidArray,"PIP3","PIP2","PIP")
  id = which(lipidArray == "TG")
  lipidArray[id] = "TRG"
  id = which(lipidArray == "DG")
  lipidArray[id] = "DRG"
  id = which(lipidArray == "MG")
  lipidArray[id] = "MRG"
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  #Set pdf dimensions
  if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
  if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(width= w, height=h, file = paste("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIP3 correlation/acyl chain/",stages[s],".pdf", sep = ""))  
    nf <- layout(matrix(seq(1:28),7,4, byrow=T), respect=T)
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      carbonAcylList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      
      if(length(listOfLipid) > 0){
        
        
        
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfCarbonLength[1] = "18"
          mtAcylLen = td/nd  
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = listOfCarbonLength
        }
        else{
          for(j in 1:length(listOfLipid)){
            acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
            listOfCarbonLength[j] = acyl[1]  
          }
          listOfCarbonLength = unique(listOfCarbonLength)
          listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
          mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = listOfCarbonLength
          
          for(j in 1:length(listOfCarbonLength)){
            
            r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
            # sum all the lipid of the same acyl carbon length
            if(length(r) > 1){
              lc1 = rowSums(td[,r])/rowSums(td)
              lc2 = rowSums(nd[,r])/rowSums(nd)  
            }else{
              lc1 = td[,r]/rowSums(td)
              lc2 = nd[,r]/rowSums(nd)  
            }
            
            mtAcylLen[,j] = lc1/lc2
          }
        }
        
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        if(all(mtAcylLen==0 | is.infinite(mtAcylLen)))
          next
        aheatmap(t(mtAcylLen), main = paste(lipidArray[i]," (T:N) ",stages[s],sep=""), Rowv = NA, Colv = TRUE, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
        
        #Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
        #          rep("brown",323))
        #heatmap.2(t(mtAcylLen), main = paste("acyl chain length (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="col",na.rm = TRUE, na.color=par("bg"),
        #          scale="column", key=T, keysize=0.5, density.info="none",
        #          trace="none",cexCol=1.2, RowSideColors=Label,
        #          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
        #          lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE)   
        #heatmap.2(t(mtDoubleBond), main = paste("Number of double bonds (",lipidArray[i],") (T:N) ",stages[s],sep=""),col=redblue(256), dendrogram="col",na.rm = TRUE, na.color=par("bg"),
        #          scale="column", key=T, keysize=0.5, density.info="none",
        #          trace="none",cexCol=1.2, RowSideColors=Label,
        #          lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
        #          lwid=c(1.5,0.2,2.5,2.5), srtCol=90,Rowv=FALSE) 
        
      }
    }
    dev.off()
  }
  graphics.off()
}
# This function is to identify lipid species in a correlation circle
extractCorrCycleQuadrant <- function(stage, res_pca,type){
  fviz_pca_var(res_pca)
  loadings = res_pca$rotation
  sdev = res_pca$sdev
  var <- get_pca_var(res_pca)
  var.coord = var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
  pc1 = var.coord[,1]
  pc2 = var.coord[,2]
  species = rownames(var.coord)
  quadrant = list()
  quadrant[[1]] = species[which(pc1 > 0 & pc2 > 0)]
  quadrant[[2]] = species[which(pc1 < 0  & pc2 > 0)]
  quadrant[[3]] = species[which(pc1 < 0 & pc2 < 0)]
  quadrant[[4]] = species[which(pc1 > 0 & pc2 < 0)]
  l[1] = length(quadrant[[1]])
  l[2] = length(quadrant[[2]])
  l[3] = length(quadrant[[3]])
  l[4] = length(quadrant[[4]])
  ll = max(l[1],l[2],l[3],l[4])
  quadrantData = matrix(nrow = ll, ncol = 4)
  for(i in 1:4){
    for(j in 1:l[i]){
      quadrantData[j,i] = quadrant[[i]][j]
    }
  }
  if(type == 0){
    write.xlsx(quadrantData, paste(stage,"_acyl_chain.xlsx",sep=""))  
  }else{
    write.xlsx(quadrantData, paste(stage,"_double_bond.xlsx",sep=""))  
  }
  
}
# This function is to compute the correlation between PIPx and other lipid species
plotAcylChainCorr <- function(tData, nData){
  setwd("D:/project/lipidomics/data/lipid analysis/pipx analysis/PIP3 correlation/acyl chain")
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = lipidArray[-which(lipidArray == "CH")]
  lipidArray = c(lipidArray,"PIP3","PIP2","PIP")
  id = which(lipidArray == "TG")
  lipidArray[id] = "TRG"
  id = which(lipidArray == "DG")
  lipidArray[id] = "DRG"
  id = which(lipidArray == "MG")
  lipidArray[id] = "MRG"
  numberOfSpecies =  length(lipidArray)/2
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  patientList = list()
  
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(s in 1:length(stages)){
    pdf(file = paste("acyl_species_cluster (",stages[s],").pdf", sep = "")) 
    dataList = list()
    for(i in 1:length(lipidArray)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      carbonAcylList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      carbonLengthToRemove = vector()
      
      if(length(listOfLipid) > 0){
        
        
        
        listOfDoubleBond = vector()
        if(length(listOfLipid) == 1){
          listOfCarbonLength[1] = "18"
          mtAcylLen = td/nd  
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = paste(listOfCarbonLength,"-",lipidArray[i],sep="")
        }
        else{
          for(j in 1:length(listOfLipid)){
            acyl = unlist(strsplit(unlist(strsplit(listOfLipid[j], split ="-"))[1], split = ":"))
            listOfCarbonLength[j] = acyl[1]  
          }
          listOfCarbonLength = unique(listOfCarbonLength)
          listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
          mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
          rownames(mtAcylLen) = rownames(td)
          colnames(mtAcylLen) = paste(listOfCarbonLength,"-",lipidArray[i],sep="")
          
          for(j in 1:length(listOfCarbonLength)){
            
            r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
            # sum all the lipid of the same acyl carbon length
            if(length(r) > 1){
              lc1 = rowSums(td[,r])/rowSums(td)
              lc2 = rowSums(nd[,r])/rowSums(nd)  
            }else{
              lc1 = td[,r]/rowSums(td)
              lc2 = nd[,r]/rowSums(nd)  
            }
            
            mtAcylLen[,j] = lc1/lc2
          }
        }
        
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        mtAcylLen = mtAcylLen[,which(!apply(mtAcylLen==0 | is.infinite(mtAcylLen) ,2,all)),drop=F]
        dataList[[lipidArray[i]]] = mtAcylLen
      }
    }
    mtData = dataList[[1]]
    for(i in 2:length(dataList)){
      mtData = cbind(mtData,dataList[[i]])
    }
    #par(mfrow=c(2,1))
    # do PCA
    #my.active = mtData
    #matrix = as.numeric(as.matrix(my.active))
    #dim(matrix) = dim(my.active)
    #rownames(matrix) = rownames(my.active)
    mtData = mtData[,apply(mtData, 2, var, na.rm=TRUE) != 0] # remove zero variance columns
    res.pca = prcomp(mtData, center = TRUE, scale = TRUE)
    extractCorrCycleQuadrant(stages[s],res.pca,0)
    #the matrix of variable loadings (columns are eigenvectors)
    #head(unclass(res.pca$rotation)[, 1:4])
    # Eigenvalues
    #eig = (res.pca$sdev)^2
    
    # Variances in percentage
    #variance = eig*100/sum(eig)
    
    # Cumulative variances
    #cumvar = cumsum(variance)
    
    #eig.matrix = data.frame(eig = eig, variance = variance,
    #                        cumvariance = cumvar)
    #head(eig.matrix)
    #eig.val = get_eigenvalue(res.pca)
    #head(eig.val)
    
    # Variable correlation/coordinates
    #fviz_pca_var(res.pca)
    #loadings = res.pca$rotation
    #sdev = res.pca$sdev
    #var <- get_pca_var(res.pca)
    #var.coord = var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
    #head(var.coord[, 1:4])
    
    # Plot the correlation circle
    #a <- seq(0, 2*pi, length = 100)
    #plot( cos(a), sin(a), type = 'l', col="gray",
    #      xlab = "PC1",  ylab = "PC2")
    
    #abline(h = 0, v = 0, lty = 2)
    
    # Add active variables
    #arrows(0, 0, var.coord[, 1], var.coord[, 2], 
    #       length = 0.1, angle = 15, code = 2)
    
    # Add labels
    #text(var.coord, labels=colnames(my.active), cex = 0.5, adj=1)
    
    # clustering based on the absolute correlation between variables
    mtData = mtData[,which(!apply(mtData==0,2,all))]
    library(psych)
    corr = cor(na.omit(mtData))
    # plot upper correlation matrix
    #corrplot(corr, type="upper", order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.5)
    M = hclust(dist(abs(corr)))
    M$call = NULL
    plot(M, xlab = "Acyl species", cex= 0.3)
  }
  graphics.off()
}
plotAcylSpeciesPatientHeatMapByStage <- function(tData, nData){
  tData = normalisePatientSpecieMatrix(tData)
  nData = normalisePatientSpecieMatrix(nData)
  lipidArray = c("PI","PIP","PIP2","PIP3")
  patientData = getPatientData()
  patientData = patientData[rownames(tData),]
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
  acylLenStageList = list()
  doubleBondStageList = list()
  for(s in 1:length(stages)){
    acylLenStageList[[stages[s]]] = list()
    doubleBondStageList[[stages[s]]] = list()
    for(i in 1:4){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      carbonAcylList_td = vector()
      doubleBondList_td = vector()
      carbonAcylList_nd = vector()
      doubleBondList_nd = vector()
      listOfLipid = colnames(td)
      listOfCarbonLength = vector()
      carbonLengthToRemove = vector()
      doubleBondToRemove = vector()
      
      if(length(listOfLipid) > 0){
        listOfDoubleBond = vector()
        listOfCarbonLength = c("36","38")
        listOfDoubleBond = c("1","2","3","4")
        listOfCarbonLength = as.character(sort(as.numeric(listOfCarbonLength)))
        listOfDoubleBond = as.character(sort(as.numeric(listOfDoubleBond)))
        mtAcylLen = matrix(nrow = nrow(td), ncol = length(listOfCarbonLength))
        mtDoubleBond = matrix(nrow = nrow(td), ncol = length(listOfDoubleBond))
        
        rownames(mtAcylLen) = rownames(td)
        colnames(mtAcylLen) = listOfCarbonLength
        rownames(mtDoubleBond) = rownames(td)
        colnames(mtDoubleBond) = listOfDoubleBond
        for(j in 1:length(listOfCarbonLength)){
          
          r = which(grepl(paste(listOfCarbonLength[j],":", sep=""),listOfLipid))
          # sum all the lipid of the same acyl carbon length
          if(length(r) > 1){
            lc1 = rowSums(td[,r])/rowSums(td)
            lc2 = rowSums(nd[,r])/rowSums(nd)  
          }else{
            lc1 = td[,r]/rowSums(td)
            lc2 = nd[,r]/rowSums(nd)  
          }
          
          mtAcylLen[,j] = lc1/lc2
        }
        for(j in 1:length(listOfDoubleBond)){
          r = which(grepl(paste(":",listOfDoubleBond[j], sep = ""),listOfLipid))
          # sum all the lipid of the same double bond number
          if(length(r) > 1){
            lb1 = rowSums(td[,r])/rowSums(td)
            lb2 = rowSums(nd[,r])/rowSums(nd)  
          }else{
            lb1 = td[,r]/rowSums(td)
            lb2 = nd[,r]/rowSums(nd)  
          }
          mtDoubleBond[,j] = lb1/lb2
        }
        
        ind = which(mtAcylLen == 0)
        mtAcylLen = log10(mtAcylLen)
        mtAcylLen[ind] = min(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.infinite(mtAcylLen))] = max(mtAcylLen[which(is.finite(mtAcylLen))])
        mtAcylLen[which(is.nan(mtAcylLen))] = 0
        ind = which(mtDoubleBond == 0)
        mtDoubleBond = log10(mtDoubleBond)
        mtDoubleBond[ind] = min(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.infinite(mtDoubleBond))] = max(mtDoubleBond[which(is.finite(mtDoubleBond))])
        mtDoubleBond[which(is.nan(mtDoubleBond))] = 0
        
        acylLenStageList[[stages[s]]][[lipidArray[i]]] = mtAcylLen
        doubleBondStageList[[stages[s]]][[lipidArray[i]]] = mtDoubleBond
        
      }
    }
  }
  
  stages = c("adenoma","dukes a","dukes b","dukes c","dukes d")
  for(s in 1:length(stages)){
    patientNum = rownames(acylLenStageList[[stages[s]]][["PI"]])
    len = length(patientNum)
    if (len <5) w<- len * 8 else w<- 32
    if (len <5) h<- 8 else h<- (((len/4)+1)*8)
    pdf(file = paste("acyl_species_patient_by_stage_heatmap_", stages[s], ".pdf",sep = ""), width= w, height=h)
    #Set graphcial grid
    par(mfrow = c(len,2))
    
    for(i in 1:length(patientNum)){
      prAcylLength = matrix(nrow = 4, ncol = 2)
      prDoubleBond = matrix(nrow = 4, ncol = 4)
      prAcylLength = rbind(rbind(rbind(acylLenStageList[[stages[s]]][["PI"]][i,],acylLenStageList[[stages[s]]][["PIP"]][i,]),acylLenStageList[[stages[s]]][["PIP2"]][i,]),acylLenStageList[[stages[s]]][["PIP3"]][i,])
      prDoubleBond = rbind(rbind(rbind(doubleBondStageList[[stages[s]]][["PI"]][i,],doubleBondStageList[[stages[s]]][["PIP"]][i,]),doubleBondStageList[[stages[s]]][["PIP2"]][i,]),doubleBondStageList[[stages[s]]][["PIP3"]][i,])
      rownames(prAcylLength) <- lipidArray
      rownames(prDoubleBond) <- lipidArray
      aheatmap(t(prAcylLength), main = paste("acyl chain length (T:N) patient (",patientNum[i],")",sep=""), Rowv = NA, Colv = NA, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL)   
      aheatmap(t(prDoubleBond), main = paste("Number of double bonds (T:N) patient (",patientNum[i],")",sep=""),Rowv = NA, Colv = NA, revC = FALSE,annLegend = TRUE, labRow = NULL, labCol = NULL) 
    }
    dev.off()
  }
  graphics.off()
}