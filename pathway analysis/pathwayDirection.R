getReactionRatio <- function(PSM){
  reactions = as.matrix(read.csv("D:/project/lipidomics/reactions.csv", header = T))
  result = list()
  reactionCount = list()
  layerCount = list()
  sReactions = matrix(ncol = 2)
  firstRun = TRUE
  
  for(i in 1:nrow(reactions)){
    reactant = reactions[i,1]
    product = reactions[i,2]  
    rSpecies = colnames(filterClassPSM(PSM,reactant))
    pSpecies = colnames(filterClassPSM(PSM,product))
    if(reactant == "DG"){
      reactant = "DRG"
    }
    
    if(product == "DG"){
      product = "DRG"
    }
    rFAs = getFAs(rSpecies)
    pFAs =  getFAs(pSpecies)
    sharedFAs = which(rFAs %in% pFAs)
    for(j in 1:length(sharedFAs)){  
      fai = sharedFAs[j]
      fa = rFAs[fai]
      if(is.null(layerCount[[fa]]))		{
        layerCount[[fa]] = 1
      }else{
        layerCount[[fa]] = layerCount[[fa]] + 1
      }
      if(is.null(reactionCount[[paste(reactions[i,], collapse = "->")]])){
        reactionCount[[paste(reactions[i,], collapse = "->")]] = 1
      }else{
        reactionCount[[paste(reactions[i,], collapse = "->")]] = reactionCount[[paste(reactions[i,], collapse = "->")]] + 1
      }
      if(firstRun){
        
        sReactions[1,] = c(paste(fa, reactant,sep = "-"),paste(pFAs[which(pFAs==fa)],product,sep = "-"))
        firstRun = FALSE
      }else{
        sReactions = rbind(sReactions, c(paste(fa,reactant,sep = "-"),paste(pFAs[which(pFAs==fa)],product,sep = "-")))
        
      }
    }
  }
  #par(ask = TRUE)
  prMat = matrix(ncol =nrow(sReactions), nrow = nrow(PSM))	
  colVec = vector(length = nrow(sReactions))
  rownames(prMat) = rownames(PSM)
  for(i in 1:nrow(sReactions)){
    reactant = sReactions[i,1]
    product = sReactions[i,2]
    #colVec[i] =  paste(reactant,"->", product, sep = "")
    s = unlist(strsplit(reactant, split="-"))
    layer = s[1]
    rec = s[2]
    s1 = unlist(strsplit(product, split="-"))
    pro = s1[2]
    colVec[i] =  paste(paste(rec,"-", pro, sep = ""), "[", layer, "]", sep = "")
    rVec = PSM[,which(colnames(PSM) == reactant)]
    pVec = PSM[,which(colnames(PSM)== product)]
    prVec = vector(length = length(pVec))
    #for(j in 1:length(pVec)){
    #  if(rVec[j] != 0 && pVec[j] != 0){
    #    prVec[j] = log10(pVec[j]/rVec[j]);
    #  }else{
    #    prVec[j] = 4;
    #  }
    #}
    prVec = pVec/rVec
    prMat[,i] =prVec
  }
  colnames(prMat) = colVec	
  result$prMat =  prMat
  result$sReactions = sReactions
  result$reactionCount = reactionCount
  result$layerCount = layerCount
  result
}