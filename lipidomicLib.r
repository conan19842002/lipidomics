#a library of regularly used functions in the babraham project.
library(RMySQL)
library(MASS)
library(caret)
library(gplots)
library(RColorBrewer)
library(igraph)
library(stringr)
#con <- dbConnect(MySQL(), user="admin", password="johnny5", dbname= "crc", host = "mysql-jfoster.ebi.ac.uk", port = 4226)
con <- dbConnect(MySQL(), user="root", password="root", dbname= "lnet", host = "127.0.0.1", port = 3306)

getPatientData <-function(){
  patientQuery = "select p.number, s.name, p.date_of_analysis, p.tumor_sample_amount_microgram, p.normal_sample_amount_microgram from patient as p, stage as s where p.l_stage_id = s.stage_id order by 1"
  patData = dbGetQuery(con, patientQuery)
  rownames(patData) = patData[,1]
  patData = patData[,-1]
  patData
}
# This function is to select more fields in patient table for PCA analysis
getAllPatientData <- function(){
  patientQuery = "select p.number, p.genre, p.age, s.name, p.date_of_analysis, p.tumor_sample_amount_microgram, p.normal_sample_amount_microgram from patient as p, stage as s where p.l_stage_id = s.stage_id and p.age IS NOT NULL order by 1"
  patData = dbGetQuery(con, patientQuery)
  rownames(patData) = patData[,1]
  patData = patData[,-1]
  patData
}
getLipidClass <- function(){
  patientQuery = "select name from lipid_class order by name"
  patData = dbGetQuery(con, patientQuery)
  #rownames(patData) = patData[,1]
  patData = patData[,1]
  patData
}
getPatientStageData <- function(){
  patientQuery = "select s.name, count(p.number) from patient as p, stage as s where p.l_stage_id = s.stage_id group by s.name"
  patData = dbGetQuery(con, patientQuery)
  rownames(patData) = patData[,1]
  patData
}
getPatientFattyAcidData <- function(){
  patientQuery = "select DISTINCT(p.number), s.name from patient_has_fatty_acid_species pt, stage s, patient p where pt.f_patient_id = p.patient_id and p.l_stage_id = s.stage_id order by 1"
  patData = dbGetQuery(con, patientQuery)
  rownames(patData) = patData[,1]
  patData
}
# tumour is dataset type
# true : tumour dataset ( TRG, DRG, MRG) -> (TG, DG, MG)
# false : other dataset ( TG, DG, MG)
getMassSpecies <- function(speciesL, tumour){
  mass = rep(1,length(speciesL))
  names(mass) = speciesL
  id1 = which(speciesL == "SG" | speciesL == "C18-SG" | speciesL == "C18 Sphingosine")
  if(length(id1) > 0){
    mass[id1] = 299.2824
  }
  
  id2 = which(speciesL == "C18-S1P")
  if(length(id2) > 0){
    mass[id2] = 379.472
  }
  id3 = which(speciesL == "CH")
  if(length(id3) > 0){
    mass[id3] = 1
  }
  id4 = which(speciesL == "C18-SPC")
  if(length(id4) > 0){
    mass[id4] = 465.3457
  }
  id = c(id1,id2,id3,id4)
  if(length(id) > 0)
    speciesL = speciesL[-id]
  df = read.csv("D:/project/lipidomics/data/lipidHomeDataExport.csv", header = F)
  l = unlist(strsplit(speciesL, split ="-"))
  if(tumour){
    l[which(l=="TRG")] = "TG"
    l[which(l=="DRG")] = "DG"
    l[which(l=="MRG")] = "MG"  
    l[which(l=="DHCer")] = "dhCer"
  }
  speciesL = paste(l[seq(2,length(l),2)],l[seq(1,length(l),2)])
  df = df[which(df$V2 %in% speciesL),]
  if(nrow(df)>0){
    for(i in 1:nrow(df)){
      name = unlist(strsplit(as.character(df[i,'V2']), split =" "))
      if(tumour){
        if(name[1]=="TG")
          name[1]="TRG"
        if(name[1]=="DG")
          name[1]="DRG"
        if(name[1]=="MG")
          name[1]="MRG"
        if(name[1]=="dhCer")
          name[1]="DHCer"
      }
      newN = paste(name[2],name[1],sep="-")
      mass[newN] = as.numeric(as.character(df[i,'V9']))
    }
  }
  mass
}
# Convert mass concentration to molar concentration
getMolarData <- function(data, tumour){
  l = colnames(data)
  mass = getMassSpecies(l, tumour)
  data = t(t(data) / mass)
  data
}
getPatientFattyAcidMatrix <- function(tumor){
  patientData = getPatientFattyAcidData()
  patients = rownames(patientData)
  specieData = getFattyAcidData()
  species = rownames(specieData)[which(!grepl('17_0-FA',specieData[,1]))]
  patientFattyAcidMatrix = matrix(ncol = length(patients), nrow =  length(species))
  rownames(patientFattyAcidMatrix) = species
  #if(tumor){
  #  ts = "t"
  #}else{
  #  ts = "n"
  #}
  #colnames(patientFattyAcidMatrix) = paste(ts,patients, sep = "")
  colnames(patientFattyAcidMatrix) = patients
  for(i in 1:length(patients)){
    quantityDataQ= paste("select h.quantity, CONCAT(s.carbons, '_', s.double_bonds, '-FA') as name from patient_has_fatty_acid_species as h, patient as p, fatty_acid_species as s where h.f_patient_id = p.patient_id and h.tumour = ",tumor," and h.f_fatty_acid_species_id = s.fatty_acid_species_id and p.number  = ",patients[i]," and h.quantity != 400 order by s.carbons", sep = "")    
    quantitydata = dbGetQuery(con, quantityDataQ)
    patientFattyAcidMatrix[,i] = quantitydata[,1]	
    #print(nrow(quantitydata))	
  }
  
  t(patientFattyAcidMatrix)
}
# this function is a new version of filterClassPSM function, which deals with all lipid including PIPx
filterClassPSM_1 <- function(data, lipidClass){
  
  lipidSpecies = vector()
  tmp = colnames(data)
  if(toupper(lipidClass) == "S1P" | toupper(lipidClass) == "SPC" | toupper(lipidClass) == "CH" | toupper(lipidClass) =="SG"){
    data1 = matrix(data[,which(grepl(lipidClass, tmp))])
    colnames(data1) = lipidClass
    rownames(data1) = rownames(data)
    data = data1
  }else if(toupper(lipidClass) == "PIP" | toupper(lipidClass) == "PIP2" | toupper(lipidClass) == "PIP3"){
    k = 1
    for(i in 1:length(tmp)){
      # only take species which are formatted as carbon:bond
      if(grepl(":",tmp[i])){
        lipidSpecies[k] = tmp[i]
        k = k + 1
      }
    }
    tSpeciesVec = unlist(strsplit(lipidSpecies, split ="-"))
    l = length(tSpeciesVec)
    data = data[,lipidSpecies]  
    data = data[,which(toupper(tSpeciesVec[seq(2,l,2)]) == toupper(lipidClass))]
  }else{
    lpClass = unlist(strsplit(tmp,"-"))
    data = data[, which(toupper(lpClass[seq(2,length(lpClass),by=2)]) == toupper(lipidClass))]
  }
  data  
}

# this is for cell data
convertToMatrix <- function(data){
  row = nrow(data[[1]])
  speciesClassNum = length(data)
  colNames = vector()
  quantity = vector()
  k = 1
  for(i in 1:speciesClassNum){
    colN = colnames(data[[i]])
    for(j in 1:length(colN)){
      colNames[k] = colN[j]
      quantity[k] = data[[i]][1,j]
      k = k + 1
    }
  }
  # in case of one row data only
  m = matrix(nrow = 1, ncol = length(colNames))
  m[1,] = quantity
  colnames(m) = colNames
  m
}
# get data for PCA analysis
getPatientDataMatrix <-function(tumor){
  patientData = getAllPatientData()
  patients = rownames(patientData)
  patients = rownames(patientData)
  specieData = getSpecieData()
  species = rownames(specieData)[which(specieData[,2] == 0)]
  #lipidClass = getLipidClass()
  #lipidName = rownames(lipidClass)
  patientSpecieMatrix = matrix(ncol = length(patients), nrow =  length(species)+3)
  if(tumor){
    ts = "t"
  }else{
    ts = "n"
  }
  colnames(patientSpecieMatrix) = paste(ts,patients, sep = "")
  #colnames(patientSpecieMatrix) = patients
  
  for(i in 1:length(patients)){
    quantityDataQ= paste("select h.quantity, s.name, c.name from patient_has_lipid_species as h, patient as p, lipid_species as s, lipid_class as c where h.l_patient_id = p.patient_id and h.tumour = ",tumor," and h.l_lipid_species_id = s.lipid_species_id and p.number  = ",patients[i]," and s.standard = false and s.l_lipid_class_id = c.lipid_class_id order by c.name", sep = "")  	
    quantitydata = dbGetQuery(con, quantityDataQ)
    # measure lipid concentrations on log10 scale
    age = as.numeric(patientData[i,2])
    if(age >= 31 && age <= 40){
      ageStr = "31-40"
    }else if(age >= 41 && age <= 49){
      ageStr = "41-49"
    }else if(age >= 50 && age <= 59){
      ageStr = "50-59"
    }else if(age >= 60 && age <= 72){
      ageStr = "60-72"
    }else if(age >= 73 && age <= 98){
      ageStr = "73-98"
    }
    
    dataVec = c(quantitydata[,1],ageStr,patientData[i,1],patientData[i,3])
    patientSpecieMatrix[,i] = dataVec	
    #print(nrow(quantitydata))	
  }
  rownames(patientSpecieMatrix) = c(species,'age','genre','stage')
  t(patientSpecieMatrix)
}
plotPCAAnalysis <- function(){
  pr = getPatientDataMatrix(TRUE)
  pr1 = data.frame(pr)
  # get quantitative varibles
  x1 = pr[,1:547]
  
  # convert to numeric
  class(x1) = "numeric"
  # get categorical variables (age, genre, stage)
  x2 = pr1[,548:550]
  scale(x2, center = TRUE, scale = TRUE)
  # do PCA for mixed data types
  library(PCAmixdata)
  obj = PCAmix(X.quanti = x1, X.quali = x2, ndim = 3)
  # factor scores for quantitative variables
  A1 = obj$quanti.cor
  # factor scores for categorical variables
  A2 = obj$categ.coord
  # factor scores for rows in the data
  head(F)
  # Component map with factor scores of the data (rows)
  setwd("D:/project/lipidomics/data/lipid analysis")
  #par(mfrow=c(2,2))
  pdf(file = "pca_component_factor_score.pdf")
  plot(obj, choice = "ind", ces = 0.6, habillage = "ind")
  dev.off()
  # Component map with factor scores of the numerical columns
  pdf(file = "pca_numerical_col_factor_score.pdf")
  plot(obj, choice = "cor", ces = 0.6, habillage = "cor")
  dev.off()
  pdf(file = "pca_categorical_col_factor_score.pdf")
  # Component map with factor scores of the categorical columns
  plot(obj, choice = "levels", ces = 0.6, habillage = "levels")
  dev.off()
  # contributions of the variables
  pdf(file = "pca_variables_contribution.pdf")
  plot(obj, choice = "sqload", ces = 0.6, habillage = "sqload")
  dev.off()
  # clustering
  # Construction of the hierarchy
  library(ClustOfVar)
  pdf(file = "pca_hierarchical_cluster.pdf")
  tree = hclustvar(X.quanti = x1, X.quali = x2)
  # Graphical representation
  plot(tree, cex = 0.5)
  # Partition in 6 clusters
  part = cutreevar(tree, 6)
  dev.off()
  plotData = obj$levels$coord
  # get variables names
  vars = rownames(plotData)
  # extract age group
  ageVars = vars[which(!is.na(as.numeric(vars)))]
  ageCoord = plotData[ageVars,1:2]
  genreVars = vars[which(vars == "female" | vars == "male")] 
  genreCoord = plotData[genreVars,1:2]
  for(i in 1:length(genreVars)){
    if(genreVars[i] == "female")
    {
      genreVars[i] = "F"
    }
    else{
      genreVars[i] = "M"
    }
  }
  
  stageVars = vars[which(vars %in% c("adenoma","dukes a","dukes b","dukes c","dukes d") )]
  stageCoord = plotData[stageVars,1:2]
  for(i in 1:length(stageVars)){
    if(stageVars[i] == "adenoma"){
      stageVars[i] = "Ade"
    }
    else if(stageVars[i] == "dukes a"){
      stageVars[i] = "A"
    }else if(stageVars[i] == "dukes b"){
      stageVars[i] = "B"
    }
    else if(stageVars[i] == "dukes c"){
      stageVars[i] = "C"
    }else {
      stageVars[i] = "D"
    }
    
  }
  # Do PCA with using non-numeric variables as suppplementary variables
  
  
  # extract numeric columns
  # take first 65 rows for PCA analysis
  # install.packages("devtools")
  devtools::install_github("kassambara/factoextra")
  
  # load
  library("factoextra")
  my.active = matrix(0, nrow = nrow(pr1[1:70,1:547]), ncol = ncol(pr1[1:70,1:547])) 
  my.active = pr1[1:70,1:547]
  matrix = as.numeric(as.matrix(my.active))
  dim(matrix) = dim(my.active)
  rownames(matrix) = rownames(my.active)
  res.pca = prcomp(matrix, center = TRUE, scale = TRUE)
  #the matrix of variable loadings (columns are eigenvectors)
  head(unclass(res.pca$rotation)[, 1:4])
  # Eigenvalues
  eig = (res.pca$sdev)^2
  
  # Variances in percentage
  variance = eig*100/sum(eig)
  
  # Cumulative variances
  cumvar = cumsum(variance)
  
  eig.matrix = data.frame(eig = eig, variance = variance,
                          cumvariance = cumvar)
  head(eig.matrix)
  eig.val = get_eigenvalue(res.pca)
  head(eig.val)
  barplot(eig.matrix[, 2], names.arg=1:nrow(eig.matrix), 
          main = "Variances",
          xlab = "Principal Components",
          ylab = "Percentage of variances",
          col ="steelblue")  
  # Add connected line segments to the plot
  lines(x = 1:nrow(eig.matrix), 
        eig.matrix[, 2], 
        type="b", pch=19, col = "red")
  # plot variance
  fviz_screeplot(res.pca, ncp=10)
  #plot eigenvalues
  fviz_screeplot(res.pca, ncp=10, choice="eigenvalue")
  # Variable correlation/coordinates
  # this is to see how each component is strongly correlated with variables
  loadings = res.pca$rotation
  sdev = res.pca$sdev
  var <- get_pca_var(res.pca)
  var.coord = var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
  head(var.coord[, 1:4])
  
  # Plot the correlation circle
  a <- seq(0, 2*pi, length = 100)
  plot( cos(a), sin(a), type = 'l', col="gray",
        xlab = "PC1 - 13.6%",  ylab = "PC2 - 9.1%")
  
  abline(h = 0, v = 0, lty = 2)
  
  # Add active variables
  arrows(0, 0, var.coord[, 1], var.coord[, 2], 
         length = 0.1, angle = 15, code = 2)
  
  # Add labels
  text(var.coord, labels=colnames(my.active), cex = 1, adj=1)
  
  # add supplmentary variables
  # plot for age
  pdf(file = "pca_suppl_age1.pdf")
  quali.sup <- as.factor(as.matrix(pr1[1:70, 548]))
  
  fviz_pca_ind(res.pca, data = matrix, axes = c(1,2),
               habillage = quali.sup, addEllipses = TRUE, ellipse.level = 0.68) +
    theme_minimal()
  
  dev.off()
  # plot for genre
  pdf(file = "pca_suppl_genre1.pdf")
  quali.sup <- as.factor(as.matrix(pr1[1:70, 549]))
  fviz_pca_ind(res.pca, data = matrix, axes = c(1,2),
               habillage = quali.sup, addEllipses = TRUE, ellipse.level = 0.68) +
    theme_minimal()
  dev.off()
  # plot for stage
  pdf(file = "pca_suppl_stage1.pdf")
  quali.sup <- as.factor(as.matrix(pr1[1:70, 550]))
  fviz_pca_ind(res.pca, data = matrix, axes = c(1,3),
               habillage = quali.sup, addEllipses = TRUE, ellipse.level = 0.68) +
    theme_minimal()
  
  # plot individual factor map
  ind.coord <- res.pca$x
  plot(ind.coord[,1], ind.coord[,2], pch = 19,  
       xlab="PC1 - 13.6%",ylab="PC2 - 9.1%")
  abline(h=0, v=0, lty = 2)
  text(ind.coord[,1], ind.coord[,2], labels=rownames(ind.coord),
       cex=0.7, pos = 3)
  # plot individual factor map using sum of squared cosines
  # The sum of squared cosines over all axes is 1
  # the closer the observation to 1, the better interpreabilty of representation
  fviz_pca_ind(res.pca, data = matrix, col.ind="cos2") +
    scale_color_gradient2(low="white", mid="blue", 
                          high="red", midpoint=0.50) + theme_minimal()
  ####################################################
  par(new=TRUE)
  
  plot(ageCoord,
       , xlab = "Dim 1 (9.061 %)"
       , ylab = "Dim 2 (6.726 %)"
       , type = "p"
       , cex = 0.8
       , pch = 17
       , col = "red"
  )
  text(ageCoord, labels=ageVars, cex= 0.7, pos=3)
  par(new=TRUE)
  plot(genreCoord,
       , xlab = ""
       , ylab = ""
       , type = "p"
       , cex = 0.8
       , pch = 19
       , col = "blue"
  )
  text(genreCoord, labels=genreVars, cex= 0.7, pos=3)
  par(new=TRUE)
  plot(stageCoord,
       , xlab = ""
       , ylab = ""
       , axes = FALSE
       , type = "p"
       , cex = 0.8
       , pch = 15
       , col = "green"
  )
  text(stageCoord, labels=stageVars, cex= 0.7, pos=3)
}
# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
getPatientSpecieMatrix <-function(tumor){
  patientData = getPatientData()
  patients = rownames(patientData)
  specieData = getSpecieData()
  species = rownames(specieData)[which(specieData[,2] == 0)]
  patientSpecieMatrix = matrix(ncol = length(patients), nrow =  length(species))
  rownames(patientSpecieMatrix) = species
  #if(tumor){
  #  ts = "t"
  #}else{
  #  ts = "n"
  #}
  #colnames(patientSpecieMatrix) = paste(ts,patients, sep = "")
  colnames(patientSpecieMatrix) = patients
  for(i in 1:length(patients)){
    quantityDataQ= paste("select h.quantity, s.name, c.name from patient_has_lipid_species as h, patient as p, lipid_species as s, lipid_class as c where h.l_patient_id = p.patient_id and h.tumour = ",tumor," and h.l_lipid_species_id = s.lipid_species_id and p.number  = ",patients[i]," and s.standard = false and s.l_lipid_class_id = c.lipid_class_id order by 2", sep = "")		
    quantitydata = dbGetQuery(con, quantityDataQ)
    patientSpecieMatrix[,i] = quantitydata[,1]	
    #print(nrow(quantitydata))	
  }
  
  t(patientSpecieMatrix)
}


normalisePatientSpecieMatrix <- function(patientSpecieMatrix){
  for(i in 1:nrow(patientSpecieMatrix)){
    patientSpecieMatrix[i,] = patientSpecieMatrix[i,]/sum(patientSpecieMatrix[i,])
  }
  patientSpecieMatrix
}
normaliseSpecieMatrix <- function(mt){
  for(i in 1:nrow(mt)){
    mt[i,] = mt[i,]/sum(mt[i,])
  }
  mt
}
filterClassPSM <- function(PSM, lclass){
  specieData = getSpecieData()
  cespecies = rownames(specieData)[which(specieData[,1] == lclass)]
  coi = which(colnames(PSM) %in% cespecies)	
  PSM[,coi,drop=FALSE]	
}
getSpecieData <- function(){
  specieInformationQ = "select s.name,c.name,s.standard,s.carbons, s.double_bonds from lipid_species as s, lipid_class as c where s.l_lipid_class_id = c.lipid_class_id and s.standard=false order by 1"
  specieData = dbGetQuery(con, specieInformationQ)	
  rownames(specieData) = specieData[,1]
  specieData = specieData[,-1]
  specieData
}
getFattyAcidData <- function(){
  fattyAcidInformationQ = "select CONCAT(s.carbons, '_', s.double_bonds, '-FA') as name from fatty_acid_species as s order by 1"
  specieData = dbGetQuery(con, fattyAcidInformationQ)	
  rownames(specieData) = specieData[,1]
  specieData
}
applyDetectionLimit<- function(patientSpecieMatrix){
  specieData = getSpecieData()
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  for(i in 1:ncol(patientSpecieMatrix)){
    specie = colnames(patientSpecieMatrix)[i]
    class = specieData[which(rownames(specieData) == specie),1]
    detectionLimit = detectionLimits[which(detectionLimits[,1] == class),2]
    patientSpecieMatrix[which(patientSpecieMatrix[,i] < detectionLimit),i] = detectionLimit
  }
  patientSpecieMatrix
}



getRatioPSM<- function(tPSM, nPSM){
  applyDetectionLimit(tPSM)/applyDetectionLimit(nPSM)
}

stageMeta <- function(patientData){
  meta = 1
  patientMeta = patientData[,meta]
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = vector()
  result
}

metastasisMeta <- function(patientData){
  meta = 13
  patientMeta = patientData[,meta]
  result = list()
  a= is.na(patientMeta)
  for(i in 1:length(a)){
    if(a[i]){
      if(patientData[i,1] == "dukes d"){
        patientMeta[i] = "M"
      }else{
        patientMeta[i] = "N"
      }
      
    }else{
      if(patientData[i,1] == "dukes d"){
        patientMeta[i] = "M"
      }else{
        patientMeta[i] = "N"
      }	
    }
  }
  result$patientMeta = patientMeta
  result$toRemove = vector()
  result
}

maleMeta <- function(patientData){
  meta = 6
  patientMeta = patientData[,meta]	
  a = which(is.na(patientMeta))
  toRemove = a
  patientMeta = as.character(patientMeta[-toRemove])
  for(i in 1:length(patientMeta)){
    if(patientMeta[i] == "0"){
      patientMeta[i] = "F"
    }else{
      patientMeta[i] = "M"
    }
  }
  
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = toRemove
  result	
}

mortalityMeta <- function(patientData){
  meta = 15
  patientMeta = patientData[,meta]	
  a = which(is.na(patientMeta))	
  patientMeta[a] = "alive"
  patientMeta[-a] = "dead"		
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = vector()
  result	
}

siteMeta <- function(patientData){
  meta = 7
  patientMeta = patientData[,meta]
  a = names(which(table(patientMeta)>5))
  toKeep = which(patientMeta %in% a)
  v = 1:length(patientMeta)
  toRemove = v[-toKeep]
  patientMeta = patientMeta[-toRemove]
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = toRemove
  result
}

#this is probably regression even though we have distinct classes
sizeMeta<- function(patientData){
  meta = 11
  patientMeta = patientData[,meta]	
  a = which(is.na(patientMeta))
  toRemove = a
  patientMeta = patientMeta[-toRemove]
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = toRemove
  result
}

spreadMeta = function(patientData){
  meta = 12
  patientMeta = patientData[,meta]	
  a = which(is.na(patientMeta))
  toRemove = a
  patientMeta = patientMeta[-toRemove]
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = toRemove
  result	
}
sampleMeta<- function(pNamesVec){
  patientMeta = vector(length = length(pNamesVec))
  for(i in 1:length(pNamesVec)){
    patientMeta[i] = substr(pNamesVec[i],1,1)
  }
  toRemove = vector()
  result = list()
  result$patientMeta = patientMeta
  result$toRemove = toRemove
  result	
}

tExpression <- function(tData, nData){
  tVec = vector(length = ncol(tData))
  for(i in 1:ncol(tData)){
    t = t.test(tData[,i], nData[,i], paired = T)
    tVec[i] = t$p.value
  }
  #bontVec = tVec * i
  #bontVec
  tVec
}

mdsPlot <- function(PSMO, nFeatures){
  d <- dist(PSMO[,1:nFeatures]) # euclidean distances between the rows
  mds <- isoMDS(d, k=2) # k is the number of dim
  mds # view results
  # plot solution
  x <- mds$points[,1]
  y <- mds$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",  main=paste("isoMDS of ", nFeatures," most important sample species",sep = ""), type="n")
  text(x, y, labels = row.names(PSMO), cex=.7) 
}

pcaPlot <- function(PSMO, nFeatures){
  pcaa <- princomp(PSMO[,1:nFeatures], cor=TRUE)
  #plot(pcaa,type="lines") # scree plot
  biplot(pcaa, main = paste("isoMDS of ", nFeatures," most important sample species",sep = "")) 
}

plotFeatureRatioDistribution<- function(logData, features){
  ldata = logData[,features]		
  for(i in 1:ncol(ldata)){
    feature = features[i]
    pdf(width = 8 , height = 8 , file = paste("featureRatio_",feature,".pdf", sep = ""))		
    plot(density(logData[,i]))
    dev.off()
  }
  
  
}

getPatientColVec <- function(patientNumbers, meta){
  a = list()
  colVec = vector(length = length(patientNumbers))
  pData = getPatientData()
  mClasses = unique(pData[,meta])
  
  colPal = rainbow(length(mClasses))
  
  
  for(i in 1:length(patientNumbers)){
    p = substr(patientNumbers[i], 2, nchar(patientNumbers[i]))
    roi = which(rownames(pData) == p)
    moi = which(mClasses == pData[roi,meta])	
    if(length(moi)==0){
      colVec[i] = "black"
      colPal[1] = "black"
      if(length(which(mClasses=="unknown"))==0){
        mClasses[which(is.na(mClasses))] = "unknown"
      }
    }else{
      colVec[i] = colPal[moi]
    }
    
  }
  a$classcols = colPal
  a$colVec = colVec
  a$classes = mClasses
  a
}

removeSparseSpecies<-function(PSM,cutOff,threshold){
  excluders = vector()
  for(i in 1:ncol(PSM)){
    toExclude = which(PSM[,i] <= threshold)		
    if(length(toExclude)>cutOff){
      excluders = c(excluders, i)
    }
  }
  PSM[,-excluders]
}
# get set of 
#you are here
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

getFAs <- function(species){
  fas = vector(length = length(species))
  for(i in 1:length(species)){
    fas[i] = as.vector(unlist(strsplit(species[i], split ="-")))[length(as.vector(unlist(strsplit(species[i], split ="-"))))-1]
  }
  fas
}

closeDBConnection <- function(){
  all_cons <- dbListConnections(MySQL())
  for(con in all_cons)
    dbDisconnect(con)
}
getFeatureSumCoefV <- function(PSMN,PSMT){
  ncoefv = vector(length = ncol(PSMN))
  tcoefv = vector(length = ncol(PSMT))
  nans = vector()
  for(i in 1:ncol(PSMN)){
    ncoefv[i] = calculateCoefv(PSMN[,i])	
    tcoefv[i] = calculateCoefv(PSMT[,i])
  }
  
  coefvM = matrix(nrow =length(ncoefv), ncol = 3)
  rownames(coefvM) = colnames(PSMN)
  coefvM[,1] = ncoefv
  coefvM[,2] = tcoefv
  for(i in 1:nrow(coefvM)){
    coefvM[i,3] = coefvM[i,1] + coefvM[i,2]	
  }
  coefvM
}


calculateCoefv <- function(aVec){
  asdev = sd(aVec)
  amean = mean(aVec)
  result = (asdev*100)/amean
  result	
}
# This function is to calculate z-score of subpathway when ignoring the layers
getSubPathwayZScore <- function(pathway, tData, nData,alt){  
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getEdgeZScore(pathway[i],pathway[i+1], tData, nData,alt)
  }
  z_score = sum(z)/sqrt(size)
  z_score
}
getFattyAcidEdgeZScore <- function(pathway, refData, controlData, alt){
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getFattyAcidZScore(pathway[i],pathway[i+1], refData, controlData, alt)
  }
  z_score = sum(z)/sqrt(size)
  z_score
  
}
# calculate z-sore of subpathway for Hbec
getHbecSubPathwayZScore <- function(pathway, refData, controlData,alt){  
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getHbecEdgeZScore(pathway[i],pathway[i+1], refData, controlData,alt)
  }
  z_score = sum(z)/sqrt(size)
  z_score
}
getHbecSubPathwayZScore0 <- function(pathway, refData, controlData,alt,type){  
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getHbecEdgeZScore0(pathway[i],pathway[i+1], refData, controlData,alt,type)
  }
  z_score = sum(z)/sqrt(size)
  z_score
}
getHbecSubPathwayZScore1 <- function(pathway, refData, controlData,alt){  
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getHbecEdgeZScore1(pathway[i],pathway[i+1], refData, controlData,alt)
  }
  z_score = sum(z)/sqrt(size)
  z_score
}
getFattyAcidSubPathwayZScore <- function(pathway, refData, controlData,alt){  
  size = length(pathway)-1
  z = vector()
  z_score = 0
  for(i in 1:size){
    # get z-score for each edge
    z[i] = getFattyAcidEdgeZScore(pathway[i],pathway[i+1], refData, controlData,alt)
  }
  z_score = sum(z)/sqrt(size)
  z_score
}
getSubPathwayZScoreByStage <- function(pathway, tData, nData,alt,sigLevel){
  # get data by stages
  #patientData = getPatientData()
  #pml = stageMeta(patientData)
  #patientMeta = pml$patientMeta
  #split into stages
  #stages = unique(patientMeta)
  stages = getPatientStageData()
  # only consider stages which have sufficiently large amount of data
  # Set the threshold to be 20 (ideally at least 30)
  index = which(stages[,2] >= 20)
  stages = stages[index,1]
  stagePatients =  list()
  size = length(pathway)-1
  z_score = 0
  s = 0
  z = vector()
  p = vector()
  n = length(stages)
  for(i in 1:n){
    poi = which(patientMeta == stages[i])  
    stagePatients = tData[poi,]
    nnData = nData[poi,]
    for(j in 1:size){
      # get z-score for each edge
      z[j] = getEdgeZScore(pathway[j],pathway[j+1], stagePatients, nnData,alt)
    }
    # convert to p-value
    p[i] = 1 - pnorm(sum(z)/sqrt(size))
  }
  p = sort(p, decreasing = FALSE)
  pp = vector()
  zz = vector()
  for(i in 1:n){
    s = 0
    for(j in i:n){
      s = s + (factorial(n) / (factorial(j) * factorial(n-j))) * p[i]^j*(1-p[i])^(n-j)
    }
    pp[i] = s
    zz[i] = qnorm(1 - pp[i])
  }
  z_score = max(zz)
  z_score
}
# this function is to plot weighted pathway graph
plotPathwayGraph <- function(reactant, product, ZScore, pdfFileName){
  edge = data.frame(from = reactant, to = product, thickness = ZScore)
  qgraph(as.matrix(edge),edge.labels=T,label.prop = 0.6, mode="direct", edge.label.cex=0.6, asize=1.5, vsize=4,fade=F,filename=pdfFileName,filetype = "pdf", height = 5, width = 10)   
  #edge = data.frame(from = reactant, to = product)
  #n = nrow(edge)
  #edge = cbind(edge, seq(min(2*ZScore),max(2*ZScore),length=n))
  #qgraph(as.matrix(edge),edge.labels=F,label.prop = 0.6, mode="direct", edge.label.cex=0.6, asize=1.5, vsize=4,fade=FALSE, filename=pdfFileName,filetype = "pdf", height = 5, width = 10)   
  
}
getEdgeZScoreByStage <- function(reactant, product, tData, nData){
  # get data by stages
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = unique(patientMeta)
  stagePatients =  list()
  z_score = 0
  s = 0
  p = vector()
  n = length(stages)
  for(i in 1:n){
    poi = which(patientMeta == stages[i])  
    stagePatients = tData[poi,]
    trData = filterClassPSM(stagePatients,reactant)
    tpData = filterClassPSM(stagePatients,product)
    nrData = filterClassPSM(nData,reactant)
    npData = filterClassPSM(nData,product)
    if(is.vector(nrData)){
      nnrData = nrData[poi]
    }else{
      nnrData = nrData[poi,]  
    }
    if(is.vector(npData)){
      nnpData = npData[poi]
    }else{
      nnpData = npData[poi,]  
    }
    tSumV = list()
    nSumV = list()
    p[i] = 0
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
    if(is.vector(nnpData)){
      nSumV[[1]] = nnpData
    }else{
      nSumV[[1]] = rowSums(nnpData)  
    }
    if(is.vector(nnrData)){
      nSumV[[2]] = nnrData
    }else{
      nSumV[[2]] = rowSums(nnrData)  
    }
    tDataTest = tSumV[[1]]/tSumV[[2]]
    nDataTest = nSumV[[1]]/nSumV[[2]]
    
    t = t.test(tDataTest,nDataTest, alternative= "greater", paired = T)
    if(is.finite(t$p.value)){
      p[i] = t$p.value
    }  
  }
  p = sort(p, decreasing = FALSE)
  pp = vector()
  zz = vector()
  for(i in 1:n){
    s = 0
    for(j in i:n){
      s = s + (factorial(n) / (factorial(j) * factorial(n-j))) * p[i]^j*(1-p[i])^(n-j)
    }
    pp[i] = s
    zz[i] = qnorm(1 - pp[i])
  }
  z_score = max(zz)
  z_score
}
getFattyAcidEdgeZScoreByStage <- function(reactant, product, tData, nData,alt){
  # get data by stages
  patientData = getPatientFattyAcidData()
  stages = unique(patientData[,2])
  stagePatients =  list()
  z_score = 0
  s = 0
  p = vector()
  n = length(stages)
  for(i in 1:n){
    poi = which(patientData[,2] == stages[i])  
    stagePatients = tData[poi,]
    trData = tData[poi, which(colnames(stagePatients) == reactant)]
    tpData = tData[poi, which(colnames(stagePatients) == product)]
    nrData = nData[poi, which(colnames(stagePatients) == reactant)]
    npData = nData[poi, which(colnames(stagePatients) == product)]
    tSumV = list()
    nSumV = list()
    p[i] = 0
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
    tDataTest = tSumV[[1]]/tSumV[[2]]
    nDataTest = nSumV[[1]]/nSumV[[2]]
    tDataTest[is.infinite(tDataTest)] = 1e6
    nDataTest[is.infinite(nDataTest)] = 1e6
    t = t.test(tDataTest,nDataTest, alternative= alt, paired = T)
    if(is.finite(t$p.value)){
      p[i] = t$p.value
    }  
  }
  p = sort(p, decreasing = FALSE)
  pp = vector()
  zz = vector()
  for(i in 1:n){
    s = 0
    for(j in i:n){
      s = s + (factorial(n) / (factorial(j) * factorial(n-j))) * p[i]^j*(1-p[i])^(n-j)
    }
    pp[i] = s
    zz[i] = qnorm(1 - pp[i])
  }
  z_score = max(zz)
  z_score
}
# This function is to find all active pathways by performing a deep search
getFattyAcidSigPathway <- function(dataL,pathway, edgeW, root, curData, preData,alt){
  
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        zscore = getFattyAcidSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        if(zscore > 1.645){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getFattyAcidSigPathway(dataL, apathway, edgeW, nodeL[i], curData, preData,alt)
        }
      }
    }
    
  }
  dataL
  
}
getSigPathway <- function(dataL,pathway, edgeW, root, curData, preData,alt,isTumour){
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        if(isTumour){
          zscore = getSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        }else{
          zscore = getHbecSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)  
        }
        
        if(zscore > 1.645){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getSigPathway(dataL, apathway, edgeW, nodeL[i], curData, preData,alt,isTumour)
        }
      }
    }
    
  }
  dataL
}
getSigPathway0 <- function(dataL,pathway, edgeW, root, curData, preData,alt,isTumour,type){
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        if(isTumour){
          zscore = getSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        }else{
          zscore = getHbecSubPathwayZScore0(c(pathway,nodeL[i]), curData, preData,alt,type)  
        }
        
        if(zscore > 1.645){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getSigPathway0(dataL, apathway, edgeW, nodeL[i], curData, preData,alt,isTumour,type)
        }
      }
    }
    
  }
  dataL
}
getSigPathway_1 <- function(dataL,pathway, edgeW, root, curData, preData,alt,sigLevel,isTumour){
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  z = qnorm(1-sigLevel)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        if(isTumour){
          zscore = getSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        }else{
          zscore = getHbecSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)  
        }
        
        if(zscore > z){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getSigPathway(dataL, apathway, edgeW, nodeL[i], curData, preData,alt,sigLevel,isTumour)
        }
      }
    }
    
  }
  dataL
}
getSigPathway1 <- function(dataL,pathway, edgeW, root, curData, preData,alt,isTumour){
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        if(isTumour){
          zscore = getSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        }else{
          zscore = getHbecSubPathwayZScore1(c(pathway,nodeL[i]), curData, preData,alt)  
        }
        
        if(zscore > 1.645){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getSigPathway1(dataL, apathway, edgeW, nodeL[i], curData, preData,alt,isTumour)
        }
      }
    }
    
  }
  dataL
}
getSigPathway1_1 <- function(dataL,pathway, edgeW, root, curData, preData,alt,sigLevel,isTumour){
  visited = dataL[[2]]
  visited[root] = TRUE
  dataL[[2]] = visited
  edgeWL = edgeW[[root]]
  # get list of its neighbors with decreasing order of ZScores
  neighborL = sort(edgeWL, decreasing = TRUE)
  # get the next node
  nodeL = names(neighborL)
  z = qnorm(1-sigLevel)
  if(length(nodeL)>0){
    for(i in 1:length(nodeL)){
      if(!visited[nodeL[i]]){
        if(isTumour){
          zscore = getSubPathwayZScore(c(pathway,nodeL[i]), curData, preData,alt)
        }else{
          zscore = getHbecSubPathwayZScore1(c(pathway,nodeL[i]), curData, preData,alt)  
        }
        
        if(zscore > z){
          pathVec = dataL[[1]]
          path = paste(pathway, collapse = "->")
          apathway = c(pathway,nodeL[i])
          newPath = paste(apathway, collapse = "->")
          if(length(pathVec)>0){
            pathN = names(pathVec)
            if(path %in% pathN){
              pathVec = pathVec[-which(pathN==path)]
            }
          }
          pathVec[newPath] = zscore
          dataL[[1]] = pathVec
          dataL = getSigPathway1(dataL, apathway, edgeW, nodeL[i], curData, preData,alt,sigLvel,isTumour)
        }
      }
    }
    
  }
  dataL
}
getFattyAcidZScore <- function(reactant, product, refData, controlData, alt){
  tmp = colnames(refData)
  rRefData = refData[,which(grepl(reactant, tmp))]
  pRefData = refData[,which(grepl(product, tmp))]
  rControlData = controlData[,which(grepl(reactant, tmp))]
  pControlData = controlData[,which(grepl(product, tmp))]
  refSumV = list()
  controlSumV = list()
  refSumV[[1]] = pRefData
  controlSumV[[1]] = pControlData
  refSumV[[2]] = rRefData
  controlSumV[[2]] = rControlData
  z_score = 0
  ind1 = which(refSumV[[2]] == 0)
  ind2 = which(controlSumV[[2]] == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(refSumV[[2]])-1){
    if(l>0){
      tDataTest = refSumV[[1]][-ind]/refSumV[[2]][-ind]  
      nDataTest = controlSumV[[1]][-ind]/controlSumV[[2]][-ind]  
    }else{
      tDataTest = refSumV[[1]]/refSumV[[2]]  
      nDataTest = controlSumV[[1]]/controlSumV[[2]]
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
  
}
getSMSEnzymeScore <- function(refData, controlData,alt){
  
  rRefData1 = filterClassPSM_1(refData,"PC")
  rRefData2 = filterClassPSM_1(refData,"Cer")
  pRefData1 = filterClassPSM_1(refData,"SM")
  pRefData2 = filterClassPSM_1(refData,"DG")
  
  rControlData1 = filterClassPSM_1(controlData,"PC")
  rControlData2 = filterClassPSM_1(controlData,"Cer")
  pControlData1 = filterClassPSM_1(controlData,"SM")
  pControlData2 = filterClassPSM_1(controlData,"DG")
  
  rRefSum1 = rowSums(rRefData1)
  rRefSum2 = rowSums(rRefData2)
  
  pRefSum1 = rowSums(pRefData1)
  pRefSum2 = rowSums(pRefData2)
  
  rControlSum1 = rowSums(rControlData1)
  rControlSum2 = rowSums(rControlData2)
  
  pControlSum1 = rowSums(pControlData1)
  pControlSum2 = rowSums(pControlData2)
  
  pRef1 = pRefSum1*pRefSum2
  pControl1 = pControlSum1*pControlSum2
  
  rRef1 = rRefSum1*rRefSum2
  rControl1 = rControlSum1*rControlSum2
    
  z_score = 0
  ind1 = which(pRef1 == 0)
  ind2 = which(pControl1 == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(pRef1)-1){
    if(l>0){
      tDataTest = pRef1[-ind]/rRef1[-ind]  
      nDataTest = pControl1[-ind]/rControl1[-ind]  
    }else{
      tDataTest = pRef1/rRef1
      nDataTest = pControl1/rControl1
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
  
}
getSMaseEnzymeScore <- function(refData, controlData,alt){
  
  rRefData1 = filterClassPSM_1(refData,"SM")
  rRefData2 = filterClassPSM_1(refData,"DG")
  pRefData1 = filterClassPSM_1(refData,"Cer")
  pRefData2 = filterClassPSM_1(refData,"PC")
  
  rControlData1 = filterClassPSM_1(controlData,"SM")
  rControlData2 = filterClassPSM_1(controlData,"DG")
  pControlData1 = filterClassPSM_1(controlData,"Cer")
  pControlData2 = filterClassPSM_1(controlData,"PC")
  
  rRefSum1 = rowSums(rRefData1)
  rRefSum2 = rowSums(rRefData2)
  
  pRefSum1 = rowSums(pRefData1)
  pRefSum2 = rowSums(pRefData2)
  
  rControlSum1 = rowSums(rControlData1)
  rControlSum2 = rowSums(rControlData2)
  
  pControlSum1 = rowSums(pControlData1)
  pControlSum2 = rowSums(pControlData2)
  
  pRef1 = pRefSum1*pRefSum2
  pControl1 = pControlSum1*pControlSum2
  
  rRef1 = rRefSum1*rRefSum2
  rControl1 = rControlSum1*rControlSum2
  
  z_score = 0
  ind1 = which(pRef1 == 0)
  ind2 = which(pControl1 == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(pRef1)-1){
    if(l>0){
      tDataTest = pRef1[-ind]/rRef1[-ind]  
      nDataTest = pControl1[-ind]/rControl1[-ind]  
    }else{
      tDataTest = pRef1/rRef1
      nDataTest = pControl1/rControl1
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
  
}
getHbecEdgeZScore0 <- function(reactant, product, refData, controlData, alt, type){
  rRefData = refData[,which(grepl(reactant,colnames(refData)))]
  pRefData = refData[,which(grepl(product,colnames(refData)))]
  rControlData = controlData[,which(grepl(reactant,colnames(controlData)))]
  pControlData = controlData[,which(grepl(product,colnames(controlData)))]
  refSumV = list()
  controlSumV = list()
  refSumV[[1]] = pRefData
  controlSumV[[1]] = pControlData
  refSumV[[2]] = rRefData
  controlSumV[[2]] = rControlData
  z_score = 0
  ind1 = which(refSumV[[2]] == 0)
  ind2 = which(controlSumV[[2]] == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(refSumV[[2]])-1){
    if(l>0){
      tDataTest = refSumV[[1]][-ind]/refSumV[[2]][-ind]  
      nDataTest = controlSumV[[1]][-ind]/controlSumV[[2]][-ind]  
    }else{
      tDataTest = refSumV[[1]]/refSumV[[2]]  
      nDataTest = controlSumV[[1]]/controlSumV[[2]]
    }
    if(type == 0){
      t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)  
    }else{
      t = t.test(tDataTest,nDataTest,alternative=alt)
    }
    
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
}
getHbecEdgeZScore <- function(reactant, product, refData, controlData, alt){
  rRefData = filterClassPSM_1(refData,reactant)
  pRefData = filterClassPSM_1(refData,product)
  rControlData = filterClassPSM_1(controlData,reactant)
  pControlData = filterClassPSM_1(controlData,product)
  refSumV = list()
  controlSumV = list()
  if(is.vector(pRefData)){
    refSumV[[1]] = pRefData
    controlSumV[[1]] = pControlData
  }else{
    refSumV[[1]] = rowSums(pRefData)  
    controlSumV[[1]] = rowSums(pControlData)
  }
  if(is.vector(rRefData)){
    refSumV[[2]] = rRefData
    controlSumV[[2]] = rControlData
  }else{
    refSumV[[2]] = rowSums(rRefData)  
    controlSumV[[2]] = rowSums(rControlData)  
  }
  z_score = 0
  ind1 = which(refSumV[[2]] == 0)
  ind2 = which(controlSumV[[2]] == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(refSumV[[2]])-1){
    if(l>0){
      tDataTest = refSumV[[1]][-ind]/refSumV[[2]][-ind]  
      nDataTest = controlSumV[[1]][-ind]/controlSumV[[2]][-ind]  
    }else{
      tDataTest = refSumV[[1]]/refSumV[[2]]  
      nDataTest = controlSumV[[1]]/controlSumV[[2]]
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
}
getHbecEdgeZScore1 <- function(reactant, product, refData, controlData, alt){
  rRefData = filterClassPSM_1(refData,reactant)
  pRefData = filterClassPSM_1(refData,product)
  rControlData = filterClassPSM_1(controlData,reactant)
  pControlData = filterClassPSM_1(controlData,product)
  refSumV = list()
  controlSumV = list()
  if(is.vector(pRefData)){
    refSumV[[1]] = pRefData
    controlSumV[[1]] = pControlData
  }else{
    refSumV[[1]] = rowSums(pRefData)  
    controlSumV[[1]] = rowSums(pControlData)
  }
  if(is.vector(rRefData)){
    refSumV[[2]] = rRefData
    controlSumV[[2]] = rControlData
  }else{
    refSumV[[2]] = rowSums(rRefData)  
    controlSumV[[2]] = rowSums(rControlData)  
  }
  z_score = 0
  ind1 = which(refSumV[[2]] == 0)
  ind2 = which(controlSumV[[2]] == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(refSumV[[2]])-1){
    if(l>0){
      tDataTest = refSumV[[1]][-ind]/refSumV[[2]][-ind]  
      nDataTest = controlSumV[[1]][-ind]/controlSumV[[2]][-ind]  
    }else{
      tDataTest = refSumV[[1]]/refSumV[[2]]  
      nDataTest = controlSumV[[1]]/controlSumV[[2]]
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt)
    z_score = qnorm(1 - t$p.value)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
  }    
  z_score
}
getFattyAcidEdgeZScore <- function(reactant, product, tData, nData,alt){
  tmp = colnames(tData)
  trData = tData[,which(grepl(reactant, tmp))]
  tpData = tData[,which(grepl(product, tmp))]
  nrData = nData[,which(grepl(reactant, tmp))]
  npData = nData[,which(grepl(product, tmp))]
 
  
  z_score = 0
  ind1 = which(trData == 0)
  ind2 = which(nrData == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(trData)-1){
    if(l>0){
      tDataTest = tpData[-ind]/trData[-ind]  
      nDataTest = npData[-ind]/nrData[-ind]  
    }else{
      tDataTest = tpData/trData
      nDataTest = npData/nrData
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
    
  }                       
  z_score
  
}
getEdgeZScore <- function(reactant, product, tData, nData,alt){
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
  z_score = 0
  ind1 = which(tSumV[[2]] == 0)
  ind2 = which(nSumV[[2]] == 0)
  ind = union(ind1,ind2)
  l = length(ind)
  if(l < length(tSumV[[2]])-1){
    if(l>0){
      tDataTest = tSumV[[1]][-ind]/tSumV[[2]][-ind]  
      nDataTest = nSumV[[1]][-ind]/nSumV[[2]][-ind]  
    }else{
      tDataTest = tSumV[[1]]/tSumV[[2]]
      nDataTest = nSumV[[1]]/nSumV[[2]]
    }
    
    t = t.test(tDataTest,nDataTest,alternative=alt,paired=T)
    if(is.nan(t$p.value) | is.infinite(t$p.value)){
      z_score = 0
    }else{
      z_score = qnorm(1 - t$p.value)  
    }
    
  }                       
  z_score
}
plotLipidHeatMapByStage <- function(tData, nData){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = as.vector(detectionLimits[,1])
  lipidArray = c(lipidArray,"PIP","PIP2","PIP3")
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
  p <- lapply(1:length(stages), 
              function(s){
                lapply(1:length(lipidArray),
                       function(i){
                         file.path = paste("D:/project/lipidomics/data/lipid analysis/lipid_class_heat_map_",stages[s],".pdf",sep="")
                         plotLipidHeatMap(stagePatients[[stages[s]]],nData[patientList[[s]],],lipidArray,file.path);
                       })
              })
  p
}
plotLipidHeatMap <- function(tData, nData, lipidArray, file.path){
  #detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  #lipidArray = as.vector(detectionLimits[,1])
  #lipidArray = c(lipidArray,"PIP","PIP2","PIP3")
  #tData = (getPatientSpecieMatrix(TRUE))
  #nData = (getPatientSpecieMatrix(FALSE))
  pr = matrix(nrow = nrow(tData), ncol = length(lipidArray))
  colnames(pr) = lipidArray
  rownames(pr) = rownames(tData)
  tSumV = vector()
  nSumV = vector()
  
  for(i in 1:length(lipidArray)){
    td = filterClassPSM_1(tData,lipidArray[i])
    nd = filterClassPSM_1(nData,lipidArray[i])
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
  #prMat[which(is.nan(prMat))] = min(prMat[which(is.finite(prMat))])-4
  pdf(file = file.path, width=15,height=8, onefile = T)
  Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
            rep("brown",323))
  heatmap.2(t(prMat), main = paste("Lipid Class Heat Map (T:N)"),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
            scale="column", key=T, keysize=0.5, density.info="none",
            trace="none",cexCol=1.2, RowSideColors=Label,
            lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
            lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
  dev.off()
}

getSpeciesNumByClass <- function(){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  lipidArray = detectionLimits[,1]
  tData = (getPatientSpecieMatrix(TRUE))
  data = matrix(nrow = length(lipidArray), ncol = 1)
  rownames(data) = lipidArray
  for(i in 1:length(lipidArray)){
    td = filterClassPSM(tData,lipidArray[i])
    if(is.vector(td)){
      data[i,1] = 1  
    }else{
      data[i,1] = ncol(td)
    }
    
  }
  barplot(data, main="Number of species in each class", 
          ylab="Number of species", xlab = "Lipid class")
}

plotNormalSpeciesValueByStage <- function(){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  setwd("D:/project/lipidomics/data/lipid analysis/Species absolute value/Normal")
  lipidArray = detectionLimits[,1]
  nData = (getPatientSpecieMatrix(FALSE))
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = c("adenoma", "dukes a","dukes b", "dukes c", "dukes d")
  stagePatients =  list()
  patientList = list()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
  }
  specieData = getSpecieData()
  xLabels = c("Adenoma\nn=4","DukesA\nn=9","DukesB\nn=26","DukesC\nn=18","DukesD\nn=14")
  colours = c("white", "darkolivegreen1", "darkolivegreen2", "darkolivegreen3","darkolivegreen4", "darkolivegreen","darkgreen")
  for(i in 1:length(lipidArray)){
    cespecies = rownames(specieData)[which(specieData[,1] == lipidArray[i])]
    numberOfSpecies =  length(cespecies)
    #Set pdf dimensions
    if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
    if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
    pdf(width= w, height=h, file = paste(lipidArray[i],"_AbsoluteValue_Normal.pdf", sep = ""))	
    vdin = c(w,h)
    #Set graphcial grid
    if (numberOfSpecies <4) r<- numberOfSpecies else r<- 4
    if(numberOfSpecies <5) c1<- 1 else c1<-(numberOfSpecies/4)+1
    par(mfcol = c(c1,r),las=2)
    for(j in 1:numberOfSpecies){
      nnd = list()
      nAbsValue = list()
      for(s in 1:length(stages)){
        nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
        if(is.vector(nd)){
          nnd[[s]] = nd
        }else{
          nnd[[s]] = nd[,j]
        }
        nAbsValue[[s]] = nnd[[s]]
      }
      boxplot(nAbsValue,main = cespecies[j], las = 2, names = xLabels, col =colours, ylabel = "Quantity", xlabel = "Stage")
    }
    dev.off()
  }
}
plotTumourSpeciesValueByStage <- function(){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  setwd("D:/project/lipidomics/data/lipid analysis/Species absolute value/Tumour")
  lipidArray = detectionLimits[,1]
  tData = (getPatientSpecieMatrix(TRUE))
  patientData = getPatientData()
  pml = stageMeta(patientData)
  patientMeta = pml$patientMeta
  #split into stages
  stages = c("adenoma", "dukes a","dukes b", "dukes c", "dukes d")
  patientList = list()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    stagePatients[[stages[i]]] = tData[poi,]
  }
  specieData = getSpecieData()
  xLabels = c("Adenoma\nn=4","DukesA\nn=9","DukesB\nn=26","DukesC\nn=18","DukesD\nn=14")
  colours = c("white", "darkolivegreen1", "darkolivegreen2", "darkolivegreen3","darkolivegreen4", "darkolivegreen","darkgreen")
  for(i in 1:length(lipidArray)){
    cespecies = rownames(specieData)[which(specieData[,1] == lipidArray[i])]
    numberOfSpecies =  length(cespecies)
    #Set pdf dimensions
    if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
    if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
    pdf(width= w, height=h, file = paste(lipidArray[i],"_AbsoluteValue_Tumour.pdf", sep = ""))  
    vdin = c(w,h)
    #Set graphcial grid
    if (numberOfSpecies <4) r<- numberOfSpecies else r<- 4
    if(numberOfSpecies <5) c1<- 1 else c1<-(numberOfSpecies/4)+1
    par(mfcol = c(c1,r),las=2)
    for(j in 1:numberOfSpecies){
      ttd = list()
      tAbsValue = list()
      for(s in 1:length(stages)){
        td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
        if(is.vector(td)){
          ttd[[s]] = td
        }else{
          ttd[[s]] = td[,j]
        }
        tAbsValue[[s]] = ttd[[s]]
      }
      boxplot(tAbsValue,main = cespecies[j], las = 2, names = xLabels, col =colours, ylabel = "Quantity", xlabel = "Stage")
    }
    dev.off()
  }
}
getLogRatioSpecies <- function(cespecies){
  numberOfSpecies =  length(cespecies)
  remove = vector()
  tData = (getPatientSpecieMatrix(TRUE))
  nData = (getPatientSpecieMatrix(FALSE))
  patientData = getPatientData()
  stages = c("adenoma", "dukes a","dukes b", "dukes c", "dukes d")
  stagePatients =  list()
  patientList = list()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  for(j in 1:numberOfSpecies){
    ttd = list()
    nnd = list()
    logRatioValue = list()
    logRatioTestPValue = vector()
    for(s in 1:length(stages)){
      td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
      nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
      logRatioValueToRemove = vector()
      if(is.vector(td)){
        ttd[[s]] = td
        nnd[[s]] = nd
      }else{
        ttd[[s]] = td[,j]
        nnd[[s]] = nd[,j]
      }
      ind = which(ttd[[s]] == 0)
      if(length(ind) > 0)
        for(i1 in 1:length(ind))
          logRatioValueToRemove[length(logRatioValueToRemove)+1] = ind[i1]
      ind = which(nnd[[s]] == 0)
      if(length(ind) > 0)
        for(i1 in 1:length(ind))
          logRatioValueToRemove[length(logRatioValueToRemove)+1] = ind[i1]
      if(length(logRatioValueToRemove) > 0){
        ttd[[s]] = ttd[[s]][-logRatioValueToRemove]
        nnd[[s]] = nnd[[s]][-logRatioValueToRemove]
      }
      if(length(ttd[[s]]) > 2){
        t = t.test(ttd[[s]],nnd[[s]], paired = T)
        logRatioValue[[s]] = log10(ttd[[s]]/nnd[[s]])
        logRatioTestPValue[s] = t$p.value
      }
    }
    if(length(logRatioValue) < length(stages)){
      remove[length(remove)+1] = j
    }
    
  }
  if(length(remove) > 0)
    cespecies = cespecies[-remove]
  cespecies
}
plotLogRatioSpeciesByStage <- function(){
  detectionLimits = read.table(file= "D:/project/lipidomics/data/detectionLimit.txt", header = F, sep = "\t")
  setwd("D:/project/lipidomics/data/lipid analysis/Species absolute value/LogRatio")
  lipidArray = detectionLimits[,1]
  lipidArray = lipidArray[-which(lipidArray == "S1P")]
  lipidArray = lipidArray[-which(lipidArray == "SPC")]
  tData = (getPatientSpecieMatrix(TRUE))
  nData = (getPatientSpecieMatrix(FALSE))
  patientData = getPatientData()
  stages = c("adenoma", "dukes a","dukes b", "dukes c", "dukes d")
  stagePatients =  list()
  patientList = list()
  for(i in 1:length(stages)){
    poi = which(patientMeta == stages[i])  
    patientList[[i]] = poi
    stagePatients[[stages[i]]] = tData[poi,]
  }
  specieData = getSpecieData()
  
  for(i in 1:length(lipidArray)){
    cespecies = rownames(specieData)[which(specieData[,1] == lipidArray[i])]
    #cespecies = getLogRatioSpecies(cespecies)
    numberOfSpecies =  length(cespecies)
    #Set pdf dimensions
    if (numberOfSpecies <5) w<- numberOfSpecies * 8 else w<- 32
    if (numberOfSpecies <5) h<- 8 else h<- (((numberOfSpecies/4)+1)*8)
    pdf(width= w, height=h, file = paste(lipidArray[i],"_AbsoluteValue_LogRatio.pdf", sep = ""))  
    vdin = c(w,h)
    #Set graphcial grid
    if (numberOfSpecies <4) r<- numberOfSpecies else r<- 4
    if(numberOfSpecies <5) c1<- 1 else c1<-(numberOfSpecies/4)+1
    par(mfcol = c(c1,r),las=2)
    n = vector()
    for(j in 1:numberOfSpecies){
      ttd = list()
      nnd = list()
      logRatioValue = list()
      logRatioTestPValue = vector()
      logRatioTestColVec = vector()
      xLabels = vector(length = 5)
      for(s in 1:length(stages)){
        td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
        nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
        logRatioValueToRemove = vector()
        if(is.vector(td)){
          ttd[[s]] = td
          nnd[[s]] = nd
        }else{
          ttd[[s]] = td[,j]
          nnd[[s]] = nd[,j]
        }
        ind = which(ttd[[s]] == 0)
        if(length(ind) > 0)
          for(i1 in 1:length(ind))
            logRatioValueToRemove[length(logRatioValueToRemove)+1] = ind[i1]
        ind = which(nnd[[s]] == 0)
        if(length(ind) > 0)
          for(i1 in 1:length(ind))
            logRatioValueToRemove[length(logRatioValueToRemove)+1] = ind[i1]
        if(length(logRatioValueToRemove) > 0){
          ttd[[s]] = ttd[[s]][-logRatioValueToRemove]
          nnd[[s]] = nnd[[s]][-logRatioValueToRemove]
        }
        logRatioTestPValue[s] = 1
        if(length(ttd[[s]]) > 2){
          t = t.test(ttd[[s]],nnd[[s]], paired = T)
          logRatioValue[[s]] = log10(ttd[[s]]/nnd[[s]])
          logRatioTestPValue[s] = t$p.value
        }
        xLabels[s] = paste(stages[s],"\nn = ", length(ttd[[s]]), sep="")
      }
      if(length(logRatioTestPValue) > 0){
        logRatioTestColVec = vector(length = length(logRatioTestPValue))
        for(s in 1:length(stages)){
          if(!is.nan(logRatioTestPValue[s]) && logRatioTestPValue[s] < 0.05){
            logRatioTestColVec[s] = "red"
          }else{
            logRatioTestColVec[s] = "white"  
          }
        }  
      }
      if(length(logRatioValue) == length(stages)){
        boxplot(logRatioValue,main = cespecies[j], las = 2, names = xLabels, col =logRatioTestColVec, ylabel = "Log ratio (T:N)", xlabel = "Stage")
        abline(h=0, col = "blue")  
      }
      
    }
    dev.off()
  }
}
plotLipidSpeciesHeatMapByStage <- function(){
  lipidArray = c("SM","DHCer","Cer","SG","S1P","SPC")
  tData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(TRUE))
  nData = normalisePatientSpecieMatrix(getPatientSpecieMatrix(FALSE))
  #patientData = getPatientData()
  #pml = stageMeta(patientData)
  #patientMeta = pml$patientMeta
  #split into stages
  #stages = unique(patientMeta)
  
  #stagePatients =  list()
  #patientList = list();
  #for(i in 1:length(stages)){
  #  poi = which(patientMeta == stages[i])  
  #  patientList[[i]] = poi
  #  stagePatients[[stages[i]]] = tData[poi,]
  #}
  # for(s in 1:length(stages)){
  for(i in 1:length(lipidArray)){
    #td = filterClassPSM(stagePatients[[stages[s]]],lipidArray[i])
    #nd = filterClassPSM(nData[patientList[[s]],],lipidArray[i])
    td = filterClassPSM(tData,lipidArray[i])
    nd = filterClassPSM(nData,lipidArray[i])
    # divide each species by its lipid class
    td = td/rowSums(td)
    nd = nd/rowSums(nd)
    if(is.vector(td)){
      next
    }else{
      pr = matrix(nrow = nrow(td), ncol = ncol(td))
      pr = td/nd
      colnames(pr) = colnames(td)
      ind = which(pr == 0)
      prMat = log10(pr)
      prMat[ind] = min(prMat[which(is.finite(prMat))])
      prMat[which(is.infinite(prMat))] = max(prMat[which(is.finite(prMat))])
      prMat[which(is.nan(prMat))] = 0
      pdf(file = paste("D:/project/lipidomics/data/lipid analysis/lipid_species_heat_map_",lipidArray[i],".pdf",sep=""), width=15,height=8, onefile = T)
      Label = c(rep("purple",250),rep("orange",250),rep("darkgreen",250),
                rep("brown",323))
      heatmap.2(t(prMat), main = paste("Lipid Species Heat Map (T:N) (", lipidArray[i],")", sep=""),col=redblue(256), dendrogram="both",na.rm = TRUE, na.color=par("bg"),
                scale="column", key=T, keysize=0.5, density.info="none",
                trace="none",cexCol=1.2, RowSideColors=Label,
                lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(2.0,5.0),
                lwid=c(1.5,0.2,2.5,2.5), srtCol=90)
      dev.off()  
    }
    
  }
  #}
  
}
getRFFeatureImportance <- function(nData, tData){
  td1 = as.data.frame(tData)
  nd1 = as.data.frame(nData)
  td1$class = 1 
  nd1$class = 0
  allData = rbind(td1,nd1)	
  #patientData = getPatientData()
  #pml = sampleMeta(rownames(allData))
  #patientMeta = pml$patientMeta
  #toRemove = pml$toRemove
  
  if(length(toRemove)==0){
    PSM = normalisePatientSpecieMatrix(allData)
  }else{
    PSM = normalisePatientSpecieMatrix(allData[-toRemove,])
  }
  finalPSM = PSM
  names(finalPSM) <- make.names(names(finalPSM))
  #split data (PART 4)
  set.seed(4)
  inTrain <- as.vector(createDataPartition(PSM$class, p = 0.75, list = FALSE))
  trainPSM <- finalPSM[inTrain,]
  testPSM <- finalPSM[-inTrain,]
  #trainClassification <- as.factor(patientMeta[inTrain])
  #testClassification <- as.factor(patientMeta[-inTrain])
  
  # training (PART 5, random forest)
  tc <- trainControl(method = "cv", number = 10 )
  set.seed(2)
  nfeatures = ncol(finalPSM)-1
  tunedGrid <- expand.grid(.mtry = seq(from = 10, to = round(nfeatures/2, -1), by = 20 ) )  # try these first and if you don't find an optimum try others
  fit <- train(class~., trainPSM, method = "rf", tuneGrid = tunedGrid, trControl = tc, importance = TRUE)
  
  prediction = predict(fit$finalModel, newdata = testPSM)
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  # run the RFE algorithm
  n = ncol(trainPSM)
  results <- rfe(trainPSM[,1:n-1], trainPSM[,n], sizes=c(1:n-1), rfeControl=control)
  plot(results, type=c("g", "o"))
  # list the chosen features
  predictors(results)
  #rbind(prediction, trainClassification)
  #b= varImp(fit)$importance
  #b = b[order(b[,1], decreasing = TRUE),]
  #colnames(b) = rep("scaled importance",2)
  #b =cbind(b,rownames(b))
  #b= b[,-2]
  #b= b[,c(2,1)]
  #b
}

specieToColID <- function(species,refMat){
  ids = vector(length = length(species))
  for(i in 1:length(species)){
    ids[i] =which(colnames(refMat)== species[i])
  }
  ids
}

BNamesToLHNames<- function(names){
  result = vector(length = length(names))
  for(i in 1:length(names)){
    elems = as.vector(unlist(strsplit(names[i],"-")))
    if(length(elems) == 2){
      result[i] = gsub(":","|",paste(c(elems[2],elems[1]), collapse = " "))
    }else{
      result[i] = names[i]
    }
  }
  result
  
}
isFinite <- function(vec){
  bool = TRUE
  if(is.null(vec))
    bool = FALSE
  if(length(vec) == 0)
    bool = FALSE
  for(i in 1:length(vec)){
    if(!is.finite(vec[i])){
      bool = FALSE
      break
    }
  }
  bool
}
plot.PCAmix <- function(x,axes = c(1, 2), choice = "ind",label=TRUE,
                        coloring.ind=NULL,col.ind=NULL, coloring.var=NULL,
                        lim.cos2.plot=0,lim.contrib.plot=0, posleg="topleft",
                        xlim=NULL,ylim=NULL, cex=1,leg=TRUE,main=NULL,cex.leg=1, ...)
{
  cl<-match.call()
  if (!inherits(x, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
  
  res.pca <-x
  p1 <- res.pca$rec$p1
  p <- res.pca$rec$p
  p2<-res.pca$rec$p2
  m<-nrow(res.pca$levels$coord)
  quanti.coord <- res.pca$quanti$coord
  n<-nrow(res.pca$ind$coord)
  
  eig.axes<-res.pca$eig[axes,1]
  
  if (max(axes) > res.pca$ndim) 
    stop(paste("axes must be between 1 and ", res.pca$ndim, sep = ""))
  
  if (!(choice %in% c("ind", "sqload", "levels", "cor"))) 
    stop("\"choice\" must be either \"ind\",\"sqload\",\"cor\" or \"levels\"")  
  
  if (lim.cos2.plot != 0 & lim.contrib.plot!=0)
    stop("use either \"lim.cos2.plot\" OR \"lim.contrib.plot\"")
  
  if (!is.null(coloring.ind))
  {
    if (choice!="ind")
      warning("use \"coloring.ind\" only if choice=\"ind\"")
  }
  
  if (!is.null(coloring.ind))
  {
    if(!is.factor(coloring.ind) | length(coloring.ind)!=n)
      warning("\"coloring.ind\" must be either NULL or a qualitative variable of length equal to the number of individuals")
  }
  
  if (!is.null(coloring.var))
  {
    if (choice=="ind")
      warning("\"coloring.var\" is not used if choice=\"ind\"")
    if (coloring.var=="type")
    {
      if (choice=="cor" | choice=="levels")
        warning("\"coloring.var\" is not used if choice=\"cor\" or choice=\"levels\"")
    }
  }
  
  if (!is.null(coloring.var))
  {
    if(coloring.var!="type")
      warning("\"coloring.var\" must be either \"NULL\" or \"type\"")
  }
  
  
  
  
  dim1 <- axes[1]
  dim2 <- axes[2]
  
  lab.x <- paste("Dim ", axes[1], " (", signif(res.pca$eig[axes[1], 
                                                           2], 4), " %)", sep = "")
  lab.y <- paste("Dim ", axes[2], " (", signif(res.pca$eig[axes[2], 
                                                           2], 4), " %)", sep = "")
  
  
  if (choice == "ind") {
    if (is.null(main)) 
      main <- "Individuals component map"
    
    coord.ind<-res.pca$ind$coord
    
    if (is.null(xlim)) 
    {
      xmin <- min(coord.ind[, dim1])
      xmax <- max(coord.ind[, dim1])
      xlim <- c(xmin, xmax) * 1.2
    }
    if (is.null(ylim)) 
    {
      ymin <- min(coord.ind[, dim2])
      ymax <- max(coord.ind[, dim2])
      ylim <- c(ymin, ymax) * 1.2
    }
    
    if(is.null(col.ind) | is.null(coloring.ind))
    {
      col.plot.ind<-rep("black",nrow(coord.ind))
    }
    
    if (is.factor(coloring.ind))
    { 
      quali<-coloring.ind
      if (!is.null(col.ind))
      { 
        levels(quali)<-col.ind
        col.plot.ind<-quali
      }
      if(is.null(col.ind))
        col.plot.ind<-as.numeric(quali)
    }
    
    col.plot.ind.total<-col.plot.ind
    
    if(lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      select.ind<-1:nrow(coord.ind)
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot <- lim.cos2.plot
      base.lim <- res.pca$ind$cos2[,axes]
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)    
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$ind$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))
      select.ind<-which(apply(base.lim[,],1,sum)>=lim.plot)
    }
    
    if(length(select.ind)==0)
      warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No individuals can be plotted")
    
    coord.ind<-coord.ind[select.ind, , drop=FALSE]
    col.plot.ind<-col.plot.ind[select.ind]
    
    
    plot(coord.ind[, axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
         ylab = lab.y, pch = 20, col = as.character(col.plot.ind), 
         cex = cex, main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    
    if(length(select.ind)!=0)
    {
      if(leg==T & is.factor(coloring.ind))
        legend(posleg, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), 
               cex =cex.leg)
      
      if (label) 
        text(coord.ind[, axes], labels = rownames(coord.ind), 
             pos = 3, col = as.character(col.plot.ind), cex = cex, 
             ...)
    }
  }
  if (choice == "sqload") {
    if (is.null(main))
      main<-"Squared loadings"
    if (is.null(xlim)) 
    {
      xmax <- max(res.pca$sqload[, dim1])
      xlim <- c(-0.1, xmax * 1.2)
    }
    if (is.null(ylim)) 
    {
      ymax <- max(res.pca$sqload[, dim2])
      ylim <- c(-0.1, ymax * 1.2)
    }
    plot(0, 0, type = "n", xlab = lab.x, ylab = lab.y, xlim = xlim, 
         ylim = ylim, cex = cex,main=main, ...)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    
    if(is.null(coloring.var))
    {
      for (j in 1:nrow(res.pca$sqload)) 
      {
        arrows(0, 0, res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], 
               length = 0.1, angle = 15, code = 2, cex = cex,...)
        if (label) 
        {
          if (res.pca$sqload[j, dim1] > res.pca$sqload[j, dim2]) 
          {
            pos <- 4
          }
          else pos <- 3
          text(res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], labels = rownames(res.pca$sqload)[j], 
               pos = pos, cex = cex, ...)
        }
      }
      
    }
    if (!is.null(coloring.var))
    {
      if (coloring.var=="type")
      {
        for (j in 1:nrow(res.pca$sqload)) 
        {
          col.sq<-rep(c("blue","red"),c(p1,p2))
          arrows(0, 0, res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], 
                 length = 0.1, angle = 15, code = 2, cex = cex, col=col.sq[j], 
                 ...)
          if (label) 
          {
            if (res.pca$sqload[j, dim1] > res.pca$sqload[j, dim2]) {
              pos <- 4
            }
            else pos <- 3
            text(res.pca$sqload[j, dim1], res.pca$sqload[j, dim2], labels = rownames(res.pca$sqload)[j], 
                 pos = pos, cex = cex, col=col.sq[j], ...)
          }
        }
        if (leg==TRUE)
          legend(posleg, legend = c("numerical","categorical"), text.col = c("blue","red"), 
                 cex = cex.leg)
      }
    }
    
  }
  if (choice == "levels") {
    if (is.null(main)) 
      main <- "Levels component map"
    
    if (lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      base.lim<-res.pca$levels$cos2[,axes]
    }
    
    if (lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot<-lim.cos2.plot
      base.lim<-res.pca$levels$cos2[,axes]
    }
    
    if (lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$levels$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))    
    }
    
    color<-rep(1,m)
    
    if (is.null(xlim)) 
    {
      xmin <- min(res.pca$levels$coord[, dim1])
      xmax <- max(res.pca$levels$coord[, dim1])
      xlim <- c(xmin, xmax) * 1.2
    }
    if (is.null(ylim)) 
    {
      ymin <- min(res.pca$levels$coord[, dim2])
      ymax <- max(res.pca$levels$coord[, dim2])
      ylim <- c(ymin, ymax) * 1.2
    }
    
    
    plot(0,0, xlim = xlim, ylim = ylim,
         xlab = lab.x, ylab = lab.y, type="n", cex = cex,main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    nrow.coord.lev <- 0
    if (!is.null(res.pca$levels$coord) ) 
    {
      coord.lev <- res.pca$levels$coord[, axes, drop = FALSE]
      nrow.coord.lev <- nrow(coord.lev)
      
      test.empty.plot<-c()
      for (v in 1:nrow(coord.lev)) 
      {
        if (sum(base.lim[v, ], na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          points(coord.lev[v, 1], coord.lev[v,2], col = color[v],pch=20,cex = cex,...)
          
          if (label) 
          {
            if (abs(coord.lev[v, 1]) > abs(coord.lev[v,2])) 
            {
              if (coord.lev[v, 1] >= 0) 
                pos <- 4
              else pos <- 2
            }
            else {
              if (coord.lev[v, 2] >= 0) 
                pos <- 3
              else pos <- 1
            }
            text(coord.lev[v, 1], y = coord.lev[v, 2], 
                 labels = rownames(coord.lev)[v], pos = pos, 
                 col = color[v], cex = cex)
          }
        }
      }
      if(is.null(test.empty.plot)){
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No level can be plotted")
        return()
      }
    }
  }
  if (choice == "cor") {
    if (is.null(main)) 
      main <- "Correlation circle"
    if(lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      base.lim<-res.pca$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot<-lim.cos2.plot
      base.lim<-res.pca$quanti$cos2[,axes]
    }
    
    if(lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot<-lim.contrib.plot
      base.lim<-res.pca$quanti$contrib[,axes]
      base.lim<-100*(base.lim/sum(eig.axes))     
    }
    
    
    if(is.null(xlim))
    {
      xlim = c(-1.1, 1.1)
    }
    if(is.null(ylim))
    {
      ylim = c(-1.1, 1.1)
    }
    
    col<-rep(1,p1)
    
    plot(0, 0, main = main, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", 
         asp = 1, cex = cex,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    nrow.coord.var <- 0
    if (!is.null(res.pca["quanti"]$quanti$coord) ) 
    {
      coord.var <- res.pca$quanti$coord[, axes, drop = FALSE]
      nrow.coord.var <- nrow(coord.var)
      
      test.empty.plot<-c()      
      for (v in 1:nrow(coord.var)) 
      {
        if (sum(base.lim[v, ] , na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2, col = col[v],cex = cex)
          
          if (label) 
          {
            if (abs(coord.var[v, 1]) > abs(coord.var[v, 
                                                     2])) 
            {
              if (coord.var[v, 1] >= 0) 
                pos <- 4
              else pos <- 2
            }
            else {
              if (coord.var[v, 2] >= 0) 
                pos <- 3
              else pos <- 1
            }
            text(coord.var[v, 1], y = coord.var[v, 2], 
                 labels = rownames(coord.var)[v], pos = pos, 
                 col = col[v], cex = cex)
          }
        }
      }
      if(is.null(test.empty.plot)){
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No variable can be plotted")
        return()
      }
    }
  }
}
