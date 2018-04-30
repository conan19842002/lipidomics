library(RMySQL)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(graph)
source("/home/data/lipidomics/data/lipidomicLib.r")
require(XLConnect)
require(ggplot2)
require(gridExtra)
require(data.table)
library(NMF)
library(xlsx)
library(qgraph)
library(igraph)
library(plyr)
library(XML)

fa_acyl = c(paste0(seq(16,24,1),":0"), "18:3", "18:4", "20:5","22:6","18:2","20:2", 
            "20:3","20:4","22:4","16:1","18:1","20:1","22:1", "22:2", "24:1","22:3","22:5","24:2")

getCdpDag <- function(){
  a = c()
  b= c()
  for(i in 1:length(fa_acyl)){
    s = unlist(strsplit(fa_acyl[i],":"))
    a = c(a, as.numeric(s[1]))
    b = c(b, as.numeric(s[2]))
  }
  a1 = expand.grid(rep(list(a),2))
  b1 = expand.grid(rep(list(b),2))
  c1 = paste0(a1,":",b1)
  c1 = unique(c1)
  c1
}
getMonoLysoCl <- function(){
  a = c()
  b= c()
  for(i in 1:length(fa_acyl)){
    s = unlist(strsplit(fa_acyl[i],":"))
    a = c(a, as.numeric(s[1]))
    b = c(b, as.numeric(s[2]))
  }
  
  a1 = apply(expand.grid(rep(list(a),3)),1,sum)
  b1 = apply(expand.grid(rep(list(b),3)),1,sum)
  c1 = paste0(a1,":",b1)
  c1 = unique(c1)
  c1
}

mono_lyso_cl <- getMonoLysoCl()
cdp_dag <- getCdpDag()
# return a dataframe which contains reactions of the same layers
getReactionByLayer <- function(lipids, reactions){
  sReactions = matrix(ncol = 2)
  firstRun = TRUE
  for(i in 1:nrow(reactions)){
    reactant = reactions[i,1]
    product = reactions[i,2]	
    rSpecies = lipids[which(grepl(paste0("-",reactant), lipids))]
    pSpecies = lipids[which(grepl(paste0("-",product), lipids))]
   
    rFAs = getFAs(rSpecies)
    pFAs =  getFAs(pSpecies)
    sharedFAs = which(rFAs %in% pFAs)
    if(length(sharedFAs)>0){
      for(j in 1:length(sharedFAs)){	
        fai = sharedFAs[j]
        fa = rFAs[fai]
        if(firstRun){
          sReactions[1,] = c(paste0(fa, "-",reactant),paste0(pFAs[which(pFAs==fa)],"-",product))
          firstRun = FALSE
        }else{
          sReactions = rbind(sReactions, c(paste0(fa,"-",reactant),paste0(pFAs[which(pFAs==fa)],"-",product)))
          
        }
      }  
    }
    
  }
  df = data.frame(reactant = sReactions[,1], product = sReactions[,2])
  df
}

getAcylTable <- function(lp){
  acyl_dgs = c()
  bond_dgs = c()
  spl = unlist(strsplit(lp,"-"))
  l = length(spl)
  acyls = unlist(strsplit(spl[seq(1,l,2)],":"))
  l1 = length(acyls)
  df = data.frame(acyl = as.numeric(acyls[seq(1,l1,2)]), bond = as.numeric(acyls[seq(2,l1,2)]))
  df
}
getReactionNotSameLayer <- function(lipids, is_dg_tg = T, is_mg_dg = T, is_pg_cl = T, is_pc_cl = T, is_lyso_pa = T,
                                    is_lyso_pe = T, is_lyso_pc = T, is_lyso_ps = T, is_lyso_pi = T, is_lyso_pg = T){
  
  reactant = c()
  product = c()
  
  # FA + DG -> TG
  if(is_dg_tg){

    dg = lipids[which(grepl("-DG",lipids))]
    dg_df = getAcylTable(dg)
    
    tg = lipids[which(grepl("-TG",lipids))]
    tg_df = getAcylTable(tg)
  
    tg_df = sort(tg_df)
    for(i in 1:nrow(dg_df)){
      acyl = dg_df[i,1]
      bond = dg_df[i,2]
      re = paste0(acyl,":",bond,"-DG")
      a = tg_df$acyl- acyl
      b = tg_df$bond - bond
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      if(length(id)>0){
        sub_tg_df = tg_df[id,]
        pro = paste0(sub_tg_df$acyl, ":",sub_tg_df$bond,"-TG")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }
  }
  
  # FA + MG -> DG
  if(is_mg_dg){
    mg = lipids[which(grepl("-MG",lipids))]
    mg_df = getAcylTable(mg)
    dg_df = sort(dg_df)
    for(i in 1:nrow(mg_df)){
      acyl = mg_df[i,1]
      bond = mg_df[i,2]
      re = paste0(acyl,":",bond,"-MG")
      a = dg_df$acyl- acyl
      b = dg_df$bond - bond
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      
      if(length(id)>0){
        sub_dg_df = dg_df[id,]
        pro = paste0(sub_dg_df$acyl, ":",sub_dg_df$bond,"-DG")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }
  }
  
  # CDP-DAG + PG -> CL
  if(is_pg_cl){
   
    pg = lipids[which(grepl("-PG",lipids))]
    pg_df = getAcylTable(pg)
    cl = lipids[which(grepl("-CL",lipids))]
    cl_df = getAcylTable(cl)
    cl_df = sort(cl_df)
    for(i in 1:nrow(pg_df)){
      acyl = pg_df[i,1]
      bond = pg_df[i,2]
      re = paste0(acyl,":",bond,"-PG")
      a = cl_df$acyl- acyl
      b = cl_df$bond - bond
      cdpdag = paste0(a,":",b)
      id = which(cdpdag %in% cdp_dag)
      
      if(length(id)>0){
        sub_cl_df = cl_df[id,]
        pro = paste0(sub_cl_df$acyl, ":",sub_cl_df$bond,"-CL")
        reactant = c(reactant, rep(re,length(id)))
        product = c(product, pro)
      }
    }
  }
  
  # PC + MonoLysoCL -> CL
  if(is_pc_cl){
    pc =lipids[which(grepl("-PC",lipids))]
    cl = lipids[which(grepl("-CL",lipids))]
    pc_df = getAcylTable(pc)
    cl_df = getAcylTable(cl)
    for(i in 1:nrow(pc_df)){
      a1 = pc_df[i,1]
      b1 = pc_df[i,2]
      re = paste0(a1,":",b1,"-PC")
      a = cl_df$acyl- acyl
      b = cl_df$bond - bond
      lysocl = paste0(a,":",b)
      id = which(lysocl %in% mono_lyso_cl)
    
      if(length(id)>0){
        sub_cl_df = cl_df[id,]
        pro = paste0(sub_cl_df$acyl, ":",sub_cl_df$bond,"-CL")
        reactant = c(reactant, rep(re,length(id)))
        product = c(product, pro)
      }
    }  
  }
  
  # LPA + acyl CoA -> PA
  if(is_lyso_pa){
    
    pa =lipids[which(grepl("-PA",lipids))]
    pa_df = getAcylTable(pa)
    lpa =lipids[which(grepl("-LPA",lipids))]
    lpa_df = getAcylTable(lpa)
    for(i in 1:nrow(lpa_df)){
      a1 = lpa_df[i,1]
      b1 = lpa_df[i,2]
      re = paste0(a1,":",b1,"-LPA")
      a = pa_df$acyl- a1
      b = pa_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_pa_df = pa_df[id,]
        pro = paste0(sub_pa_df$acyl, ":",sub_pa_df$bond,"-PA")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
  }
  # LPC + FA -> PC
  if(is_lyso_pc){
    pc =lipids[which(grepl("-PC",lipids))]
    pc_df = getAcylTable(pc)
    lpc =lipids[which(grepl("-LPC",lipids))]
    lpc_df = getAcylTable(lpc)
    for(i in 1:nrow(lpc_df)){
      a1 = lpc_df[i,1]
      b1 = lpc_df[i,2]
      re = paste0(a1,":",b1,"-LPC")
      a = pc_df$acyl- a1
      b = pc_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_pc_df = pc_df[id,]
        pro = paste0(sub_pc_df$acyl, ":",sub_pc_df$bond,"-PC")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
    
  }
  
  # LPE + FA -> PE
  if(is_lyso_pe){
    pe =lipids[which(grepl("-PE",lipids))]
    pe_df = getAcylTable(pe)
    lpe =lipids[which(grepl("-LPE",lipids))]
    lpe_df = getAcylTable(lpe)
    for(i in 1:nrow(lpe_df)){
      a1 = lpe_df[i,1]
      b1 = lpe_df[i,2]
      re = paste0(a1,":",b1,"-LPE")
      a = pe_df$acyl- a1
      b = pe_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_pe_df = pe_df[id,]
        pro = paste0(sub_pe_df$acyl, ":",sub_pe_df$bond,"-PE")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
    
  }
  
  # LPS + FA -> PS
  if(is_lyso_ps){
    ps =lipids[which(grepl("-PS",lipids))]
    ps_df = getAcylTable(ps)
    lps =lipids[which(grepl("-LPS",lipids))]
    lps_df = getAcylTable(lps)
    for(i in 1:nrow(lps_df)){
      a1 = lps_df[i,1]
      b1 = lps_df[i,2]
      re = paste0(a1,":",b1,"-LPS")
      a = ps_df$acyl- a1
      b = ps_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_ps_df = ps_df[id,]
        pro = paste0(sub_ps_df$acyl, ":",sub_ps_df$bond,"-PS")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
    
  }
  
  # LPI + FA -> PI
  if(is_lyso_pi){
    pi =lipids[which(grepl("-PI",lipids))]
    pi_df = getAcylTable(pi)
    lpi =lipids[which(grepl("-LPI",lipids))]
    lpi_df = getAcylTable(lpi)
    for(i in 1:nrow(lpi_df)){
      a1 = lpi_df[i,1]
      b1 = lpi_df[i,2]
      re = paste0(a1,":",b1,"-LPI")
      a = pi_df$acyl- a1
      b = pi_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_pi_df = pi_df[id,]
        pro = paste0(sub_pi_df$acyl, ":",sub_pi_df$bond,"-PI")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
    
  }
  
  # LPG + FA -> PG
  if(is_lyso_pg){
    pg =lipids[which(grepl("-PG",lipids))]
    pg_df = getAcylTable(pg)
    lpg =lipids[which(grepl("-LPG",lipids))]
    lpg_df = getAcylTable(lpg)
    for(i in 1:nrow(lpg_df)){
      a1 = lpg_df[i,1]
      b1 = lpg_df[i,2]
      re = paste0(a1,":",b1,"-LPG")
      a = pg_df$acyl- a1
      b = pg_df$bond - b1
      fa = paste0(a,":",b)
      id = which(fa %in% fa_acyl)
      
      if(length(id)>0){
        sub_pg_df = pg_df[id,]
        pro = paste0(sub_pg_df$acyl, ":",sub_pg_df$bond,"-PG")
        reactant = c(reactant, c(rep(re,length(id)), pro))
        product = c(product, c(pro, rep(re,length(id))))
      }
    }  
    
    
  }
  df1 = data.frame(reactant = reactant, product = product)
  df1
}
getReactionNotSameLayer0 <- function(df){
  df = read.csv("D:/project/lipidomics/data/lipid analysis/ppi/hmdb_id.csv", header = TRUE)
  # FA + DG -> TG
  all_dg_id = which(grepl("-DG",df$abbreviation))
  cdp_dg_id = which(grepl("-CDP-DG",df[all_dg_id,]$abbreviation))
  dg_id = all_dg_id[-cdp_dg_id]
  dg = as.character(df[dg_id,]$abbreviation)
  dg_df = getAcylTable(dg)
  
  tg = as.character(df[which(grepl("-TG",df$abbreviation)),]$abbreviation)
  tg_df = getAcylTable(tg)
  
  mg = as.character(df[which(grepl("-MG",df$abbreviation)),]$abbreviation)
  mg_df = getAcylTable(mg)
  
  reactant = c()
  product = c()
  tg_df = sort(tg_df)
  for(i in 1:nrow(dg_df)){
    acyl = dg_df[i,1]
    bond = dg_df[i,2]
    re = paste0(acyl,":",bond,"-DG")
    id = which(tg_df$acyl >= acyl + 12 & tg_df$bond >= bond)
    if(length(id)>0){
      sub_tg_df = tg_df[id,]
      pro = paste0(sub_tg_df$acyl, ":",sub_tg_df$bond,"-TG")
      reactant = c(reactant, c(rep(re,length(id)), pro))
      product = c(product, c(pro, rep(re,length(id))))
    }
  }
  # FA + MG -> DG
  dg_df = sort(dg_df)
  for(i in 1:nrow(mg_df)){
    acyl = mg_df[i,1]
    bond = mg_df[i,2]
    re = paste0(acyl,":",bond,"-MG")
    id = which(dg_df$acyl >= acyl + 12 & dg_df$bond >= bond)
    if(length(id)>0){
      sub_dg_df = dg_df[id,]
      pro = paste0(sub_dg_df$acyl, ":",sub_dg_df$bond,"-DG")
      reactant = c(reactant, c(rep(re,length(id)), pro))
      product = c(product, c(pro, rep(re,length(id))))
    }
  }
  # CDP-DAG + PG -> CL
  al_pg_id =  which(grepl("-PG",df$abbreviation))
  pgp_id = which(grepl("-PGP",df[al_pg_id,]$abbreviation))
  pg_id = al_pg_id[-pgp_id]
  pg = as.character(df[pg_id,]$abbreviation)
  pg_df = getAcylTable(pg)
  cl = as.character(df[which(grepl("-CL",df$abbreviation)),]$abbreviation)
  cl_df = getAcylTable(cl)
  cdp_dg = as.character(df[all_dg_id[cdp_dg_id],]$abbreviation)
  cdp_dg = gsub("CDP-DG","CDP_DG",cdp_dg)
  cdp_dg_df = getAcylTable(cdp_dg)
  cl_df = sort(cl_df)
  for(i in 1:nrow(cdp_dg_df)){
    a1 = cdp_dg_df[i,1]
    b1 = cdp_dg_df[i,2]
    for(j in 1:nrow(pg_df)){
      a2 = pg_df[j,1]
      b2 = pg_df[j,2]
      re = paste0(a2,":",b2,"-PG")
      id = which(cl_df$acyl == a1+a2 & cl_df$bond == b1 + b2)
      if(length(id)>0){
        sub_cl_df = cl_df[id,]
        pro = paste0(sub_cl_df$acyl, ":",sub_cl_df$bond,"-CL")
        reactant = c(reactant, rep(re,length(id)))
        product = c(product, pro)
      }
    }
  }
  # PC + LysoCL -> CL
  id = which(grepl("-PC",df$abbreviation))
  a1 =df[id,] 
  a = a1$abrreviation
  pc = as.character(a)
  pc_df = getAcylTable(pc)
  for(i in 1:nrow(pc_df)){
    a1 = pc_df[i,1]
    b1 = pc_df[i,2]
    re = paste0(a1,":",b1,"-PC")
    id = which(cl_df$acyl >= a1 + 12 & cl_df$bond >= b1)
    if(length(id)>0){
      sub_cl_df = cl_df[id,]
      pro = paste0(sub_cl_df$acyl, ":",sub_cl_df$bond,"-CL")
      reactant = c(reactant, rep(re,length(id)))
      product = c(product, pro)
    }
  }
  df1 = data.frame(reactant = reactant, product = product)
  write.csv(df1,"D:/project/lipidomics/data/lipid analysis/ppi/reaction_not_same_layer.csv")
}
getSwissLipidId <- function(lipids){
  hmdbIdL = read.csv("D:/project/lipidomics/data/lipid analysis/ppi/hmdb_id.csv", header = TRUE)
  sdf = read.SDFset("D:/project/lipidomics/data/lipid analysis/ppi/LMSDFDownload6Dec16/LMSDFDownload6Dec16FinalAll.sdf")
}