library(BiGGR)
library(rsbml)
library(SBMLR)
data("Recon2")
#Recon2 = rsbml_read("D:/project/lipidomics/data/lipid analysis/hela_hbec/pathway analysis/MODEL1109130000.xml") 
database <- Recon2
listOfEnzyme = c("8611","8612","1606","8694","57104","56994","10390","5337","5338","9791","23761","81490","10400","8398","79888","5320","50487","1040","10423","9489","114971","201164","122618","54675","57104","6609","6610","55512","55627","8560","123099","259230","166929","125981","10715","29956","204219","79603","91012","253782","56848","8611,8612,8613","5168","5297","55361","55300","22908","259230","166929")
listOfEnzyme = unique(listOfEnzyme)
sbml <- buildSBMLFromGenes(listOfEnzyme, database)
model <- database@model
# get all reaction identifiers
ids = sapply(model@reactions, id)
#get all pathways
pathways = getPathwaysForSBML(sbml, database)
m4 = buildSBMLFromPathways(pathways, database)
# build model from reaction ids
m1 = buildSBMLFromReactionIDs(ids, database)
#rsbml_write(sbml, file="myfile.xml")

# plot graph representing the sbml model
# M_dag_hs_c : DG
# M_pchol_hs_c : PC
# M_pe_hs_c : PE
# M_pglyc_hs_c : PG
# M_clpn_hs_c : CL
# M_tag_hs_c : TG
# M_dhcrm_hs_c : DhCer
# M_pgp_hs_c: PGP
# M_lpchol_hs_c : LPC
# M_sphings_c : SG
# M_sphs1p_c : S1P
# M_ps_hs_g : PS
# M_pa_hs_c : PA
# M_CE3481_c : LPA
# M_sphmyln_hs_c : SM
rel.sp <- c("M_12dgr120_c","M_dag_hs_c","M_pchol_hs_c","M_pe_hs_c","M_pglyc_hs_c","M_clpn_hs_c","M_tag_hs_c","M_dhcrm_hs_c","M_pchol_hs_c","M_pgp_hs_c","M_pail34p_hs_c","M_pail45p_hs_c","M_lpchol_hs_c","M_sphings_c","M_sphs1p_c","M_ps_hs_g","M_ps_hs_r","M_pa_hs_c","M_pa_hs_g","M_pa_hs_m","M_CE3481_c","M_sphmyln_hs_c")
relevant.species <- c("M_glc_DASH_D_c", "M_g6p_c", "M_f6p_c", 
                      "M_fdp_c", "M_dhap_c", "M_g3p_c", 
                      "M_13dpg_c", "M_3pg_c", "M_2pg_c", 
                      "M_pep_c", "M_pyr_c")
##Plot model with random rates
rates <- rnorm(length(m1@reactions))
names(rates) <- sapply(m1@reactions, id)
hd <- sbml2hyperdraw(m1, rates=rates, relevant.species=rel.sp, layoutType="dot")
plot(hd)