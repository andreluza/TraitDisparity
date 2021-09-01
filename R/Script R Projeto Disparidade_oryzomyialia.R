# ------------------------# 

# load packages
source ("./R/packages.R")
source ("./R/functions.R")

#----------------------------------------#
#       Load data of skull shape 
#             Shape ventral

# load '.tps' file with ventral view
tps.ventral<-readland.tps(here ("data","Sigmodontinae.ventral.dig.tps"),specID = "ID", readcurves = FALSE)
##dim(tps.ventral)

# load specimen information
ventral.listed<-read.table(here ("data","Sigmodontinae.ventral.listed.txt"),h=T)
# unique species 
species.v<-ventral.listed[,4]

# Removing species from islands
# Remove Nephelomys pirrensis & Sigmodon zanjonensis & Brucepattersonius_sp (no shapefile)
# Remove other species that are out of analyzed spatial extent: Rheomys underwoodi, Rheomys mexicanus, Oryzomys palustris, Nectomys p. palmipes, Nesoryzomys swarthi, Nesoryzomys indefessus, Nesoryzomys fernandinae, Nesoryzomys darwini, Aegialomys galapagoensis
# Ending up with 228 species
rm_spp <- c("Nephelomys_pirrensis","Sigmodon_zanjonensis","Aegialomys_galapagoensis","Brucepattersonius_sp","Rheomys_underwoodi","Nesoryzomys_darwini","Rheomys_mexicanus","Oryzomys_palustris","Nesoryzomys_fernandinae","Nectomys_palmipes","Nesoryzomys_swarthi","Nesoryzomys_indefessus")
tps.ventral<-tps.ventral[,,which(ventral.listed$Species_Patton2015 %in% rm_spp == F)] 

# GPA Ventral
## Procrustes analysis to get independent shape and size data
gpa.ventral<-gpagen(tps.ventral)
names(gpa.ventral)
size.v<-gpa.ventral$Csize
shape.v<-gpa.ventral$coords

# load GPA bilateral data
load(here ("data","shape.v.RData"))

# forma média por espécie (ventral)
ventral.listed.pruned<-read.table(here ("data","Sigmodontinae.ventral.listed.pruned.txt"),h=T)
#fix(ventral.listed.pruned)
species.v<-ventral.listed.pruned[,4]
#str(species.v)

# transform shape data into array
shape.v.2d<-two.d.array(shape.v)
shape.v.2d.means<-rowsum(shape.v.2d,species.v)/as.vector(table(species.v))
shape.means.v<-arrayspecs(shape.v.2d.means,dim(shape.v)[1],dim(shape.v)[2])
#shape.means.v
#str(shape.means.v)

# transform size data into a vector in log scale
size.means.v<-rowsum(size.v,species.v)/as.vector(table(species.v))
hist(log(size.means.v))
size.means.log.v<-log(size.means.v)

# ------------------------------- #
# working just with Oryzomyialia
non_oryz <- ventral.listed.pruned$Species_Patton2015 [which(ventral.listed.pruned$Tribe 
                                          %in% c("Ichthyomyini", "Sigmodontini"))]

shape_oryz <- shape.v.2d.means[which(rownames(shape.v.2d.means) %in% non_oryz == F),]

# -----------------------------------------------#
#           Load spatial data
#       Incidence of 228 species in 1770 cells

presab_original <-read.table(here ("data","PresAbs_228sp_Neotropical_MainDataBase_Ordenado.txt"),h=T)
# remove sites with 2 or less sp
presab <- presab_original [rowSums (presab_original) > 3,]
## exclude species absent in the remaining cells
presab <- presab [,which(colSums(presab)>0)]

# ----------------------------- # 
##  Geographic coordinates data

longlat<-read.table(here ("data","Lon-lat-Disparity.txt"),h=T)
longlat  <- longlat [rowSums (presab_original) > 3,]

# ------------------------------- #
#     matching occurrence, trait, 
#        and phylogenetic data 

# Load the fully resolved phylogenties

tree_list <- tree <- read.nexus(file=here("data","Sigmodontinae_413species100Trees.trees"))

# Adjusting the names

tree_list <- lapply (tree_list, function (i)
{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# matching between distribution, trait, and phylogenetic datasets
# phylogeny
tree.pruned<-lapply (tree_list, function (phy) 
  treedata(phy,shape_oryz)$phy)
# matching trait
trait.pruned<-lapply (tree_list, function (phy) 
  treedata(phy,shape_oryz)$data)
trait.pruned <- as.data.frame(trait.pruned[[1]])
# community
comm.pruned <- presab [,which(colnames(presab) %in% rownames(trait.pruned))]

# --------------------------------------- #
#     Empirical disparity per cell

# set the number of interations (equal to the number of simulations in evol models)
niter <- 1000

# Run the 'disparity function' that we create (implement Rao's entropy and a randomization of sp. IDs across the rows of trait datasets)
RAO_OBS <- funcao_disparidade (occ =  as.matrix(comm.pruned),
                               traits= (trait.pruned),
                               n_iterations = niter)

save (RAO_OBS, file = here("output","RAO_OBS_ORYZ.RData"))

# --------------------------------------------------- #
#         Simulated disparity per cell, 
#     as simulated by three evolution models
#                  BM, EB, OU

# Set the number of computer cores to run analyses in 
# parallel along computer cores
ncores <- detectCores()-1

# simulating traits using a Brownian motion model
# 100 simulations per phylogeny
simul_BM<-lapply (tree.pruned, function (phy) 
	mvSIM(phy,nsim=niter,model="BM1",param = list(sigma = 1)))

# Get average values

simul_BM_mean <- lapply(simul_BM, function (phy)
  apply(phy,1,mean))

# Run the disparity function that organizes data and calculate observed and null Rao's entropy

# create a cluster of 'ncores', 
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("simul_BM_mean", "comm.pruned","niter","funcao_disparidade"))
# load packages in each core
clusterEvalQ(cl,library("SYNCSA"))
# run
RAO_BM <- parLapply (cl,simul_BM_mean, function (phy)	
  funcao_disparidade (occ = ((comm.pruned)),
                      traits= (phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save(RAO_BM, file=here("output","RAO_BM_ORYZ.Rdata"))

# simulating traits using a early burst model
# 100 simulations per phylogeny
simul_EB<-lapply (tree.pruned, function (phy) 
  mvSIM(phy,nsim=niter,model="EB",param = list(sigma = 1, beta = -0.5)))

# Get average values
simul_EB_mean <- lapply(simul_EB, function (phy)
  apply(phy,1,mean))

# Run the disparity function that organizes data and calculate observed and null Rao's entropy
cl <- makeCluster (ncores)
clusterExport(cl, c("simul_EB_mean", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_EB <- parLapply (cl, simul_EB_mean, function (phy)	
  funcao_disparidade (occ = ((comm.pruned)),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save (RAO_EB, file = here("output","RAO_EB_ORYZ.RData"))

# simulating traits using a Ornstein-Uhlenbeck model
# 100 simulations per phylogeny

simul_OU <-lapply (tree.pruned, function (phy) 
  mvSIM(phy,nsim=niter,model="OU1",param = list(sigma = 1,alpha=1)))

simul_OU_mean <- lapply(simul_OU, function (phy)
  apply(phy,1,mean))

# run disparity function
cl <- makeCluster (ncores)
clusterExport(cl, c("simul_OU_mean", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_OU <- parLapply (cl,simul_OU_mean, function (phy)	
  funcao_disparidade (occ = ((comm.pruned)),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

save (RAO_OU, file = here("output","RAO_OU_ORYZ.RData"))

# ------------------------------------------------------- #
#      mean phylogenetic distance between species (MPD) 
#     and associated standardized effect size (SES.MPD)

niter<-1000

# run disparity function
cl <- makeCluster (ncores)
clusterExport(cl, c("tree.pruned", "comm.pruned","niter",
                    "mpd.function","mpd.shuff"))
clusterEvalQ(cl,library("SYNCSA"))
clusterEvalQ(cl,library("picante"))

# replicating the function niter times per phylogeny
null.mpdf <- parLapply (cl, tree.pruned, function (tree) {
  
  # observed MPD
  obs.mpd <- apply (comm.pruned,1,mpd.function,
                    dist.tree=cophenetic(tree))
  
  # replicate shuffle per phylogeny 
  rep.shuff.phy <- replicate(niter, mpd.shuff(tree,
                                              my.sample.matrix=(comm.pruned)))
  
  # get statistics  
  statistics.phy <- data.frame (
    averageMPD = apply (rep.shuff.phy,1,mean),
    sdMPD = apply (rep.shuff.phy,1,sd),
    obsMPD = obs.mpd)
  
  # calculate SES
  statistics.phy$SES.MPD <- (statistics.phy$obsMPD - statistics.phy$averageMPD)/statistics.phy$sdMPD
  
  ; # return
  
  statistics.phy
  
}
)

stopCluster (cl)

# save
save (null.mpdf, file=here ("output","mpd_results_ORYZ.RData"))

