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

# av number of ind per spp
mean(table(species.v))
sd(table(species.v))

#ind_per_spp <- table (species.v)
#write.csv (ind_per_spp,file=here("output","excel_glm","ind_per_spp.CSV"))

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

# --------------------------------------------------- #
#         Simulated disparity per cell, 
#     as simulated by three evolution models
#                  BM, EB, OU

# Load the concensus

tree <- read.nexus(file=here("data","Sigmodontinae_285speciesTree.tre"))

# Adjusting the names

tree$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",tree $tip.label)

# --------------------------------------- #
#     Empirical disparity per cell

# set the number of interations (equal to the number of simulations in evol models)
niter <- 100

# Run the 'disparity function' that we create (implement Rao's entropy and a randomization of sp. IDs across the rows of trait datasets)
trait.pruned <- size.means.log.v[which(rownames(size.means.log.v) %in% colnames(presab)),]
trait.pruned <- trait.pruned[which(names(trait.pruned) %in% tree$tip.label)]
# occ
comm.pruned <- presab [,which(colnames(presab) %in% names (trait.pruned))]

# phy pruned
# matching between distribution, trait, and phylogenetic datasets
tree.pruned<-  treedata(tree,t(presab))$phy

# empirical rao
RAO_OBS <- disparity_function (occ = comm.pruned,
                               traits= trait.pruned,
                               n_iterations = niter)


save (RAO_OBS, file = here("output_concensus_supp_S3","RAO_OBS_ALL.RData"))

# --------------------------------------------------#
#      now simulate using the complete phylogeny
#                   prune after

# Set the number of computer cores to run analyses in 
# parallel along computer cores
ncores <- 6

# simulate parameters
simul_param_BM <-   fitContinuous(phy=tree.pruned,  
                dat = scale(trait.pruned), 
                model="BM", 
                ncores = ncores)

# ancestral states for each traits
ntraits <- 1
theta<-rep(0,ntraits)

# simulating traits using a Brownian motion model
# 100 simulations per phylogeny
simul_BM <- mvSIM(tree,
	      nsim=niter,
	      model="BM1",
	      param = list(sigma = simul_param_BM$opt$sigsq,
	                   theta=theta,
	                   ntraits=ntraits))


# trait
simul_BM_mean_pruned <- simul_BM[which(rownames(simul_BM) %in% colnames((comm.pruned))),]
# one column per list
simul_BM_mean_pruned<- lapply (seq (1,ncol(simul_BM_mean_pruned)), function (i)
  
  simul_BM_mean_pruned[,i]
  
  )

# Run the disparity function that organizes data and calculate observed and null Rao's entropy

# create a cluster of 'ncores', 
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("simul_BM_mean_pruned", 
                    "comm.pruned",
                    "niter",
                    "disparity_function"))
# load packages in each core
clusterEvalQ(cl,library("SYNCSA"))
# run
RAO_BM <- parLapply (cl,simul_BM_mean_pruned, function (phy)	
  disparity_function (occ = (comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

# save
save(RAO_BM, file=here("output_concensus_supp_S3","RAO_BM_ALL.Rdata"))

## ============================================== #
# simulating traits using an early burst model
# simulate parameters
simul_param_EB <- fitContinuous(phy=tree.pruned,  
                dat = scale(trait.pruned), 
                model="EB", 
                ncores = ncores)

# simulate
simul_EB<-mvSIM(tree,
        nsim=niter,
        model="EB",
        param = list(sigma = simul_param_EB$opt$sigsq,
                     theta=theta,
                     beta = simul_param_EB$opt$a))

# trait
simul_EB_mean_pruned <- simul_EB[which(rownames(simul_EB) %in% colnames((comm.pruned))),]
# one column per list
simul_EB_mean_pruned<- lapply (seq (1,ncol(simul_EB_mean_pruned)), function (i)
  
  simul_EB_mean_pruned[,i]
  
)
# Run the disparity function that organizes data and calculate observed and null Rao's entropy

# create a cluster of 'ncores', 
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("simul_EB_mean_pruned", 
                    "comm.pruned",
                    "niter",
                    "disparity_function"))
# load packages in each core
clusterEvalQ(cl,library("SYNCSA"))
# run
RAO_EB <- parLapply (cl,simul_EB_mean_pruned, function (phy)	
  disparity_function (occ = (comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

# save
save (RAO_EB, file = here("output_concensus_supp_S3","RAO_EB_ALL.RData"))

#==================================================
# simulating traits using a Ornstein-Uhlenbeck model
# 100 simulations per phylogeny
# simulating traits using a Ornstein-Uhlenbeck model
# 100 simulations per phylogeny
simul_param_OU <- fitContinuous(phy=tree.pruned,  
                dat = scale(trait.pruned), 
                model="OU", 
                ncores = ncores)

# simulate traits
simul_OU <-mvSIM(tree,
  		  nsim=niter,
          model="OU1",
  
                param = list(sigma = simul_param_OU$opt$sigsq,
                             alpha=simul_param_OU$opt$alpha,
                             theta = theta))
 

# trait
simul_OU_mean_pruned <- simul_OU[which(rownames(simul_OU) %in% colnames((comm.pruned))),]
# one column per list
simul_OU_mean_pruned<- lapply (seq (1,ncol(simul_OU_mean_pruned)), function (i)
  
  simul_OU_mean_pruned[,i]
  
)
# run disparity function
# create a cluster of 'ncores', 
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("simul_OU_mean_pruned", 
                    "comm.pruned",
                    "niter",
                    "disparity_function"))
# load packages in each core
clusterEvalQ(cl,library("SYNCSA"))
# run
RAO_OU <- parLapply (cl,simul_OU_mean_pruned, function (phy)	
  disparity_function (occ = (comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save (RAO_OU, file = here("output_concensus_supp_S3","RAO_OU_ALL.RData"))
# save parameters
save(simul_param_BM,
     simul_param_EB,
     simul_param_OU,
     file = here("output_concensus_supp_S3","params_fitcontinuous.RData"))
# ------------------------------------------------------- #
#      mean phylogenetic distance between species (MPD) 
#     and associated standardized effect size (SES.MPD)

niter<-1000

# observed MPD
obs.mpd <- apply (comm.pruned,1,mpd.function,
                    dist.tree=cophenetic(tree.pruned))
  
# replicate shuffle per phylogeny 
rep.shuff.phy <- replicate(niter, mpd.shuff(tree.pruned,
	                                    my.sample.matrix=(comm.pruned)))
  
# get statistics  
statistics.phy <- data.frame (averageMPD = apply (rep.shuff.phy,1,mean),
			      sdMPD = apply (rep.shuff.phy,1,sd),
			       obsMPD = obs.mpd
			)
  
# calculate SES
statistics.phy$SES.MPD <- (statistics.phy$obsMPD - statistics.phy$averageMPD)/statistics.phy$sdMPD
  
# save
save (statistics.phy, file=here ("output_concensus_supp_S3","mpd_results_ALL.RData"))

