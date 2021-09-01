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

ind_per_spp <- table (species.v)
write.csv (ind_per_spp,file=here("output","excel_glm","ind_per_spp.CSV"))

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

# Load the fully resolved phylogenties

tree_list <- tree <- read.nexus(file=here("data","tree - Sigmodontinae_285speciesTree.tre"))


# --------------------------------------- #
#     Empirical disparity per cell

# set the number of interations (equal to the number of simulations in evol models)
niter <- 1000

# Run the 'disparity function' that we create (implement Rao's entropy and a randomization of sp. IDs across the rows of trait datasets)
trait.pruned <- shape.v.2d.means[which(rownames(shape.v.2d.means) %in% colnames(presab)),]
trait.pruned <- trait.pruned[which(rownames(trait.pruned) %in% tree_list$tip.label),]
# occ
comm.pruned <- presab [,which(colnames(presab) %in% rownames (trait.pruned))]

# phy pruned
# prune
# matching between distribution, trait, and phylogenetic datasets
tree.pruned<- treedata(tree_list,t(presab))$phy

# empirical rao
RAO_OBS <- funcao_disparidade (occ = comm.pruned,
                               traits= trait.pruned,
                               n_iterations = niter)

save (RAO_OBS, file = here("output_sensitivity_supp_S3","RAO_OBS_ALL.RData"))


# --------------------------------------------------#
#      now simulate using the complete phylogeny
#                   prune after

# Set the number of computer cores to run analyses in 
# parallel along computer cores
ncores <- 7

# simulating traits using a Brownian motion model
# 100 simulations per phylogeny
sigmas <- seq(0,2,0.5)
simul_BM<-lapply (sigmas, function (s) 
	mvSIM(tree_list,nsim=niter,model="BM1",param = list(sigma = s)))

# Get average values

simul_BM_mean <- lapply(simul_BM, function (phy)
  apply(phy,1,mean))

# trait
simul_BM_mean_pruned <- lapply(simul_BM_mean, function (i)
  i[which(names(i) %in% colnames((comm.pruned)))])

# Run the disparity function that organizes data and calculate observed and null Rao's entropy

# create a cluster of 'ncores', 
cl <- makeCluster (ncores)
# load data and functions in each core 
clusterExport(cl, c("simul_BM_mean_pruned", "comm.pruned","niter","funcao_disparidade"))
# load packages in each core
clusterEvalQ(cl,library("SYNCSA"))
# run
RAO_BM <- parLapply (cl,simul_BM_mean_pruned, function (phy)	
  funcao_disparidade (occ = (comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save(RAO_BM, file=here("output_sensitivity_supp_S3","RAO_BM_ALL.Rdata"))

# simulating traits using a early burst model
# 100 simulations per phylogeny
betas<-seq(-1,1,0.5)
simul_EB<-lapply (sigmas, function (s) 
  
            lapply(betas, function (b)
  
             mvSIM(tree_list,nsim=niter,model="EB",param = list(sigma = s, beta = b))
            
            )
          )

# Get average values
simul_EB_mean <- lapply(simul_EB, function (s)
  
                      lapply(s, function (b)
                    
                          apply(b,1,mean)
                          )
                      )

# trait
simul_EB_mean_pruned <- lapply(simul_EB_mean, function (s)
  
                          lapply(s, function (b)
  
                            b[which(names(b) %in% colnames((comm.pruned)))]
                          
                          )
                      )

# Run the disparity function that organizes data and calculate observed and null Rao's entropy
cl <- makeCluster (ncores)
clusterExport(cl, c("simul_EB_mean_pruned", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_EB <- parLapply (cl, simul_EB_mean_pruned, function (s)	
  
            lapply(s, function (b)
  
                  funcao_disparidade (occ = (comm.pruned),
                                  traits= as.data.frame(b),
                                  n_iterations = niter)
                  )
          )

stopCluster (cl)
# save
save (RAO_EB, file = here("output_sensitivity_supp_S3","RAO_EB_ALL.RData"))

# simulating traits using a Ornstein-Uhlenbeck model
# 100 simulations per phylogeny
alphas <- seq(0.5,3,0.5)
simul_OU <-lapply (sigmas, function (s) 
  
            lapply (alphas, function (a)
    
               mvSIM(tree_list,nsim=niter,model="OU1",param = list(sigma = s,alpha=a))
        )
  )

  # get average
simul_OU_mean <- lapply(simul_OU, function (s)
  
                      lapply(s, function (a)
  
                          apply(a,1,mean)
          )
)

# trait
simul_OU_mean_pruned <- lapply(simul_OU_mean, function (s)
  
                          lapply(s, function (a)
                            
                              a[which(names(a) %in% colnames((comm.pruned)))]
                              
                              )
)

# run disparity function
cl <- makeCluster (ncores)
clusterExport(cl, c("simul_OU_mean_pruned", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_OU <- parLapply (cl,simul_OU_mean_pruned, function (s)	
  
              lapply(s, function (a)
  
                funcao_disparidade (occ = (comm.pruned),
                      traits= as.data.frame(a),
                      n_iterations = niter)
                )
)

stopCluster (cl)

save (RAO_OU, file = here("output_sensitivity_supp_S3","RAO_OU_ALL.RData"))

# save used parameter values
save(sigmas,betas,alphas,file=here ("output_sensitivity","parameters.RData"))

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
null.mpdf.function <- function (tree,comm) {
  
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
  
  #return
  return (  statistics.phy)
  
}
# run
null.mpdf <- null.mpdf.function(tree=tree.pruned,
                                comm=comm.pruned)


# save
save (null.mpdf, file=here ("output_sensitivity_supp_S3","mpd_results_ALL.RData"))

