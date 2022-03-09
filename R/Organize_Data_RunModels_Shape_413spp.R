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
# load(here ("data","shape.v.RData"))

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

tree_list <- tree <- read.nexus(file=here("data","Sigmodontinae_413species100Trees.trees"))

# Adjusting the names
tree_list <- lapply (tree_list, function (i)
	{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# --------------------------------------- #
#     Empirical disparity per cell

# set the number of interations (equal to the number of simulations in evol models)
niter <- 100

# Run the 'disparity function' that we create (implement Rao's entropy and a randomization of sp. IDs across the rows of trait datasets)
trait.pruned <- shape.v.2d.means[which(rownames(shape.v.2d.means) %in% colnames(presab)),]
trait.pruned <- trait.pruned[which(rownames(trait.pruned) %in% tree_list[[1]]$tip.label),]
# occ
comm.pruned <- presab [,which(colnames(presab) %in% rownames (trait.pruned))]

# phy pruned
# prune
# matching between distribution, trait, and phylogenetic datasets
tree.pruned<-lapply (tree_list, function (phy) 
  treedata(phy,t(presab))$phy)

# empirical rao
RAO_OBS <- disparity_function (occ = comm.pruned,
                               traits= scale(trait.pruned),
                               n_iterations = niter)

save (RAO_OBS, file = here("output","RAO_OBS_ALL_multivariate.RData"))

# --------------------------------------------------#
#      now simulate using the complete phylogeny
#                   prune after

# Set the number of computer cores to run analyses in 
# parallel along computer cores
ncores <- 6

# simulate parameters
simul_param_BM <- lapply (tree.pruned, function (i) 
  
  mvgls(scale(trait.pruned)~1, 
        tree=i,  
        model="BM", 
        penalty="RidgeArch",
        REML=T,
        method = "Mahalanobis")
)

# ancestral states for each traits
ntraits <- ncol(trait.pruned)
theta<-rep(0,ntraits)

# simulating traits using a Brownian motion model
# 100 simulations per phylogeny
simul_BM<-lapply (seq(1,length(tree_list)), function (phy) 
  
	mvSIM(tree_list[[phy]],
	      nsim=niter,
	      model="BM1",
	      param = list(sigma = simul_param_BM[[phy]]$sigma$Pinv,
	                   theta=theta)))

# Get average values

simul_BM_mean <- lapply (simul_BM, function (i)
  
  Reduce("+",i)/length(i))

# trait
simul_BM_mean_pruned <- lapply(simul_BM_mean, function (i)
  i[which(rownames(i) %in% colnames((comm.pruned))),])

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
save(RAO_BM, file=here("output","RAO_BM_ALL_multivariate.Rdata"))
save(simul_param_BM, file=here("output","params_BM_multivariate.Rdata"))
