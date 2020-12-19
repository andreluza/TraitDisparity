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

# --------------------------------------- #
#     Empirical disparity per cell

# set the number of interations (equal to the number of simulations in evol models)
niter <- 1000

# Run the 'disparity function' that we create (implement Rao's entropy and a randomization of sp. IDs across the rows of trait datasets)
RAO_OBS <- funcao_disparidade (occ = as.matrix (presab),
                               traits= shape.v.2d.means,
                               n_iterations = niter)

save (RAO_OBS, file = here("output","RAO_OBS.RData"))

# --------------------------------------------------- #
#         Simulated disparity per cell, 
#     as simulated by three evolution models
#                  BM, EB, OU

# Load the fully resolved phylogenties

tree_list <- tree <- read.nexus(file=here("data","Sigmodontinae_413species100Trees.trees"))

# Adjusting the names

tree_list <- lapply (tree_list, function (i)
	{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# matching between distribution, trait, and phylogenetic datasets
tree.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$phy)
comm.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$data)[[1]]

# --------------------------------------------------#
#      ANALYSES

# Set the number of computer cores to run analyses in 
# parallel along computer cores
ncores <- 7

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
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save(RAO_BM, file=here("output","RAO_BM.Rdata"))

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
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)
# save
save (RAO_EB, file = here("output","RAO_EB.RData"))

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
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

save (RAO_OU, file = here("output","RAO_OU.RData"))

# ------------------------------------------------------- #
#      mean phylogenetic distance between species (MPD) 
#     and associated standardized effect size (SES.MPD)

# function to calculate mpd
mpd.shuff <- function (tree,my.sample.matrix) {
  
  shuff.tree <- tipShuffle (tree)
  mpd(my.sample.matrix, cophenetic (shuff.tree))
  
}

# replicating the function niter times per phylogeny
null.mpdf <- lapply (tree.pruned, function (tree) {
  
  # observed MPD
  obs.mpd <- mpd (t(comm.pruned),cophenetic (tree))
  
  # replicate shuffle per phylogeny 
  rep.shuff.phy <- replicate(niter, mpd.shuff(tree,
        my.sample.matrix=t(comm.pruned)))
  
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

# save
save (null.mpdf, file=here ("output","mpd_results.RData"))


#################################################################################
#### Simulações sob outros modelos evolutivos ####
#### EB
require(mvMORPH)
simul_EB<-mvSIM(tree,nsim=10000,model="EB",param = list(sigma = 1, beta = -0.5))
mean_EB<-apply(simul_EB,1,mean)
sd_EB<-apply(simul_EB,1,sd)

# EB média e desvio
RESULTADOS.RAO.EB<-matrix(NA,length(presab[,1]),dim(simul_EB)[2])
for(i in 1:dim(simul_EB)[2]){
  org<-organize.syncsa(comm=presab2,traits=as.data.frame(simul_EB[,i]))
  raoBM<-rao.diversity(org$community,traits=org$traits)
  RESULTADOS.RAO.EB[,i]<-as.vector(raoBM$FunRao)
} 
RESULTADOS.RAO.EB
mean.rao.EB<-apply(RESULTADOS.RAO.EB,1,mean)
sd.rao.EB<-apply(RESULTADOS.RAO.EB,1,sd)

cor(rao$FunRao,mean.rao.EB)
ES_EB<-rao$FunRao - mean.rao.EB
SES_EB<- ES_BM / sd.rao.EB

plot(rao$FunRao,ses.mpd$mpd.obs.z)
plot(mean.rao.BM,ses.mpd$mpd.obs.z)
plot(SES_BM,ses.mpd$mpd.obs.z)
plot(mean.rao.EB,ses.mpd$mpd.obs.z)
plot(SES_EB,ses.mpd$mpd.obs.z)


#### OU
require(mvMORPH)
simul_OU<-mvSIM(tree,nsim=10000,model="OU1",param = list(sigma = 1,alpha=1))
mean_OU<-apply(simul_OU,1,mean)
sd_OU<-apply(simul_OU,1,sd)

# OU média e desvio
RESULTADOS.RAO.OU<-matrix(NA,length(presab[,1]),dim(simul_OU)[2])
for(i in 1:dim(simul_OU)[2]){
  org<-organize.syncsa(comm=presab2,traits=as.data.frame(simul_OU[,i]))
  raoBM<-rao.diversity(org$community,traits=org$traits)
  RESULTADOS.RAO.OU[,i]<-as.vector(raoBM$FunRao)
} 



#### Save & Write ####
save.image(file = "my_work_space.RData")
load("my_work_space.RData")
write.csv(cbind(rao$FunRao,mean.rao.NULO,sd.rao.NULO,
            mean.rao.BM,sd.rao.BM,SES_BM,
            mean.rao.EB,sd.rao.EB,SES_EB,
            mean.rao.OU,sd.rao.OU,SES_OU,
            SES_NULO,SES_NULO_BM,SES_NULO_EB,SES_NULO_OU,
            ses.mpd$mpd.obs.z,ses.mpd$mpd.obs),"RESULTADOS.RAO.csv")


################################################################################
#### Gráficos dispersão ####
tabela<-read.table("SAMGrid_Neotropical_MainDataBase - Disparity - 2.txt",h=T)
names(tabela)

plot(tabela$media~tabela$Species_Richness_.228.)

# Disparidade observada e riqueza
plot(tabela$rao.obs~tabela$Species_Richness_.228.,
     xlab="Species Richness",ylab="Observed Disparity",pch=21,bg="grey",cex=1.4)
# Disparidade média NULA e riqueza
plot(tabela$mean.rao.NULO~tabela$Species_Richness_.228.,
     xlab="Species Richness",ylab="Random Disparity",pch=21,bg="grey",cex=1.4)
# SES.NULO e SES.MPD
plot(tabela$SES_NULO~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity",pch=21,bg="grey",cex=1.4)
abline(lm(tabela$SES_NULO~tabela$SES_MPD),lty="dashed")
summary(lm(tabela$SES_NULO~tabela$SES_MPD))
text(-4.8,1,"R²=0.42")
# SES Composto e SES.MPD
plot(tabela$SES_NULO~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity",pch=21,bg="grey",cex=1.4,ylim=c(-9,2))
points(tabela$SES_NULO_BM~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity Brownian",pch=21,bg="blue",cex=1.4)
points(tabela$SES_NULO_EB~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity Early-Burst",pch=21,bg="red",cex=1.4)
points(tabela$SES_NULO_OU~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity Ornstein-Uhlenback",pch=21,bg="green",cex=1.4)

# SES Composto e SES.MPD usando Linhas de tendência
plot(tabela$SES_NULO~tabela$SES_MPD,
     xlab="SES MPD",ylab="SES Disparity",pch=21,type="n",ylim=c(-9,2))
abline(lm(tabela$SES_NULO~tabela$SES_MPD),col="grey")
abline(lm(tabela$SES_NULO_BM~tabela$SES_MPD),col="blue")
abline(lm(tabela$SES_NULO_EB~tabela$SES_MPD),col="red")
abline(lm(tabela$SES_NULO_OU~tabela$SES_MPD),col="green")

# CI
mod <- lm(tabela$SES_NULO~tabela$SES_MPD)
preds <- predict(mod, interval = 'confidence')
# intervals
lines(tabela$SES_MPD, preds[ ,3], lty = 'dotted', col = 'gray')
lines(tabela$SES_MPD, preds[ ,2], lty = 'dotted', col = 'gray')





### ACABOU ####
#####################################################
# Coeficientes
summary(lm(tabela$rao.obs~tabela$mean.rao.EB))
summary(lm(tabela$mean.rao.NULO~tabela$mean.rao.EB))

# Teste CI
x<-rnorm(100)
y<-rnorm(100)
plot(x,y)
preds1<-predict(lm(y~x),interval='confidence',level=0.8)
abline(lm(y~x))
lines(x,preds1[,3])
lines(x,preds1[,2])


# Testes
CI_NULO_EB<-matrix(NA,length(RESULTADOS.RAO.EB[,1]),length(RESULTADOS.RAO.EB[1,])) 
for(i in 1:length(RESULTADOS.RAO.EB[1,])){
  SES_NULO_EB_Random<-(RESULTADOS.RAO.EB[,i] - mean.rao.NULO) / sd.rao.NULO
  CI_NULO_EB[,i]<-SES_NULO_EB_Random
} 
CI_NULO_EB[,1:5]

plot(CI_NULO_EB[,1]~ses.mpd$mpd.obs.z,type="n")
for(i in 1:length(CI_NULO_EB[1,])){
  abline(lm(CI_NULO_EB[,i]~ses.mpd$mpd.obs.z))
}

abline(lm(CI_NULO_EB[,3]~ses.mpd$mpd.obs.z),col="red")
plot(CI_NULO_EB[,3]~ses.mpd$mpd.obs.z)

