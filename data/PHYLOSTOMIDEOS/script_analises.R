#rm(list=ls())
#require (raster)
#require(clipr)# clear_clip()
#require(rgdal)
#require(sp)
#require(rNOMADS)
#require(rgeos)
#setwd("/Users/gferraz/Downloads/data")
#grb2files <- lapply(list.files(path = "./", pattern = "*.sfluxgrbf.grb2$"), readGDAL)
#grb2files <- lapply(list.files(path = "./", pattern = "cdas1.20151025.sfluxgrbf.grb2"), readGDAL)
#grb2files <- readGDAL("cdas1.20181231.sfluxgrbf.grb2")
#class (grb2files)
#image (grb2files[1])
#summary (grb2files)
#toraster <- raster(grb2files)
#values(toraster) <- seq (1,1716,1)
#image(toraster)
#require(maps)
#par (mar = c(0,0,0,0))
#map(regions="Brazil")
#points(coordinates (toraster)[,1]-360, coordinates (toraster)[,2])
## lembrar e-mail

require (sp)
require(spdep)
require(rgeos)
require(rgdal)

setwd ("C:/Users/Acer/OneDrive/Simposio_masto/dados_script/R PROJETO DISPARIDADE FINAL/PHYLOSTOMIDEOS/Distribui��o")
bat_distribution <- readOGR (dsn=getwd(),
                  layer = "gridded_new_world_only_phyllo")

#plot(bat_distribution)
#points (coordinates (bat_distribution[1] [which(bat_distribution[100]@data[,1]==1 ),]),cex=0.1,col="red")


PA_bat <- lapply(seq (1,dim(bat_distribution)[2]), function (i)
  bat_distribution [i]@data[,1])

PA_matrix<-do.call(cbind,PA_bat)[,-c(1:5)]
colnames(PA_matrix) <- names(bat_distribution)[-c(1:5)]

PA_matrix <- PA_matrix [rowSums (PA_matrix) > 3,]
## excluir entao especies que nao ocorrem nas celulas remanescentes
PA_matrix <- PA_matrix [,which(colSums(PA_matrix)>0)]

## filogenias

setwd("C:/Users/Acer/OneDrive/Simposio_masto/dados_script/R PROJETO DISPARIDADE FINAL/PHYLOSTOMIDEOS/doi_10.5061_dryad.s533p__v2/Bayesian_majority-rule_trees")

require(ape)
bat_phylogenies <- lapply (list.files (getwd()), read.nexus)

# resolver abreviacoes
abb_names <- lapply (seq (1,length(bat_phylogenies)), function (i)
	lapply (seq (1, length(bat_phylogenies[[i]]$tip)), function (k)
		paste(substr(strsplit (bat_phylogenies[[i]]$tip, "_")[[k]],1,4)[1],
			substr(strsplit (bat_phylogenies[[i]]$tip, "_")[[k]],1,4)[2],sep="_")))

abb_names <- lapply (abb_names, function (i) do.call (rbind, i))

### substituir os nomes nas filo

bat_phylogenies_names <- lapply (seq (1,length(bat_phylogenies)), function (i) {
	bat_phylogenies[[i]]$tip.label <- abb_names [[i]];
	bat_phylogenies})

bat_phylogenies_names <- lapply (seq (1,10), function (i)
	bat_phylogenies_names [[i]][[i]])

########## resolver as politomias (segundo Rangel et al. 2015)
### resolver as politomias da filogenia com mais sp - eh a quinta
#install.packages("remotes")
#remotes::install_github("davidnipperess/PDcalc",force=T)

require(phytools)
require("PDcalc")

ultrametric_bat_tree <- force.ultrametric(bat_phylogenies_names[[5]], method="nnls")
solved_ultrametric_bat_tree <- bifurcatr(ultrametric_bat_tree, runs = 100)

## pruning phylogenies
require(geiger)
tree.pruned<-lapply (solved_ultrametric_bat_tree, function (phy) 
	treedata(phy,t(PA_matrix))$phy)

tree.pruned <- lapply (tree.pruned, function (i) 
	drop.tip(i, tip =which(duplicated(i$tip)), trim.internal = TRUE))

comm.pruned <-lapply (tree.pruned, function (phy) 
	treedata(phy,t(PA_matrix))$data)[[1]]

### landmarks
# Carregar arquivo .tps com vista ventral
require(geomorph)
setwd("C:/Users/Acer/OneDrive/Simposio_masto/dados_script/R PROJETO DISPARIDADE FINAL/PHYLOSTOMIDEOS")
tps.ventral<-readland.tps("phyllo_ventral_final.TPS",specID = "ID", readcurves = FALSE)
dim(tps.ventral)

# Carregar lista (classificadores) vista ventral
ventral.listed<-read.csv("phyllo-ventral-final.csv",h=T,sep=";")
#fix(ventral.listed)
species.v<-ventral.listed[,2]
#str(species.v)

# GPA Ventral
gpa.ventral<-gpagen(tps.ventral)
names(gpa.ventral)
size.v<-gpa.ventral$Csize
shape.v<-gpa.ventral$coords

# GPA bilateral
#source("bilat.symmetry.r") # Call old bilat.symmetry function; I changed the name to: bilat.symm
#pairs.matrix<-matrix(seq(2,11),ncol=2,byrow=T)
#ind<-c(1:dim(tps.ventral)[3]) # vetor indivíduos
#b.s<-bilat.symmetry(tps.ventral,ind=ind,object.sym=T,land.pairs=pairs.matrix)
#shape.v<-b.s$symm.shape # componente simétrico da forma

pdf("morpho_bat.pdf")
plotAllSpecimens(shape.v)
dev.off()

#save(shape.v,file="bat.shape.v.RData")
#save(size.v, file="bat.size.v.RData")

load("shape.v.RData")
load("size.v.RData")

shape.v.2d<-two.d.array(shape.v)
shape.v.2d.means<-rowsum(shape.v.2d,species.v)/as.vector(table(species.v))
shape.means.v<-arrayspecs(shape.v.2d.means,dim(shape.v)[1],dim(shape.v)[2])
shape.means.v
str(shape.means.v)

size.means.v<-rowsum(size.v,species.v)/as.vector(table(species.v))
size.means.v
hist(log(size.means.v))
size.means.log.v <-log(size.means.v)

##

nomes_sub_LMK <- strsplit (rownames(shape.v.2d.means), "_")
nomes_sub_LMK <- lapply (nomes_sub_LMK, function (i)
	paste(substr(i,1,4)[1],substr(i,1,4)[2],sep="_"))

#colnames(shape.means.v[1,,]) == rownames(size.means.v)

rownames(shape.v.2d.means) <- unlist(nomes_sub_LMK)
rownames(size.means.log.v) <- unlist(nomes_sub_LMK)
#dimnames(shape.means.v)[3][[1]] [which (dimnames(shape.means.v)[3][[1]] %in% tree.pruned[[1]]$tip)]

require(picante)

prunned.size <-lapply (tree.pruned, function (phy) 
	treedata(phy,size.means.log.v)$data)[[1]]
prunned.phy <-lapply (tree.pruned, function (phy) 
	treedata(phy,size.means.log.v)$phy)
prunned.comm <-lapply (prunned.phy , function (phy) 
	treedata(phy,t(PA_matrix))$data)[[1]]
colnames (prunned.comm) <- seq(1, dim(prunned.comm)[2])

#prunned.phy <-lapply (prunned.phy, function (phy) 
#	treedata(phy,prunned.comm)$phy)
prunned.size <-lapply (prunned.phy , function (phy) 
	treedata(phy,size.means.log.v)$data)[[1]]
prunned.shape <-lapply (prunned.phy , function (phy) 
	treedata(phy,shape.v.2d.means)$data)[[1]]

#################################################################################
### Disparidade Observada Empírica por Sítio ####
require(SYNCSA)

### Disparidade Nula (permutação) por Sítio ####
## funcao serve para pegar os atributos observados e simulados e calcular rao nulo e sd nulo em uma filogenia

funcao_disparidade <- function (occ, traits,n_iterations) {
	data<- as.matrix(traits) # shape.v.2d.means
	colnames(data)<-seq(1,ncol(data))
	require(SYNCSA)
	## organizar os dados observados
	organized_obs_data <- organize.syncsa (comm=as.matrix(occ),traits=data,check.comm = T)
	# calcular Rao e extrai-lo
	rao_observado <-rao.diversity(organized_obs_data$community,traits=organized_obs_data$traits)$FunRao
	
	### SES nulo para ol observado
	
	RESULTADOS.RAO <- lapply(seq(1,n_iterations), function (i) {
		  row.names(organized_obs_data$traits) <- sample(row.names(organized_obs_data$traits))
		  org.perm <- organize.syncsa(comm=organized_obs_data$community,traits=organized_obs_data$traits)
		  rao.perm <- rao.diversity(org.perm$community,traits=org.perm$traits)
		  rao_aleat <-rao.perm$FunRao; rao_aleat
		})

	mean.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,mean)
	sd.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,sd)

	DF_RESULTADOS <- data.frame (Observado = rao_observado, med_nulo=mean.rao.NULO,sd_nulo=sd.rao.NULO,
			SES = ((rao_observado - mean.rao.NULO)/sd.rao.NULO))

	return (DF_RESULTADOS)

}


RAO_OBS_SHAPE <- funcao_disparidade (occ = t(prunned.comm),
			traits= prunned.shape,n_iterations = 500)

RAO_OBS_SIZE <- funcao_disparidade (occ = t(prunned.comm),
			traits= prunned.size,n_iterations = 500)

save (RAO_OBS_SHAPE,RAO_OBS_SIZE, file="RAO_OBS_BATS.RData")

### SIMULACOES

## carregar pacote para simular atributos
require(mvMORPH)

## definir o numero de iteracoes do taxa shuffle
n_iterations <- 500

# 100 simulações aplicadas em cada uma das 100 filogenias
simul_BM <-lapply (prunned.phy, function (phy) 
	mvSIM(phy,nsim=100,model="BM1",param = list(sigma = 1)))

simul_BM_mean <- lapply(simul_BM, function (phy)
	apply(phy,1,mean))

## organizar os dados - repetir para todas as filogenias 
require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_BM_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_BM <- parLapply (cl,simul_BM_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

# 100 simulações DE EB aplicadas em cada uma das 100 filogenias
simul_EB<-lapply (prunned.phy, function (phy) 
	mvSIM(phy,nsim=100,model="EB",param = list(sigma = 1, beta = -0.5)))

simul_EB_mean <- lapply(simul_EB, function (phy)
	apply(phy,1,mean))

require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_EB_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_EB <- parLapply (cl, simul_EB_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

# 100 simulações DE OU aplicadas em cada uma das 100 filogenias
simul_OU <-lapply (prunned.phy, function (phy) 
	mvSIM(phy,nsim=100,model="OU1",param = list(sigma = 1,alpha=1)))

simul_OU_mean <- lapply(simul_OU, function (phy)
	apply(phy,1,mean))


require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_OU_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_OU <- parLapply (cl,simul_OU_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

save (RAO_BM, RAO_EB,RAO_OU, file="resultados_mod_evolutivos_BAT.RData")
#load ("resultados_MPD.RData")

#######################
### MPD e SES.MPD ####
require(picante)
require(geiger)

## organizar os dados - repetir para todas as filogenias 
require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("prunned.phy", "prunned.comm"))
clusterEvalQ(cl,library("picante"))

est.ses.mpd<- parLapply (cl, prunned.phy, function (phy)
	ses.mpd(t(prunned.comm),cophenetic(phy),null.model="taxa.labels", runs=999))

stopCluster (cl)

#### ses para cada filogenia
### extrair o Z
mean.ses.mpd <- apply (do.call(cbind,lapply (est.ses.mpd, function (i) 
	i$mpd.obs.z)),1,mean)
sd.ses.mpd <- apply (do.call(cbind,lapply (est.ses.mpd, function (i) 
	i$mpd.obs.z)),1,sd)

save(mean.ses.mpd,sd.ses.mpd, file="resultados_MPD_bat.RData")

## SES do Rao Nulo para cranio
SES_RAO_OBS_SHAPE <- RAO_OBS_SHAPE$SES
SES_RAO_OBS_SHAPE <- SES_RAO_OBS_SHAPE  [-which(is.na(SES_RAO_OBS_SHAPE )==T)]
#SES_RAO_OBS <- SES_RAO_OBS [SES_RAO_OBS > -8000]
## SES do Rao Nulo para tamanho do cranio
SES_RAO_OBS_SIZE <- RAO_OBS_SIZE$SES
SES_RAO_OBS_SIZE <- SES_RAO_OBS_SIZE [-which(is.na(SES_RAO_OBS_SIZE)==T)]

## SES do Rao BM
SES_RAO_BM <- apply(do.call(cbind,lapply(RAO_BM, function (i) i$SES)),1,mean)
SES_RAO_BM <- SES_RAO_BM [-which(is.na(SES_RAO_BM)==T)]

## SES do Rao EB
SES_RAO_EB <- apply(do.call(cbind,lapply(RAO_EB, function (i) i$SES)),1,mean)
SES_RAO_EB <- SES_RAO_EB [-which(is.na(SES_RAO_EB)==T)]

## SES do Rao OU
SES_RAO_OU <- apply(do.call(cbind,lapply(RAO_OU, function (i) i$SES)),1,mean)
SES_RAO_OU <- SES_RAO_OU [-which(is.na(SES_RAO_OU)==T)]

mean.ses.mpd <- mean.ses.mpd[-which(is.na(mean.ses.mpd)==T)]


### SKULL SHAPE
pdf ("plot_forma_cranio_morcegos.pdf",height=6,width=6)
plot(SES_RAO_OBS_SHAPE ~ mean.ses.mpd, xlab="SES of mean phylogenetic distance\nbetween co-occuring species",pch=19,
	ylab= "SES of Rao diversity",main="Skull shape")
abline(lm(SES_RAO_OBS_SHAPE ~mean.ses.mpd),lwd=2)
points(SES_RAO_BM~mean.ses.mpd,col="blue",pch=19)
points(SES_RAO_EB~mean.ses.mpd,col="red",pch=19)
#points(SES_RAO_LB~mean.ses.mpd,col="purple",pch=19)
#points(SES_RAO_LB2~mean.ses.mpd,col="gray",pch=19)
points(SES_RAO_OU~mean.ses.mpd,col="green",pch=19)
abline (v=0,h=0,lwd=2,lty=2,col="gray50")

legend ('topleft', legend= c("Observed", "Brownian motion", "Early burst","Ornstein-Uhlenbeck"),
	pch=19,col=c("black", "blue", "red", "green"), bty='n',cex=1)

dev.off()

## SKULL SIZE

pdf ("plot_tamanho_cranio_morcegos.pdf",height=6,width=6)
plot(SES_RAO_OBS_SIZE ~ mean.ses.mpd, xlab="SES of mean phylogenetic distance\nbetween co-occuring species",pch=19,
	ylab= "SES of Rao diversity",main="Skull size")
abline(lm(SES_RAO_OBS_SIZE ~mean.ses.mpd),lwd=2)
points(SES_RAO_BM~mean.ses.mpd,col="blue",pch=19)
points(SES_RAO_EB~mean.ses.mpd,col="red",pch=19)
#points(SES_RAO_LB~mean.ses.mpd,col="purple",pch=19)
#points(SES_RAO_LB2~mean.ses.mpd,col="gray",pch=19)
points(SES_RAO_OU~mean.ses.mpd,col="green",pch=19)
abline (v=0,h=0,lwd=2,lty=2,col="gray50")

legend ('topleft', legend= c("Observed", "Brownian motion", "Early burst","Ornstein-Uhlenbeck"),
	pch=19,col=c("black", "blue", "red", "green"), bty='n',cex=1)
dev.off()


#### simular com arvores simuladas


simul_tree <- lapply (seq(1,100), function (i) 
		rcoal (n=length(prunned.phy [[1]]$tip.label),
		rooted=T,
		tip.label=prunned.phy [[1]]$tip.label,
		br="coalescent"))


### SIMULACOES

## carregar pacote para simular atributos
require(mvMORPH)

## definir o numero de iteracoes do taxa shuffle
n_iterations <- 500

# 100 simulações aplicadas em cada uma das 100 filogenias
simul_BM <-lapply (simul_tree, function (phy) 
	mvSIM(phy,nsim=100,model="BM1",param = list(sigma = 1)))

simul_BM_mean <- lapply(simul_BM, function (phy)
	apply(phy,1,mean))

## organizar os dados - repetir para todas as filogenias 
require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_BM_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_BM <- parLapply (cl,simul_BM_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

# 100 simulações DE EB aplicadas em cada uma das 100 filogenias
simul_EB<-lapply (simul_tree, function (phy) 
	mvSIM(phy,nsim=100,model="EB",param = list(sigma = 1, beta = -0.5)))

simul_EB_mean <- lapply(simul_EB, function (phy)
	apply(phy,1,mean))

require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_EB_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_EB <- parLapply (cl, simul_EB_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

# 100 simulações DE OU aplicadas em cada uma das 100 filogenias
simul_OU <-lapply (simul_tree, function (phy) 
	mvSIM(phy,nsim=100,model="OU1",param = list(sigma = 1,alpha=1)))

simul_OU_mean <- lapply(simul_OU, function (phy)
	apply(phy,1,mean))


require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_OU_mean", "prunned.comm","n_iterations","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_OU <- parLapply (cl,simul_OU_mean, function (phy)	
		funcao_disparidade (occ = t(prunned.comm),
			traits= as.data.frame(phy),
			n_iterations = n_iterations))

stopCluster (cl)

save (RAO_BM, RAO_EB,RAO_OU, file="resultados_mod_evolutivos_filo_simul.RData")


#######################
### MPD e SES.MPD ####
require(picante)
require(geiger)

## organizar os dados - repetir para todas as filogenias 
require(parallel)
cl <- makeCluster (4)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_tree", "prunned.comm"))
clusterEvalQ(cl,library("picante"))

est.ses.mpd<- parLapply (cl, simul_tree, function (phy)
	ses.mpd(t(prunned.comm),cophenetic(phy),null.model="taxa.labels", runs=999))

stopCluster (cl)

#### ses para cada filogenia
### extrair o Z
mean.ses.mpd <- apply (do.call(cbind,lapply (est.ses.mpd, function (i) 
	i$mpd.obs.z)),1,mean)
sd.ses.mpd <- apply (do.call(cbind,lapply (est.ses.mpd, function (i) 
	i$mpd.obs.z)),1,sd)

save(mean.ses.mpd,sd.ses.mpd, file="resultados_MPD_simul_tree.RData")


## SES do Rao BM
SES_RAO_BM <- apply(do.call(cbind,lapply(RAO_BM, function (i) i$SES)),1,mean)
SES_RAO_BM <- SES_RAO_BM [-which(is.na(SES_RAO_BM)==T)]

## SES do Rao EB
SES_RAO_EB <- apply(do.call(cbind,lapply(RAO_EB, function (i) i$SES)),1,mean)
SES_RAO_EB <- SES_RAO_EB [-which(is.na(SES_RAO_EB)==T)]

## SES do Rao OU
SES_RAO_OU <- apply(do.call(cbind,lapply(RAO_OU, function (i) i$SES)),1,mean)
SES_RAO_OU <- SES_RAO_OU [-which(is.na(SES_RAO_OU)==T)]

mean.ses.mpd <- mean.ses.mpd[-which(is.na(mean.ses.mpd)==T)]


plot(SES_RAO_BM ~ mean.ses.mpd, xlab="SES of mean phylogenetic distance\nbetween co-occuring species",pch=19,
	ylab= "SES of Rao diversity",main="Any trait generated by...",col="blue",
	ylim=c(-0.5, 0.5))
abline(lm(SES_RAO_BM ~mean.ses.mpd),lwd=2,col="blue")
points(SES_RAO_EB~mean.ses.mpd,col="red",pch=19)
abline(lm(SES_RAO_EB ~mean.ses.mpd),lwd=2,col="red")
points(SES_RAO_OU~mean.ses.mpd,col="green",pch=19)
abline(lm(SES_RAO_OU ~mean.ses.mpd),lwd=2,col="green")
abline (v=0,h=0,lwd=2,lty=2,col="gray50")

legend ('topleft', legend= c("Brownian motion", "Early burst","Ornstein-Uhlenbeck"),
	pch=19,col=c("blue", "red", "green"), bty='n',cex=1)



