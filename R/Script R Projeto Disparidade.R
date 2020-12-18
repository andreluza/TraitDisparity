rm(list=ls())
require(geomorph)
require(here)

#################################################################################
#### Carregar Forma do Crânio ####
#### Shape ventral
# Carregar arquivo .tps com vista ventral
tps.ventral<-readland.tps(here ("data","Sigmodontinae.ventral.dig.tps"),specID = "ID", readcurves = FALSE)
#dim(tps.ventral)

# Carregar lista (classificadores) vista ventral
ventral.listed<-read.table(here ("data","Sigmodontinae.ventral.listed.txt"),h=T)
#fix(ventral.listed)
species.v<-ventral.listed[,4]
#str(species.v)

# Eliminar Nephelomys pirrensis & Sigmodon zanjonensis (sem shapefile) & Brucepattersonius_sp 
# Eliminar 9 (fora do recorte espacial): Rheomys underwoodi, Rheomys mexicanus, Oryzomys palustris, Nectomys p. palmipes, Nesoryzomys swarthi, Nesoryzomys indefessus, Nesoryzomys fernandinae, Nesoryzomys darwini, Aegialomys galapagoensis
# Ending with 228 species
rm_spp <- c("Nephelomys_pirrensis","Sigmodon_zanjonensis","Aegialomys_galapagoensis","Brucepattersonius_sp","Rheomys_underwoodi","Nesoryzomys_darwini","Rheomys_mexicanus","Oryzomys_palustris","Nesoryzomys_fernandinae","Nectomys_palmipes","Nesoryzomys_swarthi","Nesoryzomys_indefessus")
tps.ventral<-tps.ventral[,,which(ventral.listed$Species_Patton2015 %in% rm_spp == F)] 
dim(tps.ventral)

# GPA Ventral
gpa.ventral<-gpagen(tps.ventral)
names(gpa.ventral)
size.v<-gpa.ventral$Csize
shape.v<-gpa.ventral$coords
# GPA bilateral
#source("bilat.symmetry.r") # Call old bilat.symmetry function; I changed the name to: bilat.symm
#pairs.matrix<-matrix(c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,24,22,28,23,29,19,25,20,26,21,27,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,55,56),nrow=26,ncol=2,byrow=T)
#ind<-c(1:dim(tps.ventral)[3]) # vetor indivíduos
#b.s<-bilat.symmetry(tps.ventral,ind=ind,object.sym=TRUE,land.pairs=pairs.matrix)
#shape.v<-b.s$symm.shape # componente simétrico da forma
#plotAllSpecimens(shape.v)
#save(shape.v,file="shape.v.RData")
load(here ("data","shape.v.RData"))

# forma média por espécie (ventral)
ventral.listed.pruned<-read.table(here ("data","Sigmodontinae.ventral.listed.pruned.txt"),h=T)
#fix(ventral.listed.pruned)
species.v<-ventral.listed.pruned[,4]
#str(species.v)

shape.v.2d<-two.d.array(shape.v)
shape.v.2d.means<-rowsum(shape.v.2d,species.v)/as.vector(table(species.v))
shape.means.v<-arrayspecs(shape.v.2d.means,dim(shape.v)[1],dim(shape.v)[2])
#shape.means.v
#str(shape.means.v)

size.means.v<-rowsum(size.v,species.v)/as.vector(table(species.v))
size.means.v
hist(log(size.means.v))
size.means.log.v<-log(size.means.v)


#################################################################################
#### Definido os Sítios ####
# Carregar Pres/Abs 228 sp 1770 sítios
presab_original <-read.table(here ("data","PresAbs_228sp_Neotropical_MainDataBase_Ordenado.txt"),h=T)
# Incluir corte de sítios com <3 sp
presab <- presab_original [rowSums (presab_original) > 3,]
## excluir entao especies que nao ocorrem nas celulas remanescentes
presab <- presab [,which(colSums(presab)>0)]

### resolver   inconsistencias no nome das spp
#colnames(presab) [which (colnames(presab) == "Abrothrix_jelskii")] <-"Abrothrix_jeslkii"
#colnames(presab) [which (colnames(presab) == "Abrothrix_longipilis")] <-"Abrothrix_longilipis"
#colnames(presab) [which (colnames(presab) == "Galenomys_garleppii")] <-"Galenomys_garleppi"
#colnames(presab) [which (colnames(presab) == "Neacomys_guianae")] <-"Neacomys_guiane"
#colnames(presab) [which (colnames(presab) == "Rhipidomys_macconnelli")] <-"Rhipidomys_maconnelli"
#colnames(presab) [which (colnames(presab) == "Thomasomys_ischyrus")] <-"Thomasomys_ischyurus"
#colnames(presab) [which (colnames(presab) == "Thomasomys_monochromos")] <-"Thomasomys_monochromus"
#colnames(presab) [which (colnames(presab) == "Thomasomys_paramorum")] <-"Thomasomys_paramarum"
#colnames(presab) [which (colnames(presab) == "Wiedomys_pyrrhorhinos")] <-"Wiedomys_pyrrhorhinus"
#colnames(presab) [which (colnames(presab) == "Zygodontomys_brevicauda")] <-"Zygodontomys_brevicaudata"
#colnames(presab) [which (colnames(presab) == "Abrothrix_andina")] <-"Abrothrix_andinus"
#colnames(presab) [which (colnames(presab) == "Abrothrix_lanosa")] <-"Abrothrix_lanosus"
##colnames(presab) [which (colnames(presab) == "Abrothrix_longilipis")] <- "Abrothrix_longipilis"
#colnames(presab) [which (colnames(presab) == "Abrothrix_olivacea")] <- "Abrothrix_olivaceus"
#colnames(presab) [which (colnames(presab) == "Galenomys_garleppii")] <- "Galenomys_garleppi"
#colnames(presab) [which (colnames(presab) == "Neomicroxus_bogotensis")] <- "Akodon_bogotensis"
#colnames(presab) [which (colnames(presab) == "Paynomys_macronyx")] <- "Chelemys_macronyx"
#colnames(presab) [which (colnames(presab) == "Phyllotis_gerbillus")] <- "Paralomys_gerbillus"

## load spatial data

longlat<-read.table(here ("data","Lon-lat-Disparity.txt"),h=T)
longlat  <- longlat [rowSums (presab_original) > 3,]

#################################################################################
### Disparidade Observada Empírica por Sítio ####
require(SYNCSA)

colnames(occ) [which(rownames (data) %in% colnames(occ) == F)]

## funcao serve para pegar os atributos observados e simulados e calcular rao nulo e sd nulo em uma filogenia

funcao_disparidade <- function (occ, traits,n_iterations) {
  data<- as.matrix(traits) # shape.v.2d.means
  colnames(data)<-seq(1,ncol(data))
  
  ## organizar os dados observados
  organized_obs_data <- organize.syncsa (comm=as.matrix(occ),traits=data)
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

RAO_OBS <- funcao_disparidade (occ = as.matrix (presab),traits= shape.v.2d.means,n_iterations = 2)

#data<-shape.v.2d.means
#colnames(data)=c(1:112)
#head(data)

## organizar os dados observados
organized_obs_data <- organize.syncsa (comm=as.matrix(presab),traits=data)

rao<-rao.diversity(organized_obs_data$community,traits=organized_obs_data$traits)
rao$FunRao
names(rao)

RESULTADOS.RAO.NULO<-matrix(NA,length(presab[,1]),1000)
for(i in 1:1000){
  row.names(data)<-sample(row.names(data))
  org.perm<-organize.syncsa(comm=presab,traits=data)
  rao.perm<-rao.diversity(org.perm$community,traits=org.perm$traits)
  RESULTADOS.RAO.NULO[,i]<-rao.perm$FunRao
}

RESULTADOS.RAO.NULO
mean.rao.NULO<-apply(RESULTADOS.RAO.NULO,1,mean)
sd.rao.NULO<-apply(RESULTADOS.RAO.NULO,1,sd)

### SES NULO ####
#cor(rao$FunRao,mean.rao.NULO)
ES_NULO<-rao$FunRao - mean.rao.NULO
SES_NULO<- ES_NULO / sd.rao.NULO

#SES_noNA<-na.omit(SES_NULO)
#colors<-c("white",jet.colors(max((SES_noNA+10)*50)))
#plot(longlat$LONG,longlat$LAT,col=colors[((SES_noNA+10)*50)],pch=15)
#maps::map("world",add=T)

#################################################################################
### Disparidade Simulada por Sítio ####
### Árvore Filogenética ####
require(ape)

#### 1 - Importar as Arvores (100 filogenias resolvidas)
tree_list <- tree <- read.nexus(file="Sigmodontinae_413species100Trees.trees")
## corrigir os nomes
tree_list <- lapply (tree_list, function (i)
	{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# estas especies nao estao no conjunto de  filogenias
# remover<-c("_Rattus_norvegicus","Tylomys_nudicaudus","Tylomys_watsoni",
#           "Ototylomys_phyllotis","Nyctomys_sumichrasti","Otonyctomys_hatti")
# tree<-drop.tip(tree,remover)
#plot(tree_list[[1]],'fan',cex=0.3)
#str(tree)
require(geiger)
tree.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$phy)
comm.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$data)[[1]]

#tree2<-tree.pruned$phy
#presab2<-tree.pruned$data

# Simulações
require(mvMORPH)

# 100 simulações aplicadas em cada uma das 100 filogenias
simul_BM<-lapply (tree.pruned, function (phy) 
	mvSIM(phy,nsim=2,model="BM1",param = list(sigma = 1)))

#mean_BM<-apply(simul_BM,1,mean)
#sd_BM<-apply(simul_BM,1,sd)
#mean_BM

# BM médio
#org<-organize.syncsa(comm=presab,traits=as.data.frame(mean_BM))
#presab2<-dplyr::select(presab,-c(as.character(org$list.warning$traits$spp[,1])))
#org<-organize.syncsa(comm=presab2,traits=as.data.frame(mean_BM))
#names(org)
#raoBM<-rao.diversity(org$community,traits=org$traits)
#names(raoBM)

# BM média e desvio
#RESULTADOS.RAO.BM.MATRIZ <- lapply (simul_BM, function (phy) 
#	matrix(NA,length(comm.pruned[1,]),dim(phy)[2]))
## cuidar que  a matriz est? transposta (em 'comm.pruned')

## inserindo informacao dos atributos das sp presentes nas comunidades
#RESULTADOS.RAO.BM.FILOGENIAS <- parLapply(cl, as.list(seq(1,length(RESULTADOS.RAO.BM.MATRIZ))), function (phy)
	#lapply(seq(1, dim(simul_BM[[phy]])[2]), function (i) {
#	for (i in 1:dim(simul_BM[[phy]])[2]) {
#	  org<-organize.syncsa(comm=t(comm.pruned),traits=as.data.frame(simul_BM [[phy]][,i]))
#	  raoBM<- rao.diversity(org$community,traits=org$traits)
#	  RESULTADOS.RAO.BM.MATRIZ[[phy]][,i]<-as.vector(raoBM$FunRao)
#	}
#   )
#)


## organizar os dados - repetir para todas as filogenias 
require(parallel)
cl <- makeCluster (6)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_BM", "comm.pruned"))
clusterEvalQ(cl,library("SYNCSA"))


organized_data <- parLapply (cl,simul_BM, function (phy)
	lapply (seq(1, dim(phy) [2]), function (simul)
		organize.syncsa(comm=t(comm.pruned),traits=as.data.frame(phy [,simul]))))

stopCluster (cl)

### calcular Rao, repetir para todas as filogenias

cl <- makeCluster (6)# numero de nucleos do pc para usar
clusterExport(cl, c("organized_data"))
clusterEvalQ(cl,library("SYNCSA"))

Rao <- parLapply (cl, seq(1,length(organized_data)), function (phy)
		lapply (seq(1, length(organized_data [[1]])), function (simul)
			rao.diversity (organized_data[[phy]][[simul]]$community, 
					traits=organized_data[[phy]][[simul]]$traits)$FunRao))

stopCluster (cl)

mean.rao.BM<- lapply(Rao, function (phy) 
	apply (do.call (rbind,phy),2,mean))

sd.rao.BM<- lapply(Rao, function (phy) 
	apply (do.call (rbind,phy),2,sd))


### Plot Rao Simulated Maps ####
#jet.colors <-
#  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
#                   space="Lab",interpolate="spline")
#colors<-c("white",jet.colors(max(mean.rao.BM*100)))
#plot(longlat$LONG,longlat$LAT,col=colors[mean.rao.BM*100],pch=15)
#maps::map("world",add=T)


#################################################################################
### SES BM ####
cor(rao$FunRao,mean.rao.BM)
ES_BM<-rao$FunRao - mean.rao.BM
SES_BM<- ES_BM / sd.rao.BM
write.csv(SES_BM,"SES_BM.csv")

#colors<-c("white",jet.colors(max((SES_BM+1)*50)))
#plot(longlat$LONG,longlat$LAT,col=colors[((SES_BM+1)*50)],pch=15)
#maps::map("world",add=T)


#################################################################################
### MPD e SES.MPD ####
require(picante)
require(geiger)
match.species<-treedata(tree,t(presab)) 
names(match.species)
#tree1<-match.species$phy
#mpd<-mpd(t(match.species$data),cophenetic(tree1))
#mpd
ses.mpd<-ses.mpd(t(match.species$data),cophenetic(tree),null.model="taxa.labels")
ses.mpd$mpd.obs.z
ses.mpd2<-ses.mpd(t(match.species$data),cophenetic(tree),null.model="independentswap")

plot(rao$FunRao~ses.mpd$mpd.obs.z)
plot(mean.rao.BM~ses.mpd$mpd.obs.z)
plot(SES_BM~ses.mpd$mpd.obs.z)


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
RESULTADOS.RAO.OU
mean.rao.OU<-apply(RESULTADOS.RAO.OU,1,mean)
sd.rao.OU<-apply(RESULTADOS.RAO.OU,1,sd)

cor(rao$FunRao,mean.rao.OU)
ES_OU<-rao$FunRao - mean.rao.OU
SES_OU<- ES_BM / sd.rao.OU

plot(rao$FunRao~ses.mpd$mpd.obs.z)
points(mean.rao.BM~ses.mpd$mpd.obs.z,col="red")
plot(SES_BM~ses.mpd$mpd.obs.z)
points(mean.rao.EB~ses.mpd$mpd.obs.z,col="blue")
plot(SES_EB~ses.mpd$mpd.obs.z)
points(mean.rao.OU~ses.mpd$mpd.obs.z,col="green")
plot(SES_OU~ses.mpd$mpd.obs.z)

plot(rao$FunRao~ses.mpd$mpd.obs)
points(mean.rao.BM~ses.mpd$mpd.obs,col="red")
plot(SES_BM~ses.mpd$mpd.obs)
points(mean.rao.EB~ses.mpd$mpd.obs,col="blue")
plot(SES_EB~ses.mpd$mpd.obs)
points(mean.rao.OU~ses.mpd$mpd.obs,col="green")
plot(SES_OU~ses.mpd$mpd.obs)



################################################################################
#### SES Simulados ####
SES_NULO # SES NULO com o Atributo Observado
SES_NULO_BM<-(mean.rao.BM - mean.rao.NULO) / sd.rao.NULO # SES NULO com o Atributo BM
SES_NULO_EB<-(mean.rao.EB - mean.rao.NULO) / sd.rao.NULO # SES NULO com o Atributo BM
SES_NULO_OU<-(mean.rao.OU - mean.rao.NULO) / sd.rao.NULO # SES NULO com o Atributo BM

plot(SES_NULO~ses.mpd$mpd.obs.z,ylim=c(-10,4))
points(SES_NULO_BM~ses.mpd$mpd.obs.z,col="blue")
points(SES_NULO_EB~ses.mpd$mpd.obs.z,col="red")
points(SES_NULO_OU~ses.mpd$mpd.obs.z,col="green")

plot(SES_NULO~ses.mpd$mpd.obs,ylim=c(-10,4))
points(SES_NULO_BM~ses.mpd$mpd.obs,col="blue")
points(SES_NULO_EB~ses.mpd$mpd.obs,col="red")
points(SES_NULO_OU~ses.mpd$mpd.obs,col="green")


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

