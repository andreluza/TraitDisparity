# ------------------------# 

# load packages
source ("./R/packages.R")
source ("./R/functions.R")

#----------------------------------------#
#       Carregar Forma do Crânio 
#             Shape ventral

## Carregar arquivo .tps com vista ventral
tps.ventral<-readland.tps(here ("data","Sigmodontinae.ventral.dig.tps"),specID = "ID", readcurves = FALSE)
##dim(tps.ventral)

## Carregar lista (classificadores) vista ventral
ventral.listed<-read.table(here ("data","Sigmodontinae.ventral.listed.txt"),h=T)
##fix(ventral.listed)
species.v<-ventral.listed[,4]
# unique (species.v)

## Eliminar Nephelomys pirrensis & Sigmodon zanjonensis (sem shapefile) & Brucepattersonius_sp 
## Eliminar 9 (fora do recorte espacial): Rheomys underwoodi, Rheomys mexicanus, Oryzomys palustris, Nectomys p. palmipes, Nesoryzomys swarthi, Nesoryzomys indefessus, Nesoryzomys fernandinae, Nesoryzomys darwini, Aegialomys galapagoensis
## Ending with 228 species
rm_spp <- c("Nephelomys_pirrensis","Sigmodon_zanjonensis","Aegialomys_galapagoensis","Brucepattersonius_sp","Rheomys_underwoodi","Nesoryzomys_darwini","Rheomys_mexicanus","Oryzomys_palustris","Nesoryzomys_fernandinae","Nectomys_palmipes","Nesoryzomys_swarthi","Nesoryzomys_indefessus")
tps.ventral<-tps.ventral[,,which(ventral.listed$Species_Patton2015 %in% rm_spp == F)] 
#dim(tps.ventral)
#
## GPA Ventral
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
size.means.v
hist(log(size.means.v))
size.means.log.v<-log(size.means.v)

# -----------------------------------------------#
#           Load spatial data
#       Incidence of 228 species in 1770 cells

presab_original <-read.table(here ("data","PresAbs_228sp_Neotropical_MainDataBase_Ordenado.txt"),h=T)
# Incluir corte de sítios com <3 sp
presab <- presab_original [rowSums (presab_original) > 3,]
## excluir entao especies que nao ocorrem nas celulas remanescentes
presab <- presab [,which(colSums(presab)>0)]

### resolver   inconsistencias no nome das spp - already solved

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

# ----------------------------- # 
##  Geographic coordinates data

longlat<-read.table(here ("data","Lon-lat-Disparity.txt"),h=T)
longlat  <- longlat [rowSums (presab_original) > 3,]

# --------------------------------------------------------- #
# Calculando a disparidade Observada (Empírica) por célula 

# numero de iteracoes (tb serve para as simulacoes)
niter <- 999

# rodar a funcao de disparidade observada
RAO_OBS <- funcao_disparidade (occ = as.matrix (presab),
                               traits= shape.v.2d.means,
                               n_iterations = niter)

save (RAO_OBS, file = here("output","RAO_OBS.RData"))
gc()
# --------------------------------------------------- #
#          Simular disparidade por celula, 
#     de acordo com diferentes modelos evolutivos

# Carregar filogenia
# 1 - Importar as Arvores (100 filogenias resolvidas)
tree_list <- tree <- read.nexus(file=here("data","Sigmodontinae_413species100Trees.trees"))
## corrigir os nomes
tree_list <- lapply (tree_list, function (i)
	{i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);i})

# estas especies nao estao no conjunto de  filogenias
# remover<-c("_Rattus_norvegicus","Tylomys_nudicaudus","Tylomys_watsoni",
#           "Ototylomys_phyllotis","Nyctomys_sumichrasti","Otonyctomys_hatti")
# tree<-drop.tip(tree,remover)
#plot(tree_list[[1]],'fan',cex=0.3)
#str(tree)

# matching entre observacoes, traits, e filogenia
tree.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$phy)
comm.pruned<-lapply (tree_list, function (phy) 
	treedata(phy,t(presab))$data)[[1]]

#tree2<-tree.pruned$phy
#presab2<-tree.pruned$data

# Simulações do modelo browniano
# 100 simulações aplicadas em cada uma das 100 filogenias
simul_BM<-lapply (tree.pruned, function (phy) 
	mvSIM(phy,nsim=niter,model="BM1",param = list(sigma = 1)))

## obter statistics

simul_BM_mean <- lapply(simul_BM, function (phy)
  apply(phy,1,mean))

## organizar os dados  de traits do mesmo modo que as analises empiricas,
# mas desta vez com dados de atributos simulados

# definir ncores para analises paralelas
ncores <- 6

# criar um cluster e carregar funcoes e dados
cl <- makeCluster (ncores)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_BM_mean", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_BM <- parLapply (cl,simul_BM_mean, function (phy)	
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

save(RAO_BM, file=here("output","RAO_BM.Rdata"))

# 100 simulações DE EB aplicadas em cada uma das 100 filogenias
simul_EB<-lapply (tree.pruned, function (phy) 
  mvSIM(phy,nsim=niter,model="EB",param = list(sigma = 1, beta = -0.5)))

simul_EB_mean <- lapply(simul_EB, function (phy)
  apply(phy,1,mean))

# run disparity function
cl <- makeCluster (ncores)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_EB_mean", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_EB <- parLapply (cl, simul_EB_mean, function (phy)	
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)

# 100 simulações DE OU aplicadas em cada uma das 100 filogenias
simul_OU <-lapply (tree.pruned, function (phy) 
  mvSIM(phy,nsim=niter,model="OU1",param = list(sigma = 1,alpha=1)))

simul_OU_mean <- lapply(simul_OU, function (phy)
  apply(phy,1,mean))

# run disparity function
cl <- makeCluster (ncores)# numero de nucleos do pc para usar
clusterExport(cl, c("simul_OU_mean", "comm.pruned","niter","funcao_disparidade"))
clusterEvalQ(cl,library("SYNCSA"))

RAO_OU <- parLapply (cl,simul_OU_mean, function (phy)	
  funcao_disparidade (occ = t(comm.pruned),
                      traits= as.data.frame(phy),
                      n_iterations = niter))

stopCluster (cl)


#################################################################################
### MPD e SES.MPD ####

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

