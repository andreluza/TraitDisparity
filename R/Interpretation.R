# ------------------------# 

# load packages
source ("./R/packages.R")

# relationship between empirical and simulated data sets



plot(rao$FunRao~ses.mpd$mpd.obs.z)
plot(mean.rao.BM~ses.mpd$mpd.obs.z)
plot(SES_BM~ses.mpd$mpd.obs.z)


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


# 
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



# test of whether slope of the regression between ses disparity and ses MPD
# varies across evolutionary models








# exploring phylogenetic uncertainty







# mapping deviations of empirical disparity vs. simulated disparity 
# across models







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
#longlat  <- longlat [rowSums (presab_original) > 3,]

library(sp)
library(raster)
library (rgdal)

coord_data <- longlat[,c("LONG","LAT")]

coordinates (coord_data)<- ~ LONG+LAT

# load shapafile of neotropics
neotropics <- readOGR(dsn=here ("data","Lowenberg_Neto_2014_shapefile"),layer="Lowenberg_Neto_2014")
neotropics <- spTransform (neotropics,CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# objects in the same proj
coord_data <- spTransform(coord_data,crs(neotropics))
# make grid over neotropics
grdpts <- makegrid(neotropics, cell=1)
# turn into ponts
spgrd <- SpatialPoints(grdpts, proj4string = CRS(proj4string(neotropics)))
# into pixels
spgrdWithin <- SpatialPixels(spgrd[neotropics,])
# and transform into grids
spgrdWithin <- as(spgrdWithin, "SpatialGrid")

# transforming into raster
r<-(raster(spgrdWithin))
r[]<-1:ncell(r)

# plot
plot(r)
plot(neotropics,add=T)

plot(spgrdWithin, add = T)
points(coord_data,pch=19,cex=0.1,col="red")

# mask based on coord values
subset_neo <- mask (r,coord_data )
subset_neo[] <- longlat$SES.SIZE.DISPARITY
plot(subset_neo)

points(coord_data,pch=19,cex=0.1,col="red")
