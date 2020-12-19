# ------------------------# 

# load packages
source ("./R/packages.R")

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
