# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output", "mpd_results.RData"))
load (here("output", "RAO_BM.RData"))
load (here("output", "RAO_EB.RData"))
load (here("output", "RAO_OU.RData"))
load (here("output", "RAO_OBS.RData"))

# relationship between empirical and simulated data sets

lapply(seq(1,100), function (i) {
plot(NA, ylim=c(-7,7),xlim=c(-7,4),
     xlab=expression (paste("Mean Pairwise Phylogenetic Distance", " (MPD"['SES'],")")),
     ylab=expression (paste("Morphological Disparity", "  (Disparity"['SES'],")")))
points(null.mpdf[[i]]$SES.MPD, RAO_OBS$SES, 
     col=rgb(red=0.3,green=0.3,blue=0.3, alpha=0.1),
     pch=19)
abline(lm(RAO_OBS$SES ~null.mpdf[[i]]$SES.MPD),
       col=rgb(red=0.3,green=0.3,blue=0.3, alpha=0.9),
       lwd=2)
points(null.mpdf[[i]]$SES.MPD, RAO_BM[[i]]$SES,
       col=rgb(red=0,green=0.3,blue=0.9, alpha=0.1),
       pch=19)
abline(lm(RAO_BM[[i]]$SES ~null.mpdf[[i]]$SES.MPD),
       col=rgb(red=0,green=0.3,blue=0.9, alpha=0.9),
       lwd=2)
points(null.mpdf[[i]]$SES.MPD, RAO_EB[[i]]$SES,
       col=rgb(red=1,green=0,blue=0, alpha=0.1),
       pch=19)
abline(lm(RAO_EB[[i]]$SES ~null.mpdf[[i]]$SES.MPD),
       col=rgb(red=1,green=0,blue=0, alpha=0.9),
       lwd=2)
points(null.mpdf[[i]]$SES.MPD, RAO_OU[[i]]$SES,
       col=rgb(red=0.1,green=0.9,blue=0.1, alpha=0.1),
       pch=19)
abline(lm(RAO_OU[[i]]$SES ~null.mpdf[[i]]$SES.MPD),
       col=rgb(red=0.1,green=0.9,blue=0.1, alpha=0.9),
       lwd=2)

legend ("bottomright", 
        legend = c("Empirical", "BM","EB","OU"),
        bty="n",
        pch=19,
        col = c(rgb(red=0.3,green=0.3,blue=0.3, alpha=0.7),
                rgb(red=0,green=0.3,blue=0.9, alpha=0.9),
                rgb(red=1,green=0,blue=0, alpha=0.9),
                rgb(red=0.1,green=0.9,blue=0.1, alpha=0.9)
                
        ))


})

# working with the average

avBM<-do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T))
avEB<-do.call(cbind,sapply(RAO_EB, "[","SES",simplify=T))
avOU<-do.call(cbind,sapply(RAO_OU, "[","SES",simplify=T))
avMPD<-do.call(cbind,sapply(null.mpdf, "[","SES.MPD",simplify=T))

plot(NA, ylim=c(-6,2),xlim=c(-6.5,2.5),
     xlab=expression (paste("Mean Pairwise Phylogenetic Distance", " (MPD"['SES'],")")),
     ylab=expression (paste("Morphological Disparity", "  (Disparity"['SES'],")")))
points(rowMeans (avMPD), RAO_OBS$SES, 
       col=rgb(red=0.3,green=0.3,blue=0.3, alpha=0.1),
       pch=19)
points(rowMeans (avMPD), rowMeans(avBM),
       col=rgb(red=0,green=0.3,blue=0.9, alpha=0.1),
       pch=19)
points(rowMeans (avMPD), rowMeans(avEB),
       col=rgb(red=1,green=0,blue=0, alpha=0.1),
       pch=19)
points(rowMeans (avMPD), rowMeans(avOU),
       col=rgb(red=0.1,green=0.9,blue=0.1, alpha=0.1),
       pch=19)

legend ("bottomright", 
        legend = c("Empirical", "BM","EB","OU"),
        bty="n",
        pch=19,
        col = c(rgb(red=0.3,green=0.3,blue=0.3, alpha=0.7),
                rgb(red=0,green=0.3,blue=0.9, alpha=0.9),
                rgb(red=1,green=0,blue=0, alpha=0.9),
                rgb(red=0.1,green=0.9,blue=0.1, alpha=0.9)
                
        ))


# test of whether slope of the regression between ses disparity and ses MPD
# varies across evolutionary models

# bind data to analysis

bind_data_disparity <- lapply (seq (1,length(null.mpdf)), function (i)
  
  rbind (
    data.frame (
      dataset="AOBS",
      disp=RAO_OBS$SES,
      MPD=null.mpdf[[i]]$SES.MPD),
  
    data.frame (
      dataset="BM",
      disp=RAO_BM[[i]]$SES,
      MPD=null.mpdf[[i]]$SES.MPD),
  
    data.frame (
      dataset="EB",
      disp=RAO_EB[[i]]$SES,
      MPD=null.mpdf[[i]]$SES.MPD),
  
    data.frame (
      dataset="OU",
      disp=RAO_OU[[i]]$SES,
      MPD=null.mpdf[[i]]$SES.MPD)
  )
)


# One model per phylogeny
model_slope <- lapply (bind_data_disparity, function (i) {
  
  # using GLM
  m1<-lm (disp ~ MPD*dataset,
           data=i)
  
  # extracting estimates
  # intercept (average disparity in the empirical dataset (level "AOBS"))
  res <- data.frame(EMPIRICAL=m1$coefficients[1], # EMpirical
             BM=m1$coefficients[6], # BM
             EB=m1$coefficients[7], # EB
            OU=m1$coefficients[8]#,
            #R2 = RsquareAdj(m1)$adj.r.squared
            ) # OU
  ; # return
  res
  
})

# transforming res list into df
model_slope_df <- do.call (rbind, model_slope)

# finally melt
require(reshape)
require(ggplot2)
model_slope_df <- melt(model_slope_df)

# plotting

a <- ggplot(model_slope_df, aes(x=value, color=variable, fill=variable)) +
  geom_density(size=1,alpha=0.1) +
  scale_fill_manual(values=c("EMPIRICAL"="gray", 
                             "BM"="blue",
                             "EB"="red",
                             "OU"="green")) + 
  scale_colour_manual(values=c("EMPIRICAL"="gray", 
                                        "BM"="blue",
                                        "EB"="red",
                                        "OU"="green"))+
  theme_classic()  

b <- a+ theme (
               axis.text.x = element_text(angle = 90),
               legend.title = element_blank(),
               #legend.position = "none",
               legend.position = c(.95, .95),
               legend.justification = c("right", "top"),
               legend.box.just = "right",
               legend.margin = margin(6, 6, 6, 6),
               plot.margin=unit(c(.2,1,.1,1),"cm")) +
  xlab ("Slope SES Disparity ~ SES MPD") + 
  ylab("Density") + 
  scale_x_continuous(breaks = seq(-2,2,0.15),
                     limits = c(-1.5,1.5)) 

# exploring phylogenetic uncertainty

b








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
