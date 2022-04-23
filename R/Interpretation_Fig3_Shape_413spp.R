# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output", "mpd_results_ALL.RData"))
load (here("output", "RAO_OBS_ALL_multivariate.RData"))
load (here("output", "RAO_BM_ALL_multivariate.RData"))

# relationship between empirical and simulated data sets

# working with the average

avBM<-rowMeans (do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T)))
avMPD<-rowMeans(do.call(cbind,sapply(null.mpdf, "[","SES.MPD",simplify=T)))

# -------------------------------------------------------------------- #
# Mapping deviations of empirical disparity vs. simulated disparity 
# across models#
#
## -----------------------------------------------#
##           Load spatial data
##       Incidence of 228 species in 1770 cells

presab_original <-read.table(here ("data","PresAbs_228sp_Neotropical_MainDataBase_Ordenado.txt"),h=T)
presab_original <- presab_original [order(as.numeric(rownames(presab_original))),]
# remove sites with 2 or less sp
presab <- presab_original [rowSums (presab_original) > 3,]
## exclude species absent in the remaining cells
presab <- presab [,which(colSums(presab)>0)]
presab <- presab[order(as.numeric(rownames(presab)),decreasing=F),]

# ----------------------------- # 
##  Geographic coordinates data

longlat<-read.table(here ("data","Lon-lat-Disparity.txt"),h=T)
# rm sites with efw spp
longlat  <- longlat [,1:2]
longlat <- longlat[order(as.numeric((longlat$LONG)),decreasing=F),]
longlat <- longlat[order(as.numeric((longlat$LAT)),decreasing=F),]
longlat <- longlat [rowSums (presab_original) > 3,]

table(rownames(presab) == rownames(longlat))

# now we need to generate the neutral SES
# get the average and sd of disparity under BM simulations
mean_BM <- do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T))
mean_BMm <- apply (mean_BM, 1, mean)
mean_sd <- apply (mean_BM, 1, sd)
#mean_BMev <- mean(mean_BM)
#sd_BMev <- sd(mean_BM)

# calculate neutral SES
SES_NEUTRAL <- (RAO_OBS$Observado - mean_BMm)/mean_sd

# data to MAP

data_to_map <- data.frame(longlat,
                          'SES-EMPIRICAL'= RAO_OBS$SES,
                          #'SES-BM'=rowMeans (avBM),
                          'SES-NEUTRAL' = SES_NEUTRAL)

#write.csv(data_to_map,here("output","data_to_SAM.csv"))

melt_data_to_map <- melt(data_to_map,id=c("LONG","LAT"))
colnames(melt_data_to_map)[which(colnames(melt_data_to_map) == "value")] <- 'SES'
# round values
melt_data_to_map$SES<-round(melt_data_to_map$SES,3)

#plot using ggplot
# help overlap things
# https://datacarpentry.org/r-raster-vector-geospatial/03-raster-reproject-in-r/
# https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/

rel <- raster(here ("data", "DEM_Elev.tif"))
rel<-projectRaster(rel, crs="+proj=longlat +datum=WGS84") # reproject raster to fit community data
rel_spdf <- as(rel, "SpatialPixelsDataFrame")
rel <- as.data.frame(rel_spdf)

# EMPIRICAL VS NULL MODEL
map1 <- ggplot(melt_data_to_map[which(melt_data_to_map$variable == 'SES.EMPIRICAL'),], 
               aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = SES),alpha=0.65)+
  #facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(range(melt_data_to_map$SES)[1],
                                range(melt_data_to_map$SES)[2]),
                       breaks=seq(range(melt_data_to_map$SES)[1],
                                  range(melt_data_to_map$SES)[2],
                                  1.1),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  ggtitle ("B")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size=7,angle=75,hjust = 1),
        legend.title = element_text(size=8),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        plot.margin=unit(c(0,2,2,2),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))


## ---------------------------------
common_legend_SES <- get_legend(map1) ## get the legend of a map with legend

# map without legend
map1 <- ggplot(melt_data_to_map[which(melt_data_to_map$variable == 'SES.EMPIRICAL'),], 
               aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = SES),alpha=0.65)+
  #facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(range(melt_data_to_map$SES)[1],
                                range(melt_data_to_map$SES)[2]),
                       breaks=seq(range(melt_data_to_map$SES)[1],
                                  range(melt_data_to_map$SES)[2],
                                  1.6),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  ggtitle ("B")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="none",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))


# empirical vs neutral
map2 <- ggplot(melt_data_to_map[which(melt_data_to_map$variable == 'SES.NEUTRAL'),], 
               aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = SES),alpha=0.65)+
  #facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(range(melt_data_to_map$SES)[1],
                                range(melt_data_to_map$SES)[2]),
                       breaks=seq(range(melt_data_to_map$SES)[1],
                                  range(melt_data_to_map$SES)[2],
                                  1.6),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  ggtitle ("C")+
  theme(legend.position="none",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))

# finally, mapping observed (empirical disparity). Separate function due to diff scale
data_to_map_emp <- data.frame(longlat,
                              #Richness=rowSums(presab),
                              'DISPARITY' = RAO_OBS$Observado)

melt_data_to_map_emp <- melt(data_to_map_emp,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp)[which(colnames(melt_data_to_map_emp) == "value")] <- 'RAO'

## empirical vs simulated
#plot using ggplot

map3 <- ggplot(melt_data_to_map_emp, 
               aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = RAO),alpha=0.65)+
  scale_fill_gradient2(midpoint = 0.35, 
                       limits=c(0,0.5),
                       breaks=seq(0,0.5,0.1),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  ggtitle ("A")+
  theme(legend.position="bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))

# ----------------------------- #
# significance

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    longlat,
    Lower.NULL=RAO_OBS$SES<=-1.96,
    #Lower.OU=((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU))<=-1.96,
    Lower.BM=SES_NEUTRAL <=-1.96
    #Lower.EB=((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB)<=-1.96)
  ), 
  # higher than
  higher = cbind (
    longlat,
    Higher.NULL=RAO_OBS$SES>=1.96,
    #Higher.OU=((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU))>=1.96,
    Higher.BM=SES_NEUTRAL>=1.96
    #Higher.EB=((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB))>=1.96
  )
)

# melt
# lower than
melt_data_cells_lower <- melt(count_cells_relative_to_models$lower,id=c("LONG","LAT"))
colnames(melt_data_cells_lower)[which(colnames(melt_data_cells_lower) == "value")] <- 'Disparity'
melt_data_cells_lower$Disparity <- ifelse (melt_data_cells_lower$Disparity == TRUE, -1, 0)

# higher than
melt_data_cells_higher <- melt(count_cells_relative_to_models$higher,id=c("LONG","LAT"))
colnames(melt_data_cells_higher)[which(colnames(melt_data_cells_higher) == "value")] <- 'Disparity'
melt_data_cells_higher$Disparity <- ifelse (melt_data_cells_higher$Disparity == TRUE, 1, 0)

# reunite these data
melt_data_cells_higher$Disparity <- as.factor(melt_data_cells_higher$Disparity + melt_data_cells_lower$Disparity)
melt_data_cells_higher$variable <- gsub ("Higher.","",melt_data_cells_higher$variable)
melt_data_cells_higher$variable<-factor(melt_data_cells_higher$variable,
                                        levels=c("NULL", "BM"))
# and map
#plot using ggplot

panel3_NULL <- ggplot(melt_data_cells_higher[which(melt_data_cells_higher$variable=="NULL"),], 
                 aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = Disparity),alpha=0.65)+
  #facet_wrap(~variable,scales = "fixed",ncol=4)+
  scale_fill_manual(
    values = c("-1" ="#268bd2","0"= "#eee8d5","1" ="#dc322f"),
    labels = c("SES<=-1.96","-1.96>SES<1.96", "SES>=1.96")
  ) + 
  theme_classic() + 
  theme_map() +
  ggtitle ("D")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="none",
        legend.justification = "center",
        legend.direction = "vertical",
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        plot.background = element_rect(fill="white",colour="white"),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))
# BM

panel3_BM <- ggplot(melt_data_cells_higher[which(melt_data_cells_higher$variable=="BM"),], 
                      aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = Disparity),alpha=0.65)+
  
  #facet_wrap(~variable,scales = "fixed",ncol=4)+
  scale_fill_manual(
    values = c("-1" ="#268bd2","0"= "#eee8d5","1" ="#dc322f"),
    labels = c("SES<=-1.96","-1.96>SES<1.96", "SES>=1.96")
  ) + 
  theme_classic() + 
  theme_map() +
  ggtitle ("D")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        plot.background = element_rect(fill="white",colour="white"),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))

# extract legend

common_legend_significance <- get_legend(panel3_BM) ## get the legend of a map with legend

# plot without legend
panel3_BM <- ggplot(melt_data_cells_higher[which(melt_data_cells_higher$variable=="BM"),], 
                    aes(x = LONG, y = LAT)) +
  geom_raster(data = rel, aes_string(x = "x", y = "y", 
                                     alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = Disparity),alpha=0.65)+
  
  #facet_wrap(~variable,scales = "fixed",ncol=4)+
  scale_fill_manual(
    values = c("-1" ="#268bd2","0"= "#eee8d5","1" ="#dc322f"),
    labels = c("SES<=-1.96","-1.96>SES<1.96", "SES>=1.96")
  ) + 
  theme_classic() + 
  theme_map() +
  ggtitle ("E")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="none",
        plot.background = element_rect(fill="white",colour="white"),
        plot.title=element_text(size=10,hjust = 0.5),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size=10))



#ggsave (here("output","maps.pdf"),height = 6,width=5)

## arrange into a grid
panel2 <- grid.arrange(map3,
                       map1,
                       map2,
                       common_legend_SES,
                       panel3_NULL,
                       panel3_BM,
                       common_legend_significance,
                       ncol=4,nrow = 20,
                       layout_matrix = rbind (c(NA,1,1,NA),
                                              c(NA,1,1,NA),
                                              c(NA,1,1,NA),
                                              c(NA,1,1,NA),
                                              c(2,2,3,3),
                                              c(2,2,3,3),
                                              c(2,2,3,3),
                                              c(NA,4,4,NA),
                                              c(5,5,6,6),
                                              c(5,5,6,6),
                                              c(5,5,6,6),
                                              c(NA,7,7,NA)))
# white background
panel2 <- cowplot::ggdraw(panel2) + 
  theme(plot.background = element_rect(fill="white", color = NA))

pdf(here ("output", "vectorized","Fig6_relief.pdf"),width=4,height =9)
panel2
dev.off()
## ------------------------- ##

## comparison of nulls/simulations

obsBM<-do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T))
obsBM <- apply (obsBM, 1, mean)

##
data_to_map_emp_sim <- data.frame(longlat,
                                  #'Observed Disparity' = RAO_OBS$Observado,
                                  'Null Disparity' = RAO_OBS$med_nulo,
                                  'BM Disparity' = obsBM
                                  )

# melt
melt_data_to_map_emp_sim <- melt(data_to_map_emp_sim,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp_sim)[which(colnames(melt_data_to_map_emp_sim) == "value")] <- 'Disparity'

# and map
#plot using ggplot

alternative_map1 <- ggplot(melt_data_to_map_emp_sim, 
                           aes(x = LONG, y = LAT)) +
  #geom_raster(data = rel, aes_string(x = "x", y = "y", 
  #                                   alpha = "DEM_Elev"),show.legend=T) +
  geom_tile(aes(fill = Disparity),alpha=0.7)+
  facet_wrap(~variable,scales = "fixed",ncol=2)+
  scale_fill_gradient2(midpoint = 0.3,
                       limits=c(0,0.5),
                       breaks=seq(0,0.5,0.1),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="right",
        legend.justification = "center",
        legend.direction = "vertical",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        plot.background = element_rect(fill="white",colour="white"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=15),
        strip.background = element_blank())

#pdf(here ("output", "vectorized","Fig8_relief.pdf"),width=5,height =2.5)
#alternative_map1
#dev.off()

## correlation between average null and simulated by OU

cor (data.frame (RAO_OBS$med_nulo, 
                 obsBM))

# save data to produce maps in QGIS

#colnames(melt_data_to_map_emp) <- c("LONG", "LAT", "variable", "values")
#colnames(melt_data_to_map) <- c("LONG", "LAT", "variable", "values")
#colnames(melt_data_cells_higher) <- c("LONG", "LAT", "variable", "values")
## bind data
#df_qgis <- rbind (melt_data_to_map_emp,
#                  melt_data_to_map,
#                  melt_data_cells_higher)
#write.csv (df_qgis, file = here ("output", "df_qgis_shape.csv"))
#