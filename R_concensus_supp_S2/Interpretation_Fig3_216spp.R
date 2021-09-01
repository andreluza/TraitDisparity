# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_concensus_supp_S2", "mpd_results.RData"))
load (here("output_concensus_supp_S2", "RAO_BM.RData"))
load (here("output_concensus_supp_S2", "RAO_EB.RData"))
load (here("output_concensus_supp_S2", "RAO_OU.RData"))
load (here("output_concensus_supp_S2", "RAO_OBS.RData"))


# relationship between empirical and simulated data sets

# working with the average

avBM<-RAO_BM$SES
avEB<-RAO_EB$SES
avOU<- RAO_OU$SES
avMPD<- statistics.phy$SES.MPD

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
mean_BM <- mean(RAO_BM$Observado)
sd_BM <- sd(RAO_BM$Observado)

# calculate neutral SES
SES_NEUTRAL <- (RAO_OBS$Observado - mean_BM)/sd_BM

# data to MAP

data_to_map <- data.frame(longlat,
                          #Richness=rowSums(presab),
                          'SES-EMPIRICAL'= RAO_OBS$SES,
                          'SES-BM'= (avBM),
                          'SES-EB' = (avEB),
                          'SES-OU'= (avOU),
                          'SES-NEUTRAL' = SES_NEUTRAL)

#write.csv(data_to_map,here("output","data_to_SAM.csv"))

melt_data_to_map <- melt(data_to_map,id=c("LONG","LAT"))
colnames(melt_data_to_map)[which(colnames(melt_data_to_map) == "value")] <- 'SES'

## empirical vs simulated
#plot using ggplot
# EMPIRICAL VS NULL MODEL
map1 <- ggplot(melt_data_to_map[which(melt_data_to_map$variable == 'SES.EMPIRICAL'),], 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = SES)) +
  #facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(-4.1,2),
                       breaks=seq(-4.2,2,0.9),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  ggtitle ("B) SES - Null Disparity")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="right",
        legend.justification = "center",
        legend.direction = "vertical",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        plot.title=element_text(size=15,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=15))

# empirical vs neutral
map2 <- ggplot(melt_data_to_map[which(melt_data_to_map$variable == 'SES.NEUTRAL'),], 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = SES)) +
  #facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(-4.1,2),
                       breaks=seq(-4.2,2,0.9),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  ggtitle ("C) SES - Neutral Disparity")+
  theme(legend.position="none",
        legend.justification = "center",
        legend.direction = "vertical",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        plot.title=element_text(size=15,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=15))


#ggsave (here("output","maps.pdf"),height = 6,width=5)

# finally, mapping observed (empirical disparity). Separate function due to diff scale
data_to_map_emp <- data.frame(longlat,
                              #Richness=rowSums(presab),
                              'DISPARITY' = RAO_OBS$Observado)

melt_data_to_map_emp <- melt(data_to_map_emp,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp)[which(colnames(melt_data_to_map_emp) == "value")] <- 'RAO.ENTROPY'

## empirical vs simulated
#plot using ggplot

map3 <- ggplot(melt_data_to_map_emp, 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = RAO.ENTROPY)) +
  scale_fill_gradient2(midpoint = 0.35, 
                       limits=c(0,0.5),
                       breaks=seq(0,0.5,0.15),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  ggtitle ("A) MORPHOLOGICAL DISPARITY")+
  theme(legend.position="right",
        legend.justification = "center",
        legend.direction = "vertical",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        plot.title=element_text(size=15,hjust = 0.5),
        plot.background = element_rect(fill="white",colour="white"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=15))


#ggsave (here("output","maps.pdf"),height = 6,width=5)

## arrange into a grid
panel2 <- grid.arrange(map3,
                       map1,
                       map2,
                       ncol=8,nrow = 8,
                       layout_matrix = rbind (c(NA,1,1,1,1,1,NA),
                                              c(NA,1,1,1,1,1,NA),
                                              c(NA,1,1,1,1,1,NA),
                                              c(NA,1,1,1,1,1,NA),
                                              c(NA,1,1,1,1,1,NA),
                                              c(2,2,2,2,3,3,3),
                                              c(2,2,2,2,3,3,3),
                                              c(2,2,2,2,3,3,3),
                                              c(2,2,2,2,3,3,3)
                       ))

panel2

## ------------------------- ##
## alternative panel
## show simulated values

obsBM<- RAO_BM$Observado
obsEB<- RAO_EB$Observado
obsOU <- RAO_OU$Observado

##
data_to_map_emp_sim <- data.frame(longlat,
                                  #'Observed' = RAO_OBS$Observado,
                                  'Null Disparity' = RAO_OBS$med_nulo,
                                  'BM Disparity' = obsBM,
                                  'EB Disparity' = obsEB,
                                  'OU Disparity' = obsOU)

# melt
melt_data_to_map_emp_sim <- melt(data_to_map_emp_sim,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp_sim)[which(colnames(melt_data_to_map_emp_sim) == "value")] <- 'Disparity'

# and map
#plot using ggplot

alternative_map1 <- ggplot(melt_data_to_map_emp_sim, 
                           aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = Disparity)) +
  facet_wrap(~variable,scales = "fixed",ncol=2)+
  scale_fill_gradient2(midpoint = 0.3,
                       limits=c(0,0.6),
                       breaks=seq(0,0.6,0.1),
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
        strip.text = element_text(size=15))

alternative_map1

## correlation between average null and simulated by OU

cor (data.frame (RAO_OBS$med_nulo, 
                 obsBM,
                 obsEB,
                 obsOU))

# panel 3
# plotting the cells with observed disparity either higher or lower than the disparity found for 3 evol models

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    longlat,
    Lower.NULL=RAO_OBS$SES<=-1.96,
    Lower.OU=((RAO_OBS$Observado - mean(obsOU))/sd(obsOU))<=-1.96,
    Lower.BM=((RAO_OBS$Observado - mean(obsBM))/sd(obsBM))<=-1.96,
    Lower.EB=((RAO_OBS$Observado - mean(obsEB))/sd(obsEB)<=-1.96)
  ), 
  # higher than
  higher = cbind (
    longlat,
    Higher.NULL=RAO_OBS$SES>=1.96,
    Higher.OU=((RAO_OBS$Observado - mean(obsOU))/sd(obsOU))>=1.96,
    Higher.BM=((RAO_OBS$Observado - mean(obsBM))/sd(obsBM))>=1.96,
    Higher.EB=((RAO_OBS$Observado - mean(obsEB))/sd(obsEB))>=1.96
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

# unite these data
melt_data_cells_higher$Disparity <- as.factor(melt_data_cells_higher$Disparity + melt_data_cells_lower$Disparity)
melt_data_cells_higher$variable <- gsub ("Higher.","",melt_data_cells_higher$variable)
melt_data_cells_higher$variable<-factor(melt_data_cells_higher$variable,
                                        levels=c("NULL", "BM","EB","OU"))
# and map
#plot using ggplot

panel3 <- ggplot(melt_data_cells_higher, 
                 aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = Disparity)) +
  facet_wrap(~variable,scales = "fixed",ncol=4)+
  scale_fill_manual(
    values = c("-1" ="#268bd2","0"= "#eee8d5","1" ="#dc322f"),
    labels = c("SES<=-1.96","-1.96>SES<1.96", "SES>=1.96")
  ) + 
  theme_classic() + 
  theme_map() +
  ggtitle ("Empirical relative to simulated disparity")+
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        plot.background = element_rect(fill="white",colour="white"),
        plot.title=element_text(size=15,hjust = 0.5),
        plot.margin=unit(c(2,.1,2,.1),"cm"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=15))

panel3

#

