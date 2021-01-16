# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output", "mpd_results_ALL.RData"))
load (here("output", "RAO_BM_ALL.RData"))
load (here("output", "RAO_EB_ALL.RData"))
load (here("output", "RAO_OU_ALL.RData"))
load (here("output", "RAO_OBS_ALL.RData"))

# relationship between empirical and simulated data sets

# working with the average

avBM<-do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T))
avEB<-do.call(cbind,sapply(RAO_EB, "[","SES",simplify=T))
avOU<-do.call(cbind,sapply(RAO_OU, "[","SES",simplify=T))
avMPD<-do.call(cbind,sapply(null.mpdf, "[","SES.MPD",simplify=T))

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
mean_BM <- do.call(cbind,sapply(RAO_BM, "[","med_nulo",simplify=T))
mean_BM <- apply (mean_BM, 1, mean)
sd_BM <- do.call(cbind,sapply(RAO_BM, "[","sd_nulo",simplify=T))
sd_BM <- apply (sd_BM, 1, mean)

# calculate neutral SES
SES_NEUTRAL <- (RAO_OBS$Observado - mean_BM)/sd_BM

# data to MAP

data_to_map <- data.frame(longlat,
                          #Richness=rowSums(presab),
                          'SES-EMPIRICAL'= RAO_OBS$SES,
                          'SES-BM'=rowMeans (avBM),
                          'SES-EB' =rowMeans (avEB),
                          'SES-OU'=rowMeans (avOU),
                          'SES-NEUTRAL' = SES_NEUTRAL)

#write.csv(data_to_map,here("output","data_to_SAM.csv"))

melt_data_to_map <- melt(data_to_map,id=c("LONG","LAT"))
colnames(melt_data_to_map)[which(colnames(melt_data_to_map) == "value")] <- 'SES'

## empirical vs simulated
#plot using ggplot
map1 <- ggplot(melt_data_to_map, 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = SES)) +
  facet_wrap(~variable,scales = "fixed",ncol=3)+
  scale_fill_gradient2(midpoint = 0,
                       limits=c(-4.1,1.6),
                       breaks=seq(-4.2,1.6,0.9),
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
map1

#ggsave (here("output","maps.pdf"),height = 6,width=5)

# finally, mapping observed (empirical disparity). Separate function due to diff scale
data_to_map_emp <- data.frame(longlat,
                          #Richness=rowSums(presab),
                          'DISPARITY' = RAO_OBS$Observado)

melt_data_to_map_emp <- melt(data_to_map_emp,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp)[which(colnames(melt_data_to_map_emp) == "value")] <- 'RAO.ENTROPY'

## empirical vs simulated
#plot using ggplot

map2 <- ggplot(melt_data_to_map_emp, 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = RAO.ENTROPY)) +
  scale_fill_gradient2(midpoint = 0.35, 
                       limits=c(0,0.5),
                       breaks=seq(0,0.5,0.15),
                       mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  ggtitle ("MORPHOLOGICAL DISPARITY")+
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
panel2 <- grid.arrange(map2,
                      map1,
                      ncol=11,nrow = 3,
                      layout_matrix = rbind (c(1,1,1,1,1,2,2,2,2,2,2),
                                             c(1,1,1,1,1,2,2,2,2,2,2),
                                             c(1,1,1,1,1,2,2,2,2,2,2)))


## ------------------------- ##
## alternative panel
## show simulated values

obsBM<-do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T))
obsBM <- apply (obsBM, 1, mean)
obsEB<-do.call(cbind,sapply(RAO_EB, "[","Observado",simplify=T))
obsEB <- apply (obsEB, 1, mean)
obsOU<-do.call(cbind,sapply(RAO_OU, "[","Observado",simplify=T))
obsOU <- apply (obsOU, 1, mean)

##
data_to_map_emp_sim <- data.frame(longlat,
                                  'Observed' = RAO_OBS$Observado,
                                  'Null Disparity' = RAO_OBS$med_nulo,
                                  'Simulated BM' = obsBM,
                                  'Simulated EB' = obsEB,
                                  'Simulated OU' = obsOU)

# melt
melt_data_to_map_emp_sim <- melt(data_to_map_emp_sim,id=c("LONG","LAT"))
colnames(melt_data_to_map_emp_sim)[which(colnames(melt_data_to_map_emp_sim) == "value")] <- 'Disparity'

# and map
#plot using ggplot

alternative_map1 <- ggplot(melt_data_to_map_emp_sim, 
                           aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = Disparity)) +
  facet_wrap(~variable,scales = "fixed",ncol=3)+
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
        strip.text = element_text(size=15))

alternative_map1

## correlation between average null and simulated by OU

cor (data.frame (RAO_OBS$med_nulo, obsBM,obsEB,obsOU))
