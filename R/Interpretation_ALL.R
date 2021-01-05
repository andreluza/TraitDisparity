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

# df average disparity
av_disp <- rbind(
  data.frame (value=RAO_OBS$SES,Data="Empirical"),
  data.frame (value=rowMeans (avBM),Data="BM"),
  data.frame (value=rowMeans (avEB),Data="EB"),
  data.frame (value=rowMeans (avOU),Data="OU"))
  
# bind MPD
av_disp <- cbind(av_disp, MPD = rowMeans (avMPD))

# rm NA (sites lacking oryzomyialia)
av_disp<- av_disp[which(is.na(av_disp$MPD) != T),]

#
figure1 <- ggplot (data = av_disp, aes (x=MPD, y=value,colour=Data)) + 
  geom_point(aes (col=Data),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("Empirical"="#9dab86", 
                               "BM"="#21209c",
                               "EB"="#fdb827",
                               "OU"="#23120b")) + 
  xlab("SES of Mean pairwise distance") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    legend.position = "none",
    plot.margin=unit(c(.2,1,.1,1),"cm"))

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

# revoving sites lacking oryzomyialia
bind_data_disparity <- lapply (bind_data_disparity, function (i)
  i[which(is.na(i$MPD) != T),])

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
  # load also the model results
  res <- list (coef = res,
               model = m1)
  ; # return
  res
  
})

# transforming res list into df
model_slope_df <- do.call (rbind, sapply (model_slope,"[","coef"))

# finally melt

model_slope_df <- melt(model_slope_df)

# plotting

fig2 <- ggplot(model_slope_df, aes(x=value, color=variable, fill=variable)) +
  geom_density(size=1,alpha=0.1) +
  scale_fill_manual(values=c("EMPIRICAL"="#9dab86", 
                             "BM"="#21209c",
                             "EB"="red",
                             "OU"="#23120b")) + 
  scale_colour_manual(values=c("EMPIRICAL"="#9dab86", 
                               "BM"="#21209c",
                               "EB"="#fdb827",
                               "OU"="#23120b"))+
  theme_classic()  

fig2a <- fig2+ theme (
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

boxplotOver <- ggboxplot(model_slope_df, x = "variable", y = "value",
                         color = "variable", palette =c("#9dab86", 
                                                        "#21209c",
                                                        "#fdb827",
                                                        "#23120b"),
                         fun = "mean_ci", add.params = list(alpha=0.2),size=0.5)

boxplotOver1 <- boxplotOver + theme(axis.line=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.title.y=element_blank(),
                                    legend.position="none",
                                    panel.background=element_blank(),
                                    panel.border=element_blank(),
                                    panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    plot.background=element_blank(),
                                    plot.margin=unit(c(.2,1,.1,1),"cm")) +
  scale_y_continuous(limits=c(-1.2, 1.2)) + coord_flip() 


#pdf (paste ("compTR",variaveis[i],".pdf",sep=""), width=5,height=4,family="serif")
grid.arrange(figure1,
             boxplotOver1,fig2a,
             ncol=4,nrow = 9,
             layout_matrix = rbind (c(1,1,1,1),
                                    c(1,1,1,1),
                                    c(1,1,1,1),
                                    c(1,1,1,1),
                                    c(2,2,2,2),
                                    c(3,3,3,3),
                                    c(3,3,3,3),
                                    c(3,3,3,3),
                                    c(3,3,3,3)))

#dev.off()


## Fstatistics

F_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$fstatistic)
apply(do.call(rbind, F_stat),2,mean)
apply(do.call(rbind, F_stat),2,sd)

# R2 adj
R2_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$adj.r.squared)
mean(unlist(R2_stat))
sd(unlist(R2_stat))

## average results of the models

coef_table <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$coefficients)
## mean
list.mean <- Reduce("+", coef_table) / length(coef_table)
## squared mean
list.squared.mean <- Reduce("+", lapply(coef_table, "^", 2)) / length(coef_table)
## variance
list.variance <- list.squared.mean - list.mean^2
## standard deviation
list.sd <- sqrt(list.variance)

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

# data to MAP

data_to_map <- data.frame(longlat,
                          #Richness=rowSums(presab),
                          'EMPIRICAL'= RAO_OBS$SES,
                          'BROWNIAN MOTION'=rowMeans (avBM),
                          'EARLY BURST' =rowMeans (avEB),
                          'ORNSTEIN-UHLENBECK'=rowMeans (avOU)
                          #'NEUTRAL' = RAO_OBS$SES-rowMeans (avEB)
                          #MPD=rowMeans (avMPD)
)

#write.csv(data_to_map,here("output","data_to_SAM.csv"))

melt_data_to_map <- melt(data_to_map,id=c("LONG","LAT"))
colnames(melt_data_to_map)[which(colnames(melt_data_to_map) == "value")] <- 'SES'

## empirical vs simulated
#plot using ggplot
map1 <- ggplot(melt_data_to_map, 
               aes(x = LONG, y = LAT)) +
  geom_tile(aes(fill = SES)) +
  facet_wrap(~variable,scales = "fixed",ncol=2)+
  scale_fill_gradient2(midpoint = 0, mid="#eee8d5", high="#dc322f", low="#268bd2") + 
  theme_classic() + 
  theme_map() +
  xlab("Longitude") + ylab("Latitude")+
  theme(legend.position="right",
        legend.justification = "center",
        legend.direction = "vertical",
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=10))
map1

ggsave (here("output","maps_ALL.pdf"),height = 6,width=5)

