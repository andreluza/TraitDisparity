# ------------------------# 

# load packages
source ("./R/packages.R")

# show estimated sigma
load (here("output_concensus_supp_S5", "params_BM.RData"))

hist_sigma <- hist(simul_params_BM$sigma$S,
     xlab="Sigma",
     main="Values of Sigma across the 112 traits",col="white")
abline (v= mean(simul_params_BM$sigma$S),lty=2,lwd=3, col= "black")

# -----------------------#
# load data
load (here("output", "mpd_results_ALL.RData"))
load (here("output_concensus_supp_S5", "RAO_OBS.RData"))
load (here("output_concensus_supp_S5", "RAO_BM.RData"))

# average observed disparity for each  model
avObsBM<-RAO_BM$Observado #rowMeans (do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T)))
#avObsEB<-rowMeans (do.call(cbind,sapply(RAO_EB, "[","Observado",simplify=T)))
#avObsOU<-rowMeans (do.call(cbind,sapply(RAO_OU, "[","Observado",simplify=T)))
nullRao <- RAO_BM$med_nulo

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    lowerNULL=table(RAO_OBS$SES<=-1.96)[2],
    #lowerOU=table(((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU)) <=-1.96)[2],
    lowerBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) <=-1.96)[2]#,
    #lowerEB=table(((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB)) <=-1.96)[2]
  ), 
  # higher than
  higher = cbind (
    highNULL=table(RAO_OBS$SES>=1.96)[2],
    #highOU=table(((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU)) >=1.96)[2],
    highBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) >=1.96)[2]#,
    #highEB=table(((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB)) >=1.96)[2]
  )
)
count_cells_relative_to_models

# relationship between empirical and simulated data sets

# working with uncertainty

avBM<-RAO_BM$SES# do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T))
#avEB<-do.call(cbind,sapply(RAO_EB, "[","SES",simplify=T))
#avOU<-do.call(cbind,sapply(RAO_OU, "[","SES",simplify=T))
avMPD<-do.call(cbind,sapply(null.mpdf, "[","SES.MPD",simplify=T))

# working with the average
# df average disparity
av_disp <- rbind(
  data.frame (value=RAO_OBS$SES,Data="Empirical"),
  data.frame (value= (avBM),Data="BM"),
  data.frame (value= (nullRao),Data="NULL")#,
  
  #data.frame (value=rowMeans (avEB),Data="EB"),
  #data.frame (value=rowMeans (avOU),Data="OU")
  )

  
# bind MPD
av_disp <- cbind(av_disp, MPD = rowMeans (avMPD))

# empirical as the first level
av_disp$Data <- factor (av_disp$Data,
        levels = c("Empirical",
                   "BM",
                   "NULL")#,
                   #"EB",
                   #"OU")
                   )

# coefficients for av estimates

av_coeff <- summary(lm (value ~ MPD*Data,
    data=av_disp))


# figure of the average
figure1 <- ggplot (data = av_disp, aes (x=MPD, y=value,colour=Data)) + 
  geom_point(aes (col=Data),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE,size =1.5) + 
  theme_classic() + 
  scale_colour_manual(values=c("Empirical"="#9dab86", 
                               "BM"="#21209c",
                               "NULL" = "#fdb827"
                               #"EB"="#fdb827",
                               #"OU"="#23120b"
                               )) + 
  xlab("SES of Mean Pairwise Distance") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    legend.position = "top",
    plot.margin=unit(c(.2,1,.1,1),"cm"))



