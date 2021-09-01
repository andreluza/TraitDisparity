# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_concensus_supp_S2", "mpd_results_ORYZ.RData"))
load (here("output_concensus_supp_S2", "RAO_BM_ORYZ.RData"))
load (here("output_concensus_supp_S2", "RAO_EB_ORYZ.RData"))
load (here("output_concensus_supp_S2", "RAO_OU_ORYZ.RData"))
load (here("output_concensus_supp_S2", "RAO_OBS_ORYZ.RData"))

# average observed disparity for each  model
avObsBM<- RAO_BM$Observado
avObsEB<- RAO_EB$Observado
avObsOU<- RAO_OU$Observado

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    lowerNULL=table(RAO_OBS$SES<=-1.96)[2],
    lowerOU=table(((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU)) <=-1.96)[2],
    lowerBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) <=-1.96)[2],
    lowerEB=table(((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB)) <=-1.96)[2]
  ), 
  # higher than
  higher = cbind (
    highNULL=table(RAO_OBS$SES>=1.96)[2],
    highOU=table(((RAO_OBS$Observado - mean(avObsOU))/sd(avObsOU)) >=1.96)[2],
    highBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) >=1.96)[2],
    highEB=table(((RAO_OBS$Observado - mean(avObsEB))/sd(avObsEB)) >=1.96)[2]
  )
)

# proportion of cells in each group
count_cells_relative_to_models


#--------------------------------------------------------
# relationship between empirical and simulated data sets

# working with the average

avBM<-RAO_BM$SES
avEB<-RAO_EB$SES
avOU<- RAO_OU$SES
avMPD<- statistics.phy$SES.MPD

# df average disparity
av_disp <- rbind(
  data.frame (value=RAO_OBS$SES,Data="Empirical"),
  data.frame (value= (avBM),Data="BM"),
  data.frame (value= (avEB),Data="EB"),
  data.frame (value= (avOU),Data="OU"))

# bind MPD
av_disp <- cbind(av_disp, MPD =  (avMPD))

# empirical as the first level
av_disp$Data <- factor (av_disp$Data,
                        levels = c("Empirical",
                                   "BM",
                                   "EB",
                                   "OU"))

# coefficients for av estimates

av_coeff <- summary(lm (value ~ MPD*Data,
                        data=av_disp))

# reported in the main text (Table 1)
tab_model(av_coeff)

# figure
figure1 <- ggplot (data = av_disp, aes (x=MPD, y=value,colour=Data)) + 
  geom_point(aes (col=Data),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("Empirical"="#9dab86", 
                               "BM"="#21209c",
                               "EB"="#fdb827",
                               "OU"="#23120b")) + 
  xlab("SES of Mean Pairwise Distance") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    legend.position = "top",
    plot.margin=unit(c(.2,1,.1,1),"cm"))


figure1
