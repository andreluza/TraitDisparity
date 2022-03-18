# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_concensus_supp_S3", "mpd_results_ALL.RData"))
load (here("output_concensus_supp_S3", "RAO_BM_ALL_multivariate.RData"))
load (here("output_concensus_supp_S3", "RAO_OBS_ALL_multivariate.RData"))

# average observed disparity for each  model
avObsBM<- rowMeans (do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T)))
nullRao <- rowMeans (do.call(cbind,sapply(RAO_BM, "[","med_nulo",simplify=T)))
avMPD<- statistics.phy$SES.MPD

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    lowerNULL=table(RAO_OBS$SES<=-1.96)[2],
    lowerBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) <=-1.96)[2]
  ), 
  # higher than
  higher = cbind (
    highNULL=table(RAO_OBS$SES>=1.96)[2],
    highBM=table(((RAO_OBS$Observado - mean(avObsBM))/sd(avObsBM)) >=1.96)[2]
  )
)

# proportion of cells in each group
count_cells_relative_to_models

#--------------------------------------------------------
# relationship between empirical and simulated data sets

# working with the average

avBM<-do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T))

# df average disparity
av_disp <- rbind(
  data.frame (value=RAO_OBS$SES,Data="Empirical"),
  data.frame (value=rowMeans (avBM),Data="BM"))

# bind MPD
av_disp <- cbind(av_disp, MPD =  (avMPD))

# empirical as the first level
av_disp$Data <- factor (av_disp$Data,
                        levels = c("Empirical",
                                   "BM"))

# coefficients for av estimates

av_coeff <- (lm (value ~ MPD*Data,
                 data=av_disp))

# compare slopes
m.lst.FEve <- emtrends (av_coeff,  "Data", var="MPD")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# reported in the main text (Table 1)
summary(av_coeff)

# difference in slope between simulated and empirical data
m.lst_tab.FEve

# test of whether slope of the regression between ses disparity and ses MPD
# varies across evolutionary models
# reported in the main text (Table 1)
tab_model(av_coeff)


# figure
figure1 <- ggplot (data = av_disp, aes (x=MPD, y=value,colour=Data)) + 
  geom_point(aes (col=Data),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("Empirical"="#9dab86", 
                               "BM"="#21209c")) + 
  xlab("SES of Mean Pairwise Distance") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    legend.position = "top",
    plot.margin=unit(c(.2,1,.1,1),"cm"))

# 


