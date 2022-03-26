# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_concensus_supp_S3", "mpd_results_ALL.RData"))
load (here("output_concensus_supp_S3", "RAO_BM_ALL.RData"))
load (here("output_concensus_supp_S3", "RAO_EB_ALL.RData"))
load (here("output_concensus_supp_S3", "RAO_OU_ALL.RData"))
load (here("output_concensus_supp_S3", "RAO_OBS_ALL.RData"))

# average observed disparity for each  model
avObsBM.mean<-apply (do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T)),1,mean)
avObsBM.sd<-apply (do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T)),1,sd)
avObsEB.mean<-apply (do.call(cbind,sapply(RAO_EB, "[","Observado",simplify=T)),1,mean)
avObsEB.sd<-apply (do.call(cbind,sapply(RAO_EB, "[","Observado",simplify=T)),1,sd)
avObsOU.mean<-apply (do.call(cbind,sapply(RAO_OU, "[","Observado",simplify=T)),1,mean)
avObsOU.sd<-apply (do.call(cbind,sapply(RAO_OU, "[","Observado",simplify=T)),1,sd)

# counting cells with significant neutral SES, observed relative to the models
# lower than
count_cells_relative_to_models <- list (
  lower = cbind (
    lowerNULL=table(RAO_OBS$SES<=-1.96)[2],
    lowerOU=table(((RAO_OBS$Observado - avObsOU.mean)/avObsOU.sd) <=-1.96)[2],
    lowerBM=table(((RAO_OBS$Observado - avObsBM.mean)/avObsBM.sd) <=-1.96)[2],
    lowerEB=table(((RAO_OBS$Observado - avObsEB.mean)/avObsEB.sd) <=-1.96)[2]
  ), 
  # higher than
  higher = cbind (
    highNULL=table(RAO_OBS$SES>=1.96)[2],
    highOU=table(((RAO_OBS$Observado - avObsOU.mean)/avObsOU.sd) >=1.96)[2],
    highBM=table(((RAO_OBS$Observado - avObsBM.mean)/avObsBM.sd) >=1.96)[2],
    highEB=table(((RAO_OBS$Observado - avObsEB.mean)/avObsEB.sd) >=1.96)[2]
  )
)

# proportion of cells in each group
count_cells_relative_to_models


#--------------------------------------------------------
# relationship between empirical and simulated data sets

# working with the average

avBM<-rowMeans (do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T)))
avEB<-rowMeans (do.call(cbind,sapply(RAO_EB, "[","SES",simplify=T)))
avOU<- rowMeans (do.call(cbind,sapply(RAO_OU, "[","SES",simplify=T)))
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



# parameter values (fitContinuous)

load(here ("output_concensus_supp_S3", "params_fitcontinuous.RData"))

# params
# BM
df_sigma_BM <- data.frame (Estimates = round(as.numeric(simul_param_BM$opt$sigsq),3),
                           Parameter = "sigma",
                           model = "BM")
df_sigma_EB <- data.frame (Estimates = round(as.numeric(simul_param_EB$opt$sigsq),3),
                           Parameter = "sigma",
                           model = "EB")
df_beta_EB <- data.frame (Estimates = round(as.numeric(simul_param_EB$opt$a),3),
                          Parameter = "beta", 
                          model = "EB")
df_sigma_OU <- data.frame (Estimates = round(as.numeric(simul_param_OU$opt$sigsq),3),
                           Parameter = "sigma",
                           model = "OU")
df_alpha_OU <- data.frame (Estimates = round(as.numeric(simul_param_OU$opt$alpha),3),
                           Parameter = "alpha",
                           model = "OU")
# rbind
df_density <- rbind (df_sigma_BM,
                     df_sigma_EB,
                     df_beta_EB,
                     df_sigma_OU,
                     df_alpha_OU)
# open parameters of the multivariate model

load(here ("output_concensus_supp_S3", "params_BM_multivariate.RData"))

df_multivariate <- as.numeric(simul_param_BM$sigma$Pinv)
df_multivariate<-data.frame (Estimates=df_multivariate,
                             Parameter = "sigma")
fig_params_multiv <- ggplot(df_multivariate, 
                              aes(x=Estimates,
                                  group=Parameter,
                                  color=Parameter,
                                  fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("sigma" = "red"))+
  scale_colour_manual(values=c("sigma" = "red"))+
  theme_classic()   + 
  theme(legend.position = c(0.7,0.5))
