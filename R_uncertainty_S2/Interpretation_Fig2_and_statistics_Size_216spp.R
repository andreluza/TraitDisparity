# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_uncertainty_S2", "mpd_results_216.RData"))
load (here("output_uncertainty_S2", "RAO_BM_216.RData"))
load (here("output_uncertainty_S2", "RAO_EB_216.RData"))
load (here("output_uncertainty_S2", "RAO_OU_216.RData"))
load (here("output_uncertainty_S2", "RAO_OBS_216.RData"))


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

# empirical as the first level
av_disp$Data <- factor (av_disp$Data,
                        levels = c("Empirical",
                                   "BM",
                                   "EB",
                                   "OU"))

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
  # intercept (average SES disparity in the empirical dataset (level "AOBS"))
  # coeff - difference between slope of SES ~MPD of empirical and each evo model
  res <- data.frame(EMPIRICAL=abs(m1$coefficients[2]), # EMpirical
                    BM=abs(m1$coefficients[2]+m1$coefficients[6]), # BM
                    EB=abs(m1$coefficients[2]+m1$coefficients[7]), # EB
                    OU=abs(m1$coefficients[2]+m1$coefficients[8])#,OU
                    #R2 = RsquareAdj(m1)$adj.r.squared
  )
  
  # load also the model results
  res <- list (coef = res,
               model = m1)
  ; # return
  res
  
})

####-------------------------# 
## Fstatistics
#F_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$fstatistic)
#apply(do.call(rbind, F_stat),2,mean)
#apply(do.call(rbind, F_stat),2,sd)
# R2 adj
#R2_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$adj.r.squared)
#mean(unlist(R2_stat))
#sd(unlist(R2_stat))

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

# transforming res list into df
model_slope_df <- do.call (rbind, sapply (model_slope,"[","coef"))

# finally melt
model_slope_df <- melt(model_slope_df)

# plotting

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
    legend.position = "none",
    plot.margin=unit(c(.2,1,.1,1),"cm"))

# density plot
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
  scale_x_continuous(breaks = seq(-0.7,1.6,0.15),
                     limits = c(-0.7,1.6))

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
                plot.background=element_rect(fill="white",colour="white"),
                panel.background=element_rect(fill="white",colour="white"),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.margin=unit(c(.2,1,.1,2),"cm")) +
  scale_y_continuous(limits=c(-0.7,1.6)) + coord_flip() 

# arrange

panel <- grid.arrange(figure1,
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


# show
list.mean
list.sd

# parameter values (fitContinuous)

load(here ("output_uncertainty_S2", "params_fitcontinuous_216.RData"))

# params
# BM
df_sigma_BM <- data.frame (Estimates = unlist(lapply (simul_param_BM, function (i) i$opt$sigsq)),
                           Parameter = "sigma",
                           model = "BM")
df_sigma_EB <- data.frame (Estimates = unlist(lapply (simul_param_EB, function (i) i$opt$sigsq)),
                           Parameter = "sigma",
                           model = "EB")
df_beta_EB <- data.frame (Estimates = unlist(lapply (simul_param_EB, function (i) i$opt$a)),
                          Parameter = "beta", 
                          model = "EB")
df_sigma_OU <- data.frame (Estimates = unlist(lapply (simul_param_OU, function (i) i$opt$sigsq)),
                           Parameter = "sigma",
                           model = "OU")
df_alpha_OU <- data.frame (Estimates = unlist(lapply (simul_param_OU, function (i) i$opt$alpha)),
                           Parameter = "alpha",
                           model = "OU")

# rbind
df_density <- rbind (df_sigma_BM,
                     df_sigma_EB,
                     df_beta_EB,
                     df_sigma_OU,
                     df_alpha_OU)
# plot
# density plot
fig_params_BM <- ggplot(df_density[which(df_density$model == "BM"),], 
                     aes(x=Estimates,
                         group=Parameter,
                         color=Parameter,
                         fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("sigma" = "red"))+
  scale_colour_manual(values=c("sigma" = "red"))+
  theme_classic()  

# density plot
fig_params_EB <- ggplot(df_density[which(df_density$model == "EB" & 
                                           df_density$Parameter == "sigma"),], 
                        aes(x=Estimates,
                            group=Parameter,
                            color=Parameter,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("sigma" = "red"))+
  scale_colour_manual(values=c("sigma" = "red"))+
  theme_classic()  

# beta
fig_params_EB_beta <- ggplot(df_density[which(df_density$model == "EB" & 
                                           df_density$Parameter == "beta"),], 
                        aes(x=Estimates,
                            group=Parameter,
                            color=Parameter,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("beta" = "green"))+
  scale_colour_manual(values=c("beta" = "green"))+
  theme_classic()  + 
  #theme (legend.position = c(-1.0010e-06,3e+09))+
  
  xlim (c(min(df_density[which(df_density$model == "EB" & 
                                 df_density$Parameter == "beta"),"Estimates"]),
          max(df_density[which(df_density$model == "EB" & 
                                 df_density$Parameter == "beta"),"Estimates"])))


# density plot
fig_params_OU <- ggplot(df_density[which(df_density$model == "OU" & 
                                           df_density$Parameter == "sigma"),], 
                        aes(x=Estimates,
                            group=Parameter,
                            color=Parameter,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("sigma" = "red"))+
  scale_colour_manual(values=c("sigma" = "red"))+
  theme_classic()  

# density plot
fig_params_OU_alpha <- ggplot(df_density[which(df_density$model == "OU" & 
                                           df_density$Parameter == "alpha"),], 
                        aes(x=Estimates,
                            group=Parameter,
                            color=Parameter,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.1)+
  scale_fill_manual(values=c("alpha" = "brown"))+
  scale_colour_manual(values=c("alpha" = "brown"))+
  theme_classic()   + theme(legend.position = "none")

# open parameters of the multivariate model

load(here ("output_uncertainty_S2", "params_BM_216_multivariate.RData"))

df_multivariate <- unlist(lapply(simul_param_BM, function (i) i$sigma$Pinv))
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
  theme(legend.position = "none")

# arrange
panel_params <- grid.arrange(fig_params_BM+ theme(legend.position = c(0.8,0.8),
                                  legend.title = element_blank(),
                                  axis.title.x = element_blank()),
             fig_params_EB+ theme(legend.position = c(0.8,0.8),
                                  legend.title = element_blank(),
                                  axis.title = element_blank()),
             fig_params_EB_beta+ theme(legend.position = c(0.4,0.8),
                                       legend.title = element_blank(),
                                       axis.text.x = element_text(size=7),
                                       axis.text.y = element_text(size=7),
                                       axis.title = element_blank()),
             fig_params_OU+ theme(legend.position = c(0.8,0.8),
                                  legend.title = element_blank(),
                                  axis.title = element_blank()),
             fig_params_OU_alpha+ theme(legend.position = c(0.8,0.8),
                                        legend.title = element_blank(),
                                        axis.title = element_blank()),
             fig_params_multiv + theme (legend.position = c(0.6,0.8),
                                        legend.title = element_blank()
             ),
  
  ncol=3,nrow = 3,
  layout_matrix = rbind (c(1,2,4),
                         c(NA,3,5),
                         c(6,6,6)))


