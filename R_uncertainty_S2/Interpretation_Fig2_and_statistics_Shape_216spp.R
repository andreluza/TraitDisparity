# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_uncertainty_S2", "mpd_results_216.RData"))
load (here("output_uncertainty_S2", "RAO_OBS_216_multivariate.RData"))
load (here("output_uncertainty_S2", "RAO_BM_216_multivariate.RData"))

# average observed disparity for each  model
avObsBM<-rowMeans (do.call(cbind,sapply(RAO_BM, "[","Observado",simplify=T)))
nullRao <- rowMeans (do.call(cbind,sapply(RAO_BM, "[","med_nulo",simplify=T)))

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

#--------------------------------------------------------
# relationship between empirical and simulated data sets

# working with the average

avBM<-do.call(cbind,sapply(RAO_BM, "[","SES",simplify=T))
avMPD<-do.call(cbind,sapply(null.mpdf, "[","SES.MPD",simplify=T))

# df average disparity
av_disp <- rbind(
  data.frame (value=RAO_OBS$SES,Data="Empirical"),
  data.frame (value=rowMeans (avBM),Data="BM"))
  
# bind MPD
av_disp <- cbind(av_disp, MPD = rowMeans (avMPD))

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
                    BM=abs(m1$coefficients[2]+m1$coefficients[4])
                    #R2 = RsquareAdj(m1)$adj.r.squared
  )
  
  # load also the model results
  res <- list (coef = res,
               model = m1)
  ; # return
  res
  
})

# average coefficients across simulations
#av_across_sim <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$coefficients)
#estimates_across_sim <- do.call(rbind,
 #                               lapply(av_across_sim, function (i) # extract estimates and melt
  #                                i[,1] # the column of estimates
   #                             ))

#data.frame (mean.coefficient = apply (estimates_across_sim,2,mean),
#            sd.coefficient = apply (estimates_across_sim,2,sd))


####-------------------------# 
## Fstatistics

F_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$fstatistic)
#apply(do.call(rbind, F_stat),2,mean)
#apply(do.call(rbind, F_stat),2,sd)

# R2 adj
R2_stat <- lapply(sapply(model_slope,'[',"model"), function (i) summary (i)$adj.r.squared)
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
                               "BM"="#21209c")) + 
  xlab("SES of Mean Pairwise Distance") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    legend.position = "none",
    plot.margin=unit(c(.2,1,.1,1),"cm"))

# density plot
fig2 <- ggplot(model_slope_df, aes(x=value, color=variable, fill=variable)) +
  geom_density(size=1,alpha=0.1) +
  scale_fill_manual(values=c("EMPIRICAL"="#9dab86", 
                             "BM"="#21209c")) + 
  scale_colour_manual(values=c("EMPIRICAL"="#9dab86", 
                                        "BM"="#21209c"))+
  theme_classic()  

fig2a <- fig2+ theme (
               axis.text.x = element_text(angle = 90),
               legend.title = element_blank(),
               #legend.position = "none",
               legend.position = c(1, .95),
               legend.justification = c("right", "top"),
               legend.box.just = "right",
               legend.margin = margin(6, 6, 6, 6),
               plot.margin=unit(c(.2,1,.1,1),"cm")) +
  xlab ("Slope SES Disparity ~ SES MPD") + 
  ylab("Density") + 
  scale_x_continuous(breaks = seq(-0.2,1.3,0.15),
                     limits = c(-0.2,1.3))

# exploring phylogenetic uncertainty

boxplotOver <- ggboxplot(model_slope_df, x = "variable", y = "value",
               color = "variable", palette =c("#9dab86", 
                                              "#21209c"),
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
                plot.margin=unit(c(.2,1,.1,1.5),"cm")) +
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

