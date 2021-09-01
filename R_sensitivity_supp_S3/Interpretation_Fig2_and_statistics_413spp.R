# ------------------------# 

# load packages
source ("./R/packages.R")

# -----------------------#
# load data
load (here("output_sensitivity_supp_S3", "mpd_results_ALL.RData"))
load (here("output_sensitivity_supp_S3", "RAO_BM_ALL.RData"))
load (here("output_sensitivity_supp_S3", "RAO_EB_ALL.RData"))
load (here("output_sensitivity_supp_S3", "RAO_OU_ALL.RData"))
load (here("output_sensitivity_supp_S3", "RAO_OBS_ALL.RData"))
load (here("output_sensitivity_supp_S3", "parameters.RData"))

## set the names to output (sigmas, betas, alphas)
#BM
names(RAO_BM)<-paste ("sigma=",sigmas,sep="")
# EB
names(RAO_EB)<-paste ("sigma=",sigmas,sep="")
RAO_EB <- lapply (RAO_EB, function (i) {

	names(i) <- paste ("beta=",betas,sep="");
	i

})

# OU
names(RAO_OU)<-paste ("sigma=",sigmas,sep="")
RAO_OU <- lapply (RAO_OU, function (i) {

	names(i) <- paste ("alpha=",alphas,sep="");
	i

})

# --------------------------------------------------------------- #
# relationship between empirical and simulated data sets
# ----------------------------------------------------------------- #

# remove estimates impossible to be done

RAO_BM_sub <- RAO_BM[-which(lapply (RAO_BM, function (i)  sum(is.na(i$SES)))>0)]# the calculation with sigma zero did not work
RAO_EB_sub<- RAO_EB[2:5] # the calculation with sigma zero did not work
RAO_OU_sub<- RAO_OU[2:5] # the calculation with sigma zero did not work

# also params

sigmas_sub <- sigmas [2:5]
betas_sub <- betas [2:5]
alphas_sub <- alphas [2:5]

# each simulated sigma
test_sigmas <- rbind (cbind(RAO_OBS$SES,"EMPIRICAL"),
                      cbind (RAO_BM_sub[[1]]$SES,paste("sigma",sigmas_sub[1],sep="")),
                      cbind (RAO_BM_sub[[2]]$SES,paste("sigma",sigmas_sub[2],sep="")),
                      cbind (RAO_BM_sub[[3]]$SES,paste("sigma",sigmas_sub[3],sep="")),
                      cbind (RAO_BM_sub[[4]]$SES,paste("sigma",sigmas_sub[4],sep="")))
#bind MPD
test_sigmas<-cbind(test_sigmas, null.mpdf$SES.MPD)
test_sigmas <- data.frame(test_sigmas)# turn into df
colnames(test_sigmas)<-c("SES.Disparity","Dataset","SES.MPD")
test_sigmas$SES.Disparity <-as.numeric(test_sigmas$SES.Disparity)# change class
test_sigmas$Dataset <-as.factor(test_sigmas$Dataset)
test_sigmas$SES.MPD <-as.numeric(test_sigmas$SES.MPD)# change class

# plot
figure1 <- ggplot (data = test_sigmas, aes (x=SES.MPD, 
                                            y=SES.Disparity,
                                            colour=Dataset)) + 
  geom_point(aes (col=Dataset),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("EMPIRICAL"="black", 
                               "sigma0.5"="#BECA5C",# blue
                               "sigma1"="#EF8D32",
                               "sigma1.5"="#CC561E",
                               "sigma2"="#AA2B1D")) + 
  xlab("") + 
  ylab ("") + 
  theme (
    title = element_text(size=7),
    axis.title = element_text(size=11),
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    plot.margin=unit(c(.2,1,.1,1),"cm")) #+ 
  #ggtitle(label = paste ("sigma=", sigmas_sub[s], ", beta=",betas[b],", alpha=", alphas[a]),
  
# each simulated beta, with sigma=1and varied beta
test_betas <- rbind (cbind(RAO_OBS$SES,"EMPIRICAL"),
                      cbind (RAO_EB_sub[[1]][[1]]$SES,paste("beta",betas_sub[1],sep="")),
                      cbind (RAO_EB_sub[[1]][[2]]$SES,paste("beta",betas_sub[2],sep="")),
                      cbind (RAO_EB_sub[[1]][[3]]$SES,paste("beta",betas_sub[3],sep="")),
                      cbind (RAO_EB_sub[[1]][[4]]$SES,paste("beta",betas_sub[4],sep="")))
#bind MPD
test_betas<-cbind(test_betas, null.mpdf$SES.MPD)
test_betas <- data.frame(test_betas)# turn into df
colnames(test_betas)<-c("SES.Disparity","Dataset","SES.MPD")
test_betas$SES.Disparity <-as.numeric(test_betas$SES.Disparity)# change class
test_betas$Dataset <-as.factor(test_betas$Dataset)
test_betas$SES.MPD <-as.numeric(test_betas$SES.MPD)# change class

# plot
figure2 <- ggplot (data = test_betas, aes (x=SES.MPD, 
                                            y=SES.Disparity,
                                            colour=Dataset)) + 
  geom_point(aes (col=Dataset),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("EMPIRICAL"="black", 
                               "beta-0.5"="#BECA5C",# blue
                               "beta0"="#EF8D32",
                               "beta0.5"="#CC561E",
                               "beta1"="#AA2B1D")) + 
  xlab("") + 
  ylab ("SES of Morphological Disparity") + 
  theme (
    title = element_text(size=7),
    axis.title = element_text(size=11),
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    plot.margin=unit(c(.2,1,.1,1),"cm")) #+ 
#ggtitle(label = paste ("sigma=", sigmas_sub[s], ", beta=",betas[b],", alpha=", alphas[a]),

# each simulated beta, with sigma=1and varied alpha
test_alphas <- rbind (cbind(RAO_OBS$SES,"EMPIRICAL"),
                     cbind (RAO_OU_sub[[1]][[1]]$SES,paste("alpha",alphas_sub[1],sep="")),
                     cbind (RAO_OU_sub[[1]][[2]]$SES,paste("alpha",alphas_sub[2],sep="")),
                     cbind (RAO_OU_sub[[1]][[3]]$SES,paste("alpha",alphas_sub[3],sep="")),
                     cbind (RAO_OU_sub[[1]][[4]]$SES,paste("alpha",alphas_sub[4],sep="")))
#bind MPD
test_alphas<-cbind(test_alphas, null.mpdf$SES.MPD)
test_alphas <- data.frame(test_alphas)# turn into df
colnames(test_alphas)<-c("SES.Disparity","Dataset","SES.MPD")
test_alphas$SES.Disparity <-as.numeric(test_alphas$SES.Disparity)# change class
test_alphas$Dataset <-as.factor(test_alphas$Dataset)
test_alphas$SES.MPD <-as.numeric(test_alphas$SES.MPD)# change class

# plot
figure3 <- ggplot (data = test_alphas, aes (x=SES.MPD, 
                                           y=SES.Disparity,
                                           colour=Dataset)) + 
  geom_point(aes (col=Dataset),alpha=0.4,size=1.5) + 
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) + 
  theme_classic() + 
  scale_colour_manual(values=c("EMPIRICAL"="black", 
                               "alpha1"="#BECA5C",# blue
                               "alpha1.5"="#EF8D32",
                               "alpha2"="#CC561E",
                               "alpha2.5"="#AA2B1D")) + 
  xlab("SES of Mean Pairwise Distance") + 
  ylab ("") + 
  theme (
    title = element_text(size=7),
    axis.title = element_text(size=11),
    legend.text = element_text(size=10),
    legend.title = element_text(size=11),
    plot.margin=unit(c(.2,1,.1,1),"cm")) #+ 
#ggtitle(label = paste ("sigma=", sigmas_sub[s], ", beta=",betas[b],", alpha=", alphas[a]),


# create a panel


panel<-grid.arrange (figure1,
              figure2,
              figure3,
              ncol=1,
              nrow=3)
panel


# end