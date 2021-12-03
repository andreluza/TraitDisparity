# ------------------------# 

# Run script of results

# ---------------------------------------------------------- #
# these results are produced by simulation traits across all 413 spp
# MAIN RESULTS

source ("./R/Interpretation_Fig2_and_statistics_413spp.R")

# to obtain the number of cells with  empirical disparity either higher or lower than
# the null/evol model disparity
# (significance in terms of to have values more extreme than 1.96)

count_cells_relative_to_models

# get the panel with the relationship between MPD and DISPARITY

# to save run
pdf (here ("output","vectorized","Fig2_all.pdf"), width=6,height=6,family="serif")

plot(panel)

dev.off()

# png format
png (here ("output","vectorized","Fig2_all.png"), 
     width = 12, height = 12, units = "cm",
     res=300,family="serif")

plot(panel)

dev.off()

# one model with average values across all phylogenies
# or a table (only one each time)
tab_model (av_coeff)

# R2
av_coeff$adj.r.squared

# get statistics (model estimates )
## average across phylogenies

list.mean

## associated sd

list.sd

## Fstatistics

apply(do.call(rbind, F_stat),2,mean)
apply(do.call(rbind, F_stat),2,sd)

# R2 adj
mean(unlist(R2_stat))
sd(unlist(R2_stat))

# get the maps
source ("./R/Interpretation_Fig3_413spp.R")

# if you want to save then run
pdf (here ("output","vectorized","Fig3_all.pdf"), 
     width=7,height=10,family="serif")

plot(panel2)

dev.off()


# if you want to save then run
png (here ("output","vectorized","Fig3_all.png"), 
     width = 13, height = 25, units = "cm",
     res=300, family="serif")

plot(panel2)

dev.off()


# OBSERVED VERSUS SIMULATED PATTERNS
pdf (here ("output","vectorized","Fig4_all_obs_simulated.pdf"), 
     width=8,height=8,family="serif")

plot(alternative_map1)

dev.off()

## correlation between average null and simulated by OU

cor (data.frame (RAO_OBS$med_nulo, obsBM,obsEB,obsOU))
