# ------------------------# 

# Run scripts of results

# these  results consider the 216 species in trait, comm, and phy datasets

source ("./R/Interpretation_Fig2_and_statistics_216_spp.R")

# get the panel with the relationship between MPD and DISPARITY

# to save run
# pdf (here ("output","Fig2_216spp.pdf"), width=6,height=7,family="serif")

plot(panel)

# dev.off()

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
source ("./R/Interpretation_Fig3_216_spp.R")

# if you want to save then run
# windows() # on windows
# quartz() #MAC
# X11() #Ubuntu

# pdf (here ("output","Fig3_216spp.pdf"), width=15,height=7,family="serif")

plot(panel2)

# dev.off()

# ---------------------------------------------------------- #
## these results consider only Oryzomyialia (~200 spp)


source ("./R/Interpretation_Fig2_and_statistics_oryzomyialia.R")

# get the panel with the relationship between MPD and DISPARITY

# to save run
# pdf (here ("output","Fig2_oryz.pdf"), width=6,height=7,family="serif")

plot(panel)

# dev.off()

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
source ("./R/Interpretation_Fig3_oryzomyialia.R")

# if you want to save then run
# windows() # on windows
# quartz() #MAC
# X11() #Ubuntu

# pdf (here ("output","Fig3_oryz.pdf"), width=15,height=7,family="serif")

plot(panel2)

# dev.off()

# ---------------------------------------------------------- #
# these results are produced by simulation traits across all 413 spp


source ("./R/Interpretation_Fig2_and_statistics_all_spp.R")

# get the panel with the relationship between MPD and DISPARITY

# to save run
# pdf (here ("output","Fig2_all.pdf"), width=6,height=7,family="serif")

plot(panel)

# dev.off()

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
source ("./R/Interpretation_Fig3_all_spp.R")

# if you want to save then run
# windows() # on windows
# quartz() #MAC
# X11() #Ubuntu

# pdf (here ("output","Fig3_all.pdf"), width=15,height=7,family="serif")

plot(panel2)

# dev.off()

