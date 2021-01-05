# ------------------------# 

# Run scripts of results

# these  results consider the 216 species in trait, comm, and phy datasets

source ("./R/Interpretation.R")

# get the image

panel

# get model estimates 
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

## these results consider only Oryzomyialia (~200 spp)

source ("./R/Interpretation_oryz.R")


# get the image

panel

## these results consider the complete phylogeny of 413 species to simulate traits
## after that we pruned the simulated trait dataset to match community / phy datasets

source ("./R/Interpretation_ALL.R")