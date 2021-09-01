
# function designed to use observed and simulated data to
# calculate SES disparitu

funcao_disparidade <- function (occ, traits,n_iterations) {
  data<- as.matrix(traits) # shape.v.2d.means
  colnames(data)<-seq(1,ncol(data))
  
  ## organize data using a function from SYNCSA package
  organized_obs_data <- organize.syncsa (comm=as.matrix(occ),traits=data)
  # calculate and extract Rao entropy
  rao_observado <-rao.diversity(organized_obs_data$community,traits=organized_obs_data$traits)$FunRao
  
  # arguments to calculate the SES
  # transform into data frame at the end
  RESULTADOS.RAO <- lapply(seq(1,n_iterations), function (i) {
    row.names(organized_obs_data$traits) <- sample(row.names(organized_obs_data$traits)) # shuffling
    org.perm <- organize.syncsa(comm=organized_obs_data$community,traits=organized_obs_data$traits)
    rao.perm <- rao.diversity(org.perm$community,traits=org.perm$traits)
    rao_aleat <-rao.perm$FunRao; # return
    rao_aleat
  })
  # average disparity
  mean.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,mean)
  # standard deviation of disparity
  sd.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,sd)
  # data frame with results (and SES already calculated)
  DF_RESULTADOS <- data.frame (Observado = rao_observado, med_nulo=mean.rao.NULO,sd_nulo=sd.rao.NULO,
                               SES = ((rao_observado - mean.rao.NULO)/sd.rao.NULO))
  
  return (DF_RESULTADOS)
  
}
 
# ---------------------------------
# function to calculate mpd (from Swenson 2014)
mpd.function <- function (x,dist.tree) {
  
  # get the names of species present in the comms
  com.names <- names (x[which(x>0)])
  
  # calcuate mpd
  mean (as.dist (dist.tree[which(rownames (dist.tree) %in% com.names),
                           which(colnames (dist.tree) %in% com.names)]))
  
}


# function to calculate mpd after shuffling the phylogeny
mpd.shuff <- function (tree,my.sample.matrix) {
  
  shuff.tree <- tipShuffle (tree)
  apply(my.sample.matrix,1, mpd.function,
        dist.tree=cophenetic (shuff.tree))
  
}