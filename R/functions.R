
## funcao serve para pegar os atributos observados e simulados e calcular rao nulo e sd nulo em uma filogenia

funcao_disparidade <- function (occ, traits,n_iterations) {
  data<- as.matrix(traits) # shape.v.2d.means
  colnames(data)<-seq(1,ncol(data))
  
  ## organizar os dados observados
  organized_obs_data <- organize.syncsa (comm=as.matrix(occ),traits=data)
  # calcular Rao e extrai-lo
  rao_observado <-rao.diversity(organized_obs_data$community,traits=organized_obs_data$traits)$FunRao
  
  ### SES nulo para ol observado
  
  RESULTADOS.RAO <- lapply(seq(1,n_iterations), function (i) {
    row.names(organized_obs_data$traits) <- sample(row.names(organized_obs_data$traits))
    org.perm <- organize.syncsa(comm=organized_obs_data$community,traits=organized_obs_data$traits)
    rao.perm <- rao.diversity(org.perm$community,traits=org.perm$traits)
    rao_aleat <-rao.perm$FunRao; rao_aleat
  })
  
  mean.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,mean)
  sd.rao.NULO <- apply(do.call(cbind,RESULTADOS.RAO),1,sd)
  
  DF_RESULTADOS <- data.frame (Observado = rao_observado, med_nulo=mean.rao.NULO,sd_nulo=sd.rao.NULO,
                               SES = ((rao_observado - mean.rao.NULO)/sd.rao.NULO))
  
  return (DF_RESULTADOS)
  
}
