# load needed packages

require("geomorph")
require ("SYNCSA")

#### simulate some data
# species composition (10 communities of 10 spp, from a pool of 100)
pool<-100
ncomm<-10
# random sample from the pool
comm_data <- replicate (ncomm,
                        sample (seq (1,pool, 1),ncomm))

#  matrix with composition
composition_data <- matrix(rbinom(ncomm*pool,1,0.2),
                           nrow=ncomm,
                           ncol=pool)
colnames(composition_data)<-paste ("s",seq(1,nrow(traits_pool)),sep="")
rownames(composition_data)<- paste ("c",seq(1,nrow(composition_data)),sep="")
# pool traits, two traits taken from a normal distribution
traits_pool <- cbind (t1=rnorm(100,0),
                      t2=rnorm(100,0))

rownames(traits_pool) <- paste ("s",seq(1,nrow(traits_pool)),sep="")

# organize data
organized_data <- organize.syncsa (comm=as.matrix(composition_data),
                                   traits=traits_pool)


# calculate Rao from syncsa
rao_observado <-rao.diversity(organized_data$community,
                              traits=organized_data$traits)$FunRao

# calculate disparity from geomorph
sp_to_get <- lapply (seq(1,nrow(organized_data$community)), function (i) 
  which(organized_data$community[i,]>0))
# trait subset
trait_per_site_subset <- lapply (sp_to_get, function (i) 
  traits_pool[i,] 
)

# what a M**F** function
geom_res <- rbind(morphol.disparity(trait_per_site_subset[[1]]~1),
      morphol.disparity(trait_per_site_subset[[2]]~1),
      morphol.disparity(trait_per_site_subset[[3]]~1),
      morphol.disparity(trait_per_site_subset[[4]]~1),
      morphol.disparity(trait_per_site_subset[[5]]~1),
      morphol.disparity(trait_per_site_subset[[6]]~1),
      morphol.disparity(trait_per_site_subset[[7]]~1),
      morphol.disparity(trait_per_site_subset[[8]]~1),
      morphol.disparity(trait_per_site_subset[[9]]~1),
      morphol.disparity(trait_per_site_subset[[10]]~1))
# are they correlated?
(res<- data.frame(rao = rao_observado,
                 morphol.disp = geom_res))

cor(res$rao,res$morphol.disp) # No !

