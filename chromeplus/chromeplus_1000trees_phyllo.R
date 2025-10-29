library(ape)
library(phytools)
library(dplyr)
library(chromePlus)
library(diversitree)
library(doMC)

tetrapodtraits.herp.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_herp_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)
tetrapodtraits.phyllo.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_phyllo_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)
tetrapodtraits.sorex.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_sorex_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)

herpestidae.tree <- read.nexus(file = "individual_trees/herpestidae_tree.nex")
phyllo.tree <- read.nexus(file = "individual_trees/phyllostom_tree.nex")
sorex.tree <- read.nexus(file = "individual_trees/sorex_tree.nex")

herp.multiphylo <- read.nexus(file = "multiphylo/herp_multiphylo.nex")
phyllo.multiphylo <- read.nexus(file = "multiphylo/phyllostom_multiphylo.nex")
sorex.multiphylo <- read.nexus(file = "multiphylo/sorex_multiphylo.nex")

class(herp.multiphylo) <- "multiPhylo"
class(phyllo.multiphylo) <- "multiPhylo"
class(sorex.multiphylo) <- "multiPhylo"

## modify following to fit the family being analyzed as necessary
phyllo.tree.known.karyo <- keep.tip(phyllo.tree, tetrapodtraits.phyllo.known.karyo$Scientific.Name)
chromeplus.phyllo.df <- tetrapodtraits.phyllo.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.phyllo.matrix <- datatoMatrix(chromeplus.phyllo.df, hyper = T)

phyllo.keep.tip <- function(tree) {
  keep.tip(tree, tetrapodtraits.phyllo.known.karyo$Scientific.Name)
  return(tree)
}

phyllo.multiphylo <- lapply(phyllo.multiphylo, keep.tip, tetrapodtraits.phyllo.known.karyo$Scientific.Name)
phyllo.multiphylo <- lapply(phyllo.multiphylo, force.ultrametric)

## this will (hopefully) cycle through all 1000 trees in each multiphylo
registerDoMC(128)


#results <- list()
iter <- 10000

results <- foreach(i = 1:1000, .errorhandling="remove", .packages = c("chromePlus", "diversitree"), .combine = rbind) %dopar% {
  # make the initial likelihood function
  lk.mk <- make.mkn(phyllo.multiphylo[[i]], states = chromeplus.phyllo.matrix, 
                    k = ncol(chromeplus.phyllo.matrix), strict=F, control=list(method="ode"))
  # constrain to a biologically realistic model of chrom evolution
  con.lk.mk<-constrainMkn(data = chromeplus.phyllo.matrix, lik = lk.mk, hyper = T, 
                          polyploidy = F, verbose = F, 
			  constrain = list(drop.demi = T, drop.poly = T))

  prior <- make.prior.exponential(1/2)

  # estimate for w
  temp <- mcmc(con.lk.mk, x.init=runif(6,0,10), w=1, prior=prior, nsteps=100)
  w <- diff(sapply(temp[2:7], quantile, c(.05, .95)))
  # run the MCMC
 # results[[i]]
 par_results <- diversitree::mcmc(con.lk.mk, x.init = colMeans(temp)[2:7], 
                       w = w, prior=prior, nsteps = iter)
 par_results <- par_results[1001:nrow(par_results),]
 return(par_results)
}

saveRDS(results, "phyllo1000.rds")
write.table(results, "phyllo1000.txt")
