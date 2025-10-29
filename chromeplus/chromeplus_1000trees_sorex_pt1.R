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
phyllostom.tree <- read.nexus(file = "individual_trees/phyllostom_tree.nex")
sorex.tree <- read.nexus(file = "individual_trees/sorex_tree.nex")

herp.multiphylo <- read.nexus(file = "multiphylo/herp_multiphylo.nex")
phyllostom.multiphylo <- read.nexus(file = "multiphylo/phyllostom_multiphylo.nex")
sorex.multiphylo <- read.nexus(file = "multiphylo/sorex_multiphylo.nex")

class(herp.multiphylo) <- "multiPhylo"
class(phyllostom.multiphylo) <- "multiPhylo"
class(sorex.multiphylo) <- "multiPhylo"

## modify following to fit the family being analyzed as necessary
sorex.tree.known.karyo <- keep.tip(sorex.tree, tetrapodtraits.sorex.known.karyo$Scientific.Name)
chromeplus.sorex.df <- tetrapodtraits.sorex.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.sorex.matrix <- datatoMatrix(chromeplus.sorex.df, hyper = T)

sorex.keep.tip <- function(tree) {
  keep.tip(tree, tetrapodtraits.sorex.known.karyo$Scientific.Name)
  return(tree)
}

sorex.multiphylo <- lapply(sorex.multiphylo, keep.tip, tetrapodtraits.sorex.known.karyo$Scientific.Name)
sorex.multiphylo <- lapply(sorex.multiphylo, force.ultrametric)

## this will (hopefully) cycle through all 1000 trees in each multiphylo
registerDoMC(256)


#results <- list()
iter <- 10000

results <- foreach(i = 1:500, .errorhandling="remove", .packages = c("chromePlus", "diversitree"), .combine = rbind) %dopar% {
  # make the initial likelihood function
  lk.mk <- make.mkn(sorex.multiphylo[[i]], states = chromeplus.sorex.matrix, 
                    k = ncol(chromeplus.sorex.matrix), strict=F, control=list(method="ode"))
  # constrain to a biologically realistic model of chrom evolution
  con.lk.mk<-constrainMkn(data = chromeplus.sorex.matrix, lik = lk.mk, hyper = T, 
                          polyploidy = F, verbose = F, 
			  constrain = list(drop.demi = T, drop.poly = T))

  prior <- make.prior.exponential(1/2)

  # estimate for w
  temp <- mcmc(con.lk.mk, x.init=runif(6,0,10), w=1, prior=prior, nsteps=iter/200)
  w <- diff(sapply(temp[2:7], quantile, c(.05, .95)))
  # run the MCMC
 # results[[i]]
 par_results <- diversitree::mcmc(con.lk.mk, x.init = colMeans(temp)[2:7], 
                       w = w, prior=prior, nsteps = iter)
 par_results <- par_results[1001:nrow(par_results),]
 return(par_results)
}

saveRDS(results, "sorex1000_pt1.rds")
write.table(results, "sorex1000_pt1.txt")
