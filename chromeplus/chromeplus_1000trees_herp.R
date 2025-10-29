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
herp.tree.known.karyo <- keep.tip(herpestidae.tree, tetrapodtraits.herp.known.karyo$Scientific.Name)
chromeplus.herp.df <- tetrapodtraits.herp.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.herp.matrix <- datatoMatrix(chromeplus.herp.df, hyper = T)

herp.keep.tip <- function(tree) {
  keep.tip(tree, tetrapodtraits.herp.known.karyo$Scientific.Name)
  return(tree)
}

herp.multiphylo <- lapply(herp.multiphylo, keep.tip, tetrapodtraits.herp.known.karyo$Scientific.Name)
herp.multiphylo <- lapply(herp.multiphylo, force.ultrametric)

## this will (hopefully) cycle through all 1000 trees in each multiphylo
registerDoMC(32)


#results <- list()
iter <- 10000

results <- foreach(i = 1:1000, .packages = c("chromePlus", "diversitree"), .combine = rbind) %dopar% {
  # make the initial likelihood function
  lk.mk <- make.mkn(herp.multiphylo[[i]], states = chromeplus.herp.matrix, 
                    k = ncol(chromeplus.herp.matrix), strict=F, control=list(method="ode"))
  # constrain to a biologically realistic model of chrom evolution
  con.lk.mk<-constrainMkn(data = chromeplus.herp.matrix, lik = lk.mk, hyper = T, 
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

saveRDS(results, "herp1000.rds")
write.table(results, "herp1000.txt")
