library(ape)
library(phytools)
library(dplyr)
library(chromePlus)
library(diversitree)

## other two karyotype datasets also included for convenience, script only includes run for Herpestidae
tetrapodtraits.herp.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_herp_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)
tetrapodtraits.phyllo.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_phyllo_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)
tetrapodtraits.sorex.known.karyo <- read.csv(file = "data_known_karyo/tetrapodtraits_sorex_known_karyo.csv", stringsAsFactors = FALSE, row.names = 1)

## phylogenies
herpestidae.tree <- read.nexus(file = "individual_trees/herpestidae_tree.nex")
phyllostom.tree <- read.nexus(file = "individual_trees/phyllostom_tree.nex")
sorex.tree <- read.nexus(file = "individual_trees/sorex_tree.nex")

## modify following to fit the family being analyzed as necessary
phyllo.tree.known.karyo <- keep.tip(phyllostom.tree, tetrapodtraits.phyllo.known.karyo$Scientific.Name)
chromeplus.phyllo.df <- tetrapodtraits.phyllo.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.phyllo.matrix <- datatoMatrix(chromeplus.phyllo.df, hyper = T)

iter <- 1000
phyllo.tree.known.karyo <- force.ultrametric(phyllo.tree.known.karyo)
phyllo.lk.mk <- make.mkn(tree = phyllo.tree.known.karyo, states = chromeplus.phyllo.matrix, k = ncol(chromeplus.phyllo.matrix), strict=F, control=list(method="ode"))
phyllo.con.lk.mk<-constrainMkn(data=chromeplus.phyllo.matrix, lik=phyllo.lk.mk, hyper=T, 
                             polyploidy=F, verbose=F, 
                             constrain=list(drop.demi=T, drop.poly=T))
temp <- mcmc(phyllo.con.lk.mk, x.init=runif(6,0,10), w=1, nsteps=iter/10)
w <- diff(sapply(temp[1:6], quantile, c(.05, .95)))
write.table(w, "phyllo_w.txt")
phyllo.mcmc <- mcmc(phyllo.con.lk.mk, x.init=runif(6,0,10), w=w, nsteps=iter)
write.table(phyllo.mcmc, "phyllo_mcmc.txt")

