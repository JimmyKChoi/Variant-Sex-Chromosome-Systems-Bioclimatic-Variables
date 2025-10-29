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
sorex.tree.known.karyo <- keep.tip(sorex.tree, tetrapodtraits.sorex.known.karyo$Scientific.Name)
chromeplus.sorex.df <- tetrapodtraits.sorex.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.sorex.matrix <- datatoMatrix(chromeplus.sorex.df, hyper = T)

iter <- 1000
sorex.tree.known.karyo <- force.ultrametric(sorex.tree.known.karyo)
sorex.lk.mk <- make.mkn(tree = sorex.tree.known.karyo, states = chromeplus.sorex.matrix, k = ncol(chromeplus.sorex.matrix), strict=F, control=list(method="ode"))
sorex.con.lk.mk<-constrainMkn(data=chromeplus.sorex.matrix, lik=sorex.lk.mk, hyper=T, 
                             polyploidy=F, verbose=F, 
                             constrain=list(drop.demi=T, drop.poly=T))
temp <- mcmc(sorex.con.lk.mk, x.init=runif(6,0,10), w=1, nsteps=iter/10)
w <- diff(sapply(temp[1:6], quantile, c(.05, .95)))
write.table(w, "sorex_w.txt")
sorex.mcmc <- mcmc(sorex.con.lk.mk, x.init=runif(6,0,10), w=w, nsteps=iter)
write.table(sorex.mcmc, "sorex_mcmc.txt")

