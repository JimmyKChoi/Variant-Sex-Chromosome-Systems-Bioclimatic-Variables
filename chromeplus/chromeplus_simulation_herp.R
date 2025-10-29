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
herp.tree.known.karyo <- keep.tip(herpestidae.tree, tetrapodtraits.herp.known.karyo$Scientific.Name)
chromeplus.herp.df <- tetrapodtraits.herp.known.karyo %>% select(Scientific.Name, auto.hap, XYorvariant)
chromeplus.herp.matrix <- datatoMatrix(chromeplus.herp.df, hyper = T)

iter <- 1000
herp.tree.known.karyo <- force.ultrametric(herp.tree.known.karyo)
herp.lk.mk <- make.mkn(tree = herp.tree.known.karyo, states = chromeplus.herp.matrix, k = ncol(chromeplus.herp.matrix), strict=F, control=list(method="ode"))
herp.con.lk.mk<-constrainMkn(data=chromeplus.herp.matrix, lik=herp.lk.mk, hyper=T, 
                             polyploidy=F, verbose=F, 
                             constrain=list(drop.demi=T, drop.poly=T))
temp <- mcmc(herp.con.lk.mk, x.init=runif(6,0,10), w=1, nsteps=iter/10)
w <- diff(sapply(temp[1:6], quantile, c(.05, .95)))
write.table(w, "herp_w.txt")
herp.mcmc <- mcmc(herp.con.lk.mk, x.init=runif(6,0,10), w=w, nsteps=iter)
write.table(herp.mcmc, "herp_mcmc.txt")

