This is a project for determining potential relationships between various bioclimatic variables and mammal sex chromosome systems that differ from the "classic" XX/XY system, through phylogenetic comparative methods using various R packages.

Other tangentially related stuff includes ancestral state reconstruction, etc.

**updated herpestidae stuff.Rmd** - Original (and most complete) script for this project, analyses of species in family Herpestidae with variant sex chromosome systems, plus data visualizations

**batproject.Rmd** - Analyses of species in family Phyllostomidae with variant sex chromosome system

**asr stuff.Rmd** - Ancestral state reconstruction script for variant Herpestidae species

**sex chromosome project updated.Rmd** - Master script used for analysis and testing, new & updated version of analysis scripts integrating all three mammal clades of interest, using data from Moura _et al._ (2024). Includes single-tree analyses using _phylolm_ (phylo. log regression) and multiPhylo analyses using _sensiPhy_. Currently being split into three separate scripts, with the goal of being able to run these with no modification.

**shapefile_stuff.Rmd** - A testbench script, mostly meant as an experiment in plotting shapefiles with ggplot to build Figure 2 in our paper. Allows for visualization of distributions for our mammal families of interest.

<ins>Ready-to-run scripts (run these in sequence!)</ins>

**initial_phyloglm_analyses.Rmd** - The first in a series of scripts acting as polished versions of the master script, meant to allow easy reproduction of results. This one focuses on the initial analyses performed with the consensus tree using _phylolm_.

**sensiphy_analyses.Rmd** - The second in a series of scripts acting as polished versions of the master script. This one focuses on the _sensiPhy_ analyses performing phylo. log regression over sets of 1000 trees from the posterior distribution for each mammal family of interest.

to do: 
- make more figures, including with plotTree.barplot!
- separate master script into individual sections for analysis and polish these up to allow for easy reproduction of results
- write this paper already...

Citations:
K. T. David, Global gradients in the distribution of animal polyploids. Proc. Natl. Acad. Sci. U.S.A. 119, e2214070119 (2022). Code modified from original scripts.

Moura, M. R. et al. A phylogeny-informed characterisation of global tetrapod traits addresses data gaps and biases. PLOS Biology 22, e3002658 (2024). Ecological & species data collected from here.


