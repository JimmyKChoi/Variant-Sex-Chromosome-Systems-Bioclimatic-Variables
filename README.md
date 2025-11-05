This is a project for determining potential relationships between various bioclimatic variables and mammal sex chromosome systems that differ from the "classic" XX/XY system, through phylogenetic comparative methods using various R packages.

Other tangentially related stuff includes ancestral state reconstruction, etc.

<ins>Archived scripts</ins>

**updated herpestidae stuff.Rmd** - Original script for this project, analyses of species in family Herpestidae with variant sex chromosome systems, plus data visualizations

**batproject.Rmd** - Analyses of species in family Phyllostomidae with variant sex chromosome system

**asr stuff.Rmd** - Ancestral state reconstruction script for variant Herpestidae species

**sex chromosome project updated.Rmd** - Master script used for analysis and testing, new & updated version of analysis scripts integrating all three mammal clades of interest, using data from Moura _et al._ (2024). Includes single-tree analyses using _phylolm_ (phylo. log regression) and multiPhylo analyses using _sensiPhy_. Currently being split into three separate scripts, with the goal of being able to run these with no modification.

**shapefile_stuff.Rmd** - A testbench script, mostly meant as an experiment in plotting shapefiles with ggplot to build Figure 2 in our paper. Allows for visualization of distributions for our mammal families of interest.

<ins>Input and output data</ins>

**sexreferencedfupdated.csv** - A dataset containing every mammal species considered for this study, their associated sex chromosome systems, whether these sex chromosomes are known or imputed, citation if known, and miscellaneous comments. Also **Supplemental Table 1** for our paper. Table caption: "Comprehensive list of mammal species included in our analysis, sampled from three mammalian families. Of 34 total species from Herpestidae included in our analysis (excluding Urva auropunctata), 9 species were coded as possessing a Y-autosome fusion. Of 413 total species from Soricidae (excluding Sorex kozlovi), 10 species were coded as possessing an X-autosome fusion. Of 205 total species from Phyllostomidae, 15 species were coded as possessing an X-autosome fusion and 14 species were coded as possessing both an X- and Y-autosome fusion, summing to 29 variant species total."

**sexreferencedf.csv** - Simplified version of above dataset that gets used in our analyses.

**multiphylo** - Nexus files for each of our three mammal families included in this study, sets of 1000 credible trees for use in **sensiphy_analyses.Rmd**.

**output** - Output files from code, essentially subsets of TetrapodTraits 1.0.0 containing only our species of interest for each of our three mammal families. Some extra variables from other datasets (e.g. "chroms.csv") have been added.

<ins>Ready-to-run scripts (run these in sequence!)</ins>

**initial_phyloglm_analyses.Rmd** - The first in a series of scripts acting as polished versions of the master script, meant to allow easy reproduction of results. This one focuses on the initial analyses performed with the consensus tree using _phylolm_.

**sensiphy_analyses.Rmd** - The second in a series of scripts acting as polished versions of the master script. This one focuses on the _sensiPhy_ analyses performing phylo. log regression over sets of 1000 trees from the posterior distribution for each mammal family of interest.

**chromeplus_simulation_herp.R** & **chromeplus_1000trees_x.R** - The third/fourth in a series of scripts acting as polished versions of the master script. These are the _chromePlus_ simulations, with the 1000trees one performing the simulation (you guessed it) over the set of 1000 trees we use for _sensiPhy_.

to do: 
- submit this paper already...

Citations:
Blackmon, H. "Meiotic drive shapes rates of karyotype evolution in mammals." Evolution 73, 511-523 (2019). Chromosome morphology data ("chroms.csv") collected from here.

David, K. T. "Global gradients in the distribution of animal polyploids." Proc. Natl. Acad. Sci. U.S.A. 119, e2214070119 (2022). Code modified from original scripts.

Hughes, J. J. et al. "The role of conflict in the formation and maintenance of variant sex chromosome systems in mammals." J. Hered. 115, 601-624 (2024). Sex chromosome coding data collected from here.

Moura, M. R. et al. "A phylogeny-informed characterisation of global tetrapod traits addresses data gaps and biases." PLOS Biology 22, e3002658 (2024). Ecological & species data collected from here.

Upham, N. S. et al. "Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation." PLOS Biology 17, e3000494 (2019). Mammal tree collected from here.




