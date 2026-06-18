## run_sa_fusion.R
## Analysis script for our sex-autosome (SA) fusion rate test in mammal clades 
## with neo-sex chromosomes (Herpestidae, Phyllostomidae), adapted from 
## Chien & Blackmon's Scarabaeoidea analysis
## paper DOI: 10.1093/jeb/voag025
## data DOI: 10.5281/zenodo.19457362

## Written (ideally) for CLI use (cluster runs), usage:
## > Rscript run_sa_fusion.R [family] [n_trees] [nsim] [n_cores] [method]
##     family:  herp, phyllo, or (default) both
##     n_trees: number of posterior trees to use (default 1000)
##     nsim:    stochastic maps per tree (default 100)
##     n_cores: number of CPU cores (default 1) used by mclapply
##     method:  pertree (estimate Q per tree, default) or fixedQ (Q estimated 
##              once and reused across trees)

## e.g. > Rscript run_sa_fusion.R both 15 2 (test run)
##      > Rscript run_sa_fusion.R both 1000 100 (full thing)

## paths
PROJ       <- "./"
FUSION_DIR <- "/input"
OUTDIR     <- "/output/sa_fusion_results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

source("sa_fusion_functions.R")

## args
args    <- commandArgs(trailingOnly = TRUE)
FAMILY  <- if (length(args) >= 1) args[1] else "both"
N_TREES <- if (length(args) >= 2) as.integer(args[2]) else 1000
NSIM    <- if (length(args) >= 3) as.integer(args[3]) else 100
N_CORES <- if (length(args) >= 4) as.integer(args[4]) else
           max(1L, as.integer(Sys.getenv("N_CORES", "1")))
METHOD  <- if (length(args) >= 5) args[5] else Sys.getenv("METHOD", "pertree")
RNGkind("L'Ecuyer-CMRG")   # independent streams for mclapply
set.seed(42)

cat(sprintf("Config: family=%s  n_trees=%d  nsim=%d  n_cores=%d  method=%s\n",
            FAMILY, N_TREES, NSIM, N_CORES, METHOD))

## taxonomy fixes (as always...)
TREE_RENAME <- list(
  Herpestidae = c(
    "Herpestes_brachyurus"   = "Urva_brachyura",
    "Herpestes_edwardsii"    = "Urva_edwardsii",
    "Herpestes_fuscus"       = "Urva_fusca",
    "Herpestes_semitorquatus"= "Urva_semitorquata",
    "Herpestes_smithii"      = "Urva_smithii",
    "Herpestes_urva"         = "Urva_urva",
    "Herpestes_vitticollis"  = "Urva_vitticollis",
    "Herpestes_javanicus"    = "Urva_javanica",
    "Herpestes_naso"         = "Xenogale_naso"
  ),
  Phyllostomidae = c(
    "Dermanura_aztecus"      = "Dermanura_azteca",
    "Dermanura_cinereus"     = "Dermanura_cinerea",
    "Dermanura_glaucus"      = "Dermanura_glauca",
    "Dermanura_gnomus"       = "Dermanura_gnoma",
    "Dermanura_rosenbergii"  = "Dermanura_rosenbergi",
    "Dermanura_toltecus"     = "Dermanura_tolteca",
    "Lophostoma_aequatorialis"= "Lophostoma_occidentale",
    "Lophostoma_silvicolum"  = "Lophostoma_silvicola",
    "Mimon_crenulatum"       = "Gardnerycteris_crenulata",
    "Mimon_koepckeae"        = "Gardnerycteris_koepckeae",
    "Vampyressa_bidens"      = "Vampyriscus_bidens",
    "Vampyressa_brocki"      = "Vampyriscus_brocki",
    "Vampyressa_nymphaea"    = "Vampyriscus_nymphaeus",
    "Diaemus_youngi"         = "Diaemus_youngii"
  )
)
## taxonomic reconciliation from Mammal Diversity Database v1.11
## data DOI: 10.5281/zenodo.4139722
MDD_UPDATE <- character(0)

## tree tip reconciliation
rename_tree_tips <- function(trees, map) {
  if (length(map) == 0) return(trees)
  for (i in seq_along(trees)) {
    tl  <- trees[[i]]$tip.label
    safe <- names(map)[!(unname(map) %in% tl)]
    dropped <- setdiff(names(map), safe)
    if (length(dropped) && i == 1)
      warning("rename_tree_tips: skipped rename(s) that would duplicate an ",
              "existing tip: ", paste(dropped, collapse = ", "), call. = FALSE)
    hit <- tl %in% safe
    tl[hit] <- map[tl[hit]]
    trees[[i]]$tip.label <- tl
  }
  trees
}

## XY -> 0, single fusion (XAfusion/YAfusion) -> 1, XAYAfusion -> 2
sex_to_nfus <- function(sx) {
  sx <- trimws(sx)
  out <- rep(1L, length(sx))
  out[sx == "XY"] <- 0L
  out[grepl("XAYA", sx, ignore.case = TRUE)] <- 2L
  out
}

load_family_df <- function(csv) {
  d <- read.csv(csv, stringsAsFactors = FALSE, row.names = 1)
  nm <- gsub(" ", "_", d$Scientific.Name)
  hit <- nm %in% names(MDD_UPDATE); nm[hit] <- MDD_UPDATE[nm[hit]]
  d <- data.frame(
    name        = nm,
    auto.hap    = suppressWarnings(as.integer(round(as.numeric(d$auto.hap)))),
    XYorvariant = suppressWarnings(as.integer(d$XYorvariant)),
    nfus        = sex_to_nfus(d$sex),
    stringsAsFactors = FALSE
  )
  d <- d[stats::complete.cases(d), ]
  d <- d[!duplicated(d$name), ]
  d
}

run_and_save <- function(tag, csv, nex, rename = NULL) {
  cat(sprintf("\n========== %s ==========\n", tag))
  df <- load_family_df(csv)
  cat(sprintf("known-karyo species: %d  (0-fusion=%d, 1-fusion=%d, 2-fusion=%d)  auto.hap %d-%d  K=%d\n",
              nrow(df), sum(df$nfus == 0), sum(df$nfus == 1), sum(df$nfus == 2),
              min(df$auto.hap), max(df$auto.hap), max(df$nfus) + 1L))

  trees <- read.nexus(nex)
  if (inherits(trees, "phylo")) trees <- structure(list(trees), class = "multiPhylo")
  if (!is.null(rename)) trees <- rename_tree_tips(trees, rename)
  trees <- rename_tree_tips(trees, MDD_UPDATE)

  ## report taxon matching after tip reconciliation (like phylolm does)
  miss <- setdiff(df$name, trees[[1]]$tip.label)
  cat(sprintf("taxon match: %d/%d species found in trees%s\n",
              sum(df$name %in% trees[[1]]$tip.label), nrow(df),
              if (length(miss)) paste0("  | MISSING: ", paste(miss, collapse = ", ")) else ""))
  cat(sprintf("trees available: %d  (using %d)\n", length(trees), min(N_TREES, length(trees))))

  res <- run_family(trees, df, nsim = NSIM, n_trees = N_TREES, n_cores = N_CORES,
                    method = METHOD)

  write.csv(data.frame(obs = res$obs), file.path(OUTDIR, paste0(tag, "_obs.csv")),
            row.names = FALSE)
  write.csv(data.frame(exp = res$exp), file.path(OUTDIR, paste0(tag, "_exp.csv")),
            row.names = FALSE)

  ## summary + plot
  summ <- summarise_overlap(res$obs, res$exp)
  summ <- cbind(family = tag, method = METHOD, n_trees = res$n_trees_run,
                nsim = NSIM, summ)
  write.csv(summ, file.path(OUTDIR, paste0(tag, "_summary.csv")), row.names = FALSE)
  plot_overlap(res$obs, res$exp, main = tag,
               file = file.path(OUTDIR, paste0(tag, "_density.pdf")))

  cat("\n--- ", tag, " result ---\n", sep = "")
  print(summ, row.names = FALSE)
  summ
}

jobs <- list(
  herp   = list(tag = "Herpestidae",
                csv = file.path(PROJ, "chromeplus/data_known_karyo/tetrapodtraits_herp_known_karyo.csv"),
                nex = file.path(PROJ, "multiphylo/herpoutput.nex")),
  phyllo = list(tag = "Phyllostomidae",
                csv = file.path(PROJ, "chromeplus/data_known_karyo/tetrapodtraits_phyllo_known_karyo.csv"),
                nex = file.path(PROJ, "multiphylo/phyllostomoutput.nex"))
)

todo <- if (FAMILY == "both") names(jobs) else FAMILY
all_summ <- list()
for (j in todo) all_summ[[j]] <- with(jobs[[j]], run_and_save(tag, csv, nex, TREE_RENAME[[tag]]))

if (length(all_summ)) {
  combined <- do.call(rbind, all_summ)
  write.csv(combined, file.path(OUTDIR, "ALL_summary.csv"), row.names = FALSE)
  cat("\n================ COMBINED ================\n")
  print(combined, row.names = FALSE)
  cat("\nOutputs written to: ", OUTDIR, "\n")
}
