## sa_fusion_functions.R
## Functions for our sex-autosome (SA) fusion rate test, adapted from 
## Chien & Blackmon's Scarabaeoidea analysis (SA_fusion_prob.R & functions.R) 
## for our mammal clades with neo-sex chromosomes (Herpestidae, Phyllostomidae)

library(ape)
library(phytools)
library(coda)

if (!exists("FUSION_DIR")) {
  FUSION_DIR <- "./input/"
}
.oldwd <- getwd()
setwd(FUSION_DIR)
source("functions.R")
setwd(.oldwd)

## call Pfsa2() function, except the null is the symmetric no-drive expectation
if (!exists("Pfsa")) 
  Pfsa <- function(Da, scs) Pfsa2(Da = Da, scs = scs, mud = 0.5)

## rate matrix building (counts single SAfusion (i.e. XA, YAfusion) and 
## dual SAfusion (i.e. XAYAfusion) as separate levels with their own states)

## rate classes: 
## fission                   = 2f+1
## auto-auto fusion          = 2f+2
## SA fusion (same auto.hap) = 2K+1
## SA fusion reversal        = 2K+2
## XAYAfusion automatically counts as two state transitions (two events)
## because 0->1->2 transition has to happen sequentially

build_model_matrix <- function(lo, hi, K = 2) {
  L <- hi - lo + 1
  S <- K * L
  m <- matrix(0L, S, S)
  ## within-level
  for (f in 0:(K - 1)) {
    off <- f * L
    asc_cls  <- 2L * f + 1L
    desc_cls <- 2L * f + 2L
    if (L >= 2) for (j in 1:(L - 1)) {
      m[off + j,     off + j + 1] <- asc_cls # fission
      m[off + j + 1, off + j    ] <- desc_cls # AAfusion
    }
  }
  ## SA fusions between consecutive levels
  sa_cls  <- 2L * K + 1L
  rev_cls <- 2L * K + 2L
  if (K >= 2) for (f in 0:(K - 2)) {
    off <- f * L
    for (j in 1:L) {
      m[off + j,     off + L + j] <- sa_cls # SAfusion
      m[off + L + j, off + j    ] <- rev_cls # reversal
    }
  }
  rownames(m) <- colnames(m) <- as.character(1:S)
  m
}

## build matrix as df with proper columns
build_tip_matrix <- function(df, lo, hi, K) {
  L <- hi - lo + 1
  S <- K * L
  x <- matrix(0, nrow = nrow(df), ncol = S)
  rownames(x) <- df$name
  colnames(x) <- as.character(1:S)
  for (k in seq_len(nrow(df))) {
    j   <- df$auto.hap[k] - lo + 1 # autosome index within level
    col <- df$nfus[k] * L + j # offset by level
    x[k, col] <- 1
  }
  x
}

## estimate the rate matrix Q once for fast "Qfix/fixedQ" estimation runs 
## (estimating Q per tree takes its sweet time in Phyllostomidae but leads to 
## more robust results)

fit_Q <- function(tree, x, model, pi = "fitzjohn") {
  keep <- intersect(tree$tip.label, rownames(x))
  tree <- keep.tip(tree, keep)
  tree <- force.ultrametric(tree, method = "extend")
  xx   <- x[tree$tip.label, , drop = FALSE]
  fit  <- fitMk(tree, xx, model = model, pi = pi)
  Q    <- unclass(as.Qmatrix(fit))
  ## reorder rows/cols to match state order
  ord  <- match(colnames(x), colnames(Q))
  if (anyNA(ord)) stop("fit_Q: state label mismatch between Q and x")
  Q    <- Q[ord, ord, drop = FALSE]
  dimnames(Q) <- list(colnames(x), colnames(x))
  Q
}

## single-tree observed/expected SA-fusion proportions analysis
## Qfix = NULL -> bespoke Q estimation on single tree
## Qfix = matrix -> utilize fit_Q
## either returns a list or NULL on failure
run_one_tree <- function(tree, x, model, lo, hi, K = 2, nsim = 50,
                         rejmax = 1e6, rejint = 1e5, Qfix = NULL) {
  L <- hi - lo + 1
  S <- K * L

  ## tree & data congruence
  keep <- intersect(tree$tip.label, rownames(x))
  if (length(keep) < 4) return(NULL)
  tree <- keep.tip(tree, keep)
  tree <- force.ultrametric(tree, method = "extend")
  xx   <- x[tree$tip.label, , drop = FALSE]

  ## stochastic maps (Qfix vs not Qfix)
  maps <- tryCatch(
    if (is.null(Qfix)) {
      make.simmap2(tree, x = xx, model = model, pi = "fitzjohn",
                   nsim = nsim, rejmax = rejmax, rejint = rejint, monitor = FALSE)
    } else {
      make.simmap2(tree, x = xx, model = model, Q = Qfix, pi = "fitzjohn",
                   nsim = nsim, rejmax = rejmax, rejint = rejint, monitor = FALSE)
    },
    error = function(e) { message("\nsimmap failed: ", conditionMessage(e)); NULL })
  if (is.null(maps)) return(NULL)
  if (inherits(maps, "phylo")) maps <- list(maps)
  class(maps) <- "multiPhylo"
  counts <- describe.simmap2(maps)$count
  if (is.null(dim(counts))) counts <- matrix(counts, nrow = 1,
                                             dimnames = list(NULL, names(counts)))

  ## transition column names are "from,to"
  ## AA fusion: descending intralevel dysploidy
  aa_cols <- unlist(lapply(0:(K - 1), function(f) {
    off <- f * L
    if (L >= 2) paste(off + (2:L), off + (1:(L - 1)), sep = ",") else character(0)
  }))
  ## SA fusion: interlevel f -> f+1 at matching auto.hap index (counts each
  ## X- and Y-autosome fusion separately)
  sa_cols <- if (K >= 2) unlist(lapply(0:(K - 2), function(f) {
    off <- f * L
    paste(off + (1:L), off + L + (1:L), sep = ",")
  })) else character(0)

  grab <- function(cn, cols) {
    hit <- intersect(cols, colnames(cn))
    if (!length(hit)) return(rep(0, nrow(cn)))
    rowSums(cn[, hit, drop = FALSE])
  }
  SA <- grab(counts, sa_cols)
  AA <- grab(counts, aa_cols)
  obs <- SA / (SA + AA)

  ## time-weighted mean of Pfsa over observed states
  autohap_of_state <- function(s) lo + ((s - 1) %% L)
  pfsa_state <- vapply(1:S,
                       function(s) Pfsa(Da = 2 * autohap_of_state(s), scs = "XY"),
                       numeric(1))
  expv <- numeric(length(maps))
  for (i in seq_along(maps)) {
    tm <- describe.simmap(maps[[i]])$times["prop", ]
    tm <- tm[names(tm) != "total"]
    w  <- setNames(numeric(S), as.character(1:S))
    w[names(tm)] <- tm
    expv[i] <- sum(pfsa_state * w)
  }

  list(obs = obs, exp = expv, ntip = length(keep))
}

## multiPhylo across-family analysis (run_one_tree but 1000 times)
run_family <- function(trees, df, nsim = 50, n_trees = NULL,
                       rejmax = 1e6, rejint = 1e5, verbose = TRUE,
                       n_cores = 1, method = c("pertree", "fixedQ")) {
  method <- match.arg(method)
  lo <- min(df$auto.hap); hi <- max(df$auto.hap)
  K  <- max(df$nfus) + 1L
  model <- build_model_matrix(lo, hi, K)
  x     <- build_tip_matrix(df, lo, hi, K)

  if (is.null(n_trees)) n_trees <- length(trees)
  n_trees <- min(n_trees, length(trees))

  ## using pre-esimated Q if specified
  Qfix <- NULL
  if (method == "fixedQ") {
    if (verbose) cat("\nfitting Q once by ML (fixedQ method)...\n")
    Qfix <- fit_Q(trees[[1]], x, model, pi = "fitzjohn")
    if (verbose) cat("\nQ fitted.\n")
  }

  one <- function(t) {
    res <- tryCatch(run_one_tree(trees[[t]], x, model, lo, hi, K = K,
                                 nsim = nsim, rejmax = rejmax, rejint = rejint,
                                 Qfix = Qfix),
                    error = function(e) NULL)
    if (verbose) {
      if (is.null(res)) cat(sprintf("\ntree %d/%d - skipped\n", t, n_trees))
      else cat(sprintf("\ntree %d/%d - obs=%.3f exp=%.3f (ntip=%d)\n", t, n_trees,
                       mean(res$obs, na.rm = TRUE), mean(res$exp, na.rm = TRUE),
                       res$ntip))
    }
    res
  }

  if (n_cores > 1) {
    results <- parallel::mclapply(seq_len(n_trees), one,
                                  mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    results <- lapply(seq_len(n_trees), one)
  }

  obs_list <- lapply(results, function(r) if (is.null(r)) NULL else r$obs)
  exp_list <- lapply(results, function(r) if (is.null(r)) NULL else r$exp)
  list(obs = unlist(obs_list), exp = unlist(exp_list),
       lo = lo, hi = hi, K = K, model = model,
       n_trees_run = sum(!vapply(results, is.null, logical(1))))
}

## plotting code in style of beetle paper
summarise_overlap <- function(obs, exp) {
  obs <- obs[is.finite(obs)]; exp <- exp[is.finite(exp)]
  oh <- HPDinterval(as.mcmc(obs))
  eh <- HPDinterval(as.mcmc(exp))
  no_overlap <- (eh[2] < oh[1]) || (oh[2] < eh[1])
  data.frame(
    obs_mean = mean(obs), obs_lo = oh[1], obs_hi = oh[2],
    exp_mean = mean(exp), exp_lo = eh[1], exp_hi = eh[2],
    no_overlap = no_overlap,
    obs_greater = mean(obs) > mean(exp)
  )
}

plot_overlap <- function(obs, exp, main = "", file = NULL) {
  obs <- obs[is.finite(obs)]; exp <- exp[is.finite(exp)]
  cols   <- viridis::viridis(2, begin = 0.5, alpha = 0.55)
  hpdcol <- viridis::viridis(2, begin = 0.5)
  if (!is.null(file)) pdf(file, width = 6, height = 6)
  de <- density(exp, bw = 0.02); do <- density(obs, bw = 0.02)
  plot(de, xlim = range(0, de$x, do$x, 0.7), ylim = c(-1, max(de$y, do$y) * 1.1),
       main = main, xlab = "Proportion of sex-autosome fusion", lwd = 1)
  polygon(de, col = cols[1]); polygon(do, col = cols[2])
  legend("topright", fill = cols, legend = c("Expected (null)", "Observed"),
         bty = "n", cex = 0.9)
  he <- HPDinterval(as.mcmc(exp)); ho <- HPDinterval(as.mcmc(obs))
  lines(y = c(-0.4, -0.4), x = he[1:2], lwd = 5, col = hpdcol[1])
  lines(y = c(-0.8, -0.8), x = ho[1:2], lwd = 5, col = hpdcol[2])
  if (!is.null(file)) dev.off()
}
