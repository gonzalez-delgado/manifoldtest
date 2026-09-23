# ================================================================
# trinucleotides.R
# -----------------------------------------------------------------
# Independence testing on RNA trinucleotide torsion angles, using
# manifold.dcov.test() from the manifoldtest package.
#
# For each trinucleotide type present in the data, the sugar-pucker
# angles X = (nu0, nu1, nu2) are tested for independence against
# Y = (nu3, nu4). Angles are in degrees in [0, 360) and
# are rescaled to [0,1)^d.
#
# Each test is run twice via manifold.dcov.test(): once with
# center_outward = FALSE (raw tangent-space coordinates) and once with
# center_outward = TRUE (center-outward rank transform), to compare the
# two testing approaches. p-values are Holm-Bonferroni corrected across
# trinucleotides, then written out as a LaTeX table (one column per
# trinucleotide, one row per method).
#
# Usage:
#   cd examples && Rscript trinucleotides.R
# ================================================================

# ---- Install missing packages ----
for (pkg in c("parallel", "pbapply", "reticulate")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("manifoldtest", quietly = TRUE)) {
    # Not on CRAN: install from this local source tree (run from examples/, so
    # the package root is one level up)
    install.packages("..", repos = NULL, type = "source")
}

library(manifoldtest)
library(parallel)
library(pbapply)
pboptions(type = "timer")  # show the progress bar under Rscript too (off by default when non-interactive)

data_path <- 'data/nucleotides.npy'
results_dir <- 'results'
alpha <- 0.05
R_perm <- 500
seed <- 2024  # global random seed
n_cores <- 3 #max(1, detectCores() - 1, na.rm = TRUE)  # workers for the parallel loop
set.seed(seed)

dir.create(results_dir, showWarnings = FALSE)
pv_file <- file.path(results_dir, 'trinucleotide_pvalues.Rda')

# ================================================================
# Helpers
# ================================================================

#' Rescale angle columns (in degrees, [0, 360)) to [0,1)^d
to_torus <- function(data_trinuc, cols) {
    (data_trinuc[, cols] %% 360) / 360
}

#' Run manifold.dcov.test() in both modes and return a named p-value vector
test_both_modes <- function(X, Y) {

    raw <- manifold.dcov.test(X=X,
                             Y=Y,
                             manifold="torus",
                             R=R_perm,
                             seed=seed,
                             center_outward=FALSE,
                             verbose=FALSE)

    co  <- manifold.dcov.test(X=X,
                             Y=Y,
                             manifold="torus",
                             R=R_perm,
                             seed=seed,
                             center_outward=TRUE,
                             verbose=FALSE)

    c(raw = raw$dcov_test$p.value, co = co$dcov_test$p.value)
}

if (file.exists(pv_file)) {

    # ================================================================
    # Cached results found: load them instead of recomputing, so plots
    # can be reformatted without rerunning the (expensive) tests.
    # ================================================================

    cat(sprintf("Found cached results in '%s'; loading instead of recomputing.\n", results_dir))
    pvalues <- readRDS(pv_file)

} else {

    # ================================================================
    # 1. Load trinucleotide data: one data frame per trinucleotide type
    # ================================================================

    cat("Loading trinucleotide data...\n")
    np <- reticulate::import("numpy")
    rna_data <- np$load(data_path, allow_pickle = TRUE)
    rna_data <- rna_data[[1]]

    trinuc_list <- names(rna_data)

    # ================================================================
    # 2. Test X = (nu0, nu1, nu2) vs Y = (nu3, nu4), for every trinucleotide
    # ================================================================

    cat(sprintf("Performing independence test on %d cores (%d trinucleotides)...\n",
                n_cores, length(trinuc_list)))

    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, library(manifoldtest))
    clusterExport(cl, c("rna_data", "to_torus", "test_both_modes", "R_perm", "seed"))
    clusterSetRNGStream(cl, seed)

    results <- pblapply(trinuc_list, function(trinuc) {

        data_trinuc <- rna_data[[trinuc]]
        X <- to_torus(data_trinuc, c('nu0', 'nu1', 'nu2'))
        Y <- to_torus(data_trinuc, c('nu3', 'nu4'))
        c(n = nrow(data_trinuc), test_both_modes(X, Y))

    }, cl = cl)
    stopCluster(cl)

    res_mat <- do.call(rbind, results)
    pvalues <- data.frame(trinucleotide = trinuc_list,
                             n = res_mat[, "n"],
                             pv = res_mat[, "raw"],
                             pv_co = res_mat[, "co"])

    saveRDS(pvalues, pv_file)
}

# ================================================================
# 3. Multiplicity correction and summary
# ================================================================

pvalues$pv_adj    <- p.adjust(pvalues$pv,    method = "holm")
pvalues$pv_co_adj <- p.adjust(pvalues$pv_co, method = "holm")

cat(sprintf("Rejection rate at alpha = %.2f (Holm-adjusted):\n", alpha))
cat(sprintf("  raw:             %.1f%%\n", mean(pvalues$pv_adj    < alpha, na.rm = TRUE) * 100))
cat(sprintf("  center-outward:  %.1f%%\n", mean(pvalues$pv_co_adj < alpha, na.rm = TRUE) * 100))

# ================================================================
# 4. LaTeX table: one column per trinucleotide, one row per method
# ================================================================

fmt_p <- function(p) formatC(p, format = "f", digits = 4)
col_spec <- paste0("l", strrep("c", nrow(pvalues)))

tex_lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    paste("Method &", paste(pvalues$trinucleotide, collapse = " & "), "\\\\"),
    "\\midrule",
    paste("raw &", paste(fmt_p(pvalues$pv_adj), collapse = " & "), "\\\\"),
    paste("center-outward &", paste(fmt_p(pvalues$pv_co_adj), collapse = " & "), "\\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\caption{Holm-Bonferroni adjusted $p$-values for independence between sugar-pucker and backbone torsion angles, by trinucleotide.}",
    "\\label{tab:trinucleotide-pvalues}",
    "\\end{table}"
)

tex_file <- file.path(results_dir, 'trinucleotide_pvalues.txt')
writeLines(tex_lines, tex_file)
cat(sprintf("LaTeX table written to '%s'\n", tex_file))
