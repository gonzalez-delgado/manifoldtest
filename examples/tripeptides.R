# ================================================================
# tripeptides.R
# -----------------------------------------------------------------
# Independence testing on tripeptide backbone dihedral angles, using
# torus.dcov.test() from the manifoldtest package.
#
# For each central amino acid, (Left, Central, Right) tripeptides
# are formed from the surrounding residues. Backbone dihedral angles
# (Phi, Psi) are rescaled from (-pi, pi] to [0,1)^2.
#
#   * Control test (H0): central-residue structure of two different
#     tripeptides are expected to be independent.
#   * H1 test: within the same tripeptide, the central residue's structure
#     is tested against its left and right neighbors' structure, where
#     dependence is expected.
#
# Each test is run twice via torus.dcov.test(): once with
# center_outward = FALSE (raw tangent-space coordinates) and once with
# center_outward = TRUE (center-outward rank transform), to compare the
# two testing approaches.
#
# Usage:
#   cd examples && Rscript tripeptides.R
# ================================================================

library(manifoldtest)
library(ggplot2)
library(ggpubr)

set_path <- 'data/tripeptides/'
results_dir <- 'results'
alpha <- 0.05
R_perm <- 500
n_top <- 50  # number of most abundant tripeptides used for the control test
seed <- 2024  # global random seed
set.seed(seed)

dir.create(results_dir, showWarnings = FALSE)
h0_file <- file.path(results_dir, 'tripeptide_pvalues_H0.Rda')
h1_file <- file.path(results_dir, 'tripeptide_pvalues_H1.Rda')

amino_acid_list <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS",
                      "ILE", "LEU", "LYS", "MET", "PHE", "SER", "THR", "TRP",
                      "TYR", "VAL")

trip_list <- expand.grid(amino_acid_list, amino_acid_list, amino_acid_list)
colnames(trip_list) <- c('Left', 'Central', 'Right')

# ================================================================
# Helpers
# ================================================================

#' Rows of a central-residue data set matching a given (Left, Central, Right)
subset_triplet <- function(data_central, left, central, right) {
    data_central[data_central$Res1 == left &
                 data_central$Res2 == central &
                 data_central$Res3 == right, ]
}

#' Rescale a residue's (Phi, Psi) dihedral angles from (-pi, pi] to [0,1)^2
to_torus <- function(data_trip, residue) {
    cols <- c(paste0("Phi_res_", residue), paste0("Psi_res_", residue))
    (data_trip[, cols] + pi) / (2 * pi)
}

#' Run torus.dcov.test() in both modes and return a named p-value vector
test_both_modes <- function(X, Y) {

    raw <- torus.dcov.test(X=X,
                          Y=Y,
                          R=R_perm,
                          seed=seed,
                          center_outward=FALSE,
                          verbose=FALSE)

    co  <- torus.dcov.test(X=X,
                          Y=Y,
                          R=R_perm,
                          seed=seed,
                          center_outward=TRUE,
                          verbose=FALSE)

    c(raw = raw$dcov_test$p.value, co = co$dcov_test$p.value)
}

if (file.exists(h0_file) && file.exists(h1_file)) {

    # ================================================================
    # Cached results found: load them instead of recomputing, so plots
    # can be reformatted without rerunning the (expensive) tests.
    # ================================================================

    cat(sprintf("Found cached results in '%s'; loading instead of recomputing.\n", results_dir))
    pvalues_control <- readRDS(h0_file)
    pvalues <- readRDS(h1_file)

} else {

    # ================================================================
    # 1. Load per-central-residue tripeptide data
    # ================================================================

    all_data <- list()
    for (central_name in amino_acid_list) {

        data_central <- get(load(file.path(set_path, paste0(central_name, "_angles.RData"))))

        # Remove tripeptides with PRO or GLY as a neighbor
        data_central <- data_central[data_central$Res1 != 'PRO' & data_central$Res3 != 'PRO' &
                                      data_central$Res1 != 'GLY' & data_central$Res3 != 'GLY', ]

        # Keep trans conformations only
        data_central <- data_central[data_central$CIS.0..TRANS.1. == 1, ]
        all_data[[central_name]] <- data_central
    }

    # ================================================================
    # 2. Count samples per tripeptide and filter to those with at least 100 samples
    # ================================================================

    pvalues <- data.frame(trip_list,
                          n = NA,
                          pv_left = NA,
                          pv_left_co = NA,
                          pv_right = NA,
                          pv_right_co = NA)

    for (i_trip in seq_len(nrow(trip_list))) {

        central_name <- as.character(trip_list$Central[i_trip])
        data_trip <- subset_triplet(all_data[[central_name]],
                                    as.character(trip_list$Left[i_trip]),
                                    central_name,
                                    as.character(trip_list$Right[i_trip]))

        pvalues$n[i_trip] <- nrow(data_trip)
    }

    keep <- pvalues$n > 100
    pvalues <- pvalues[keep, ]
    trip_list <- trip_list[keep, ]

    # ================================================================
    # 3. Control test (H0): (Left', Central', Right') independent of (Left, Central, Right)
    # ================================================================

    pvalues <- pvalues[order(pvalues$n, decreasing = TRUE), ]
    first_pv <- pvalues[1:n_top, ]
    trip_pairs <- t(combn(1:n_top, 2))

    pvalues_control <- data.frame(trip_pairs,
                                  pv=NA,
                                  pv_co=NA)

    colnames(pvalues_control)[1:2] <- c('trip_A', 'trip_B')
    n_sub <- min(min(first_pv$n), 1000)

    for (trip_pair in seq_len(nrow(trip_pairs))) {

        trip_A <- first_pv[trip_pairs[trip_pair, 1], ]
        trip_B <- first_pv[trip_pairs[trip_pair, 2], ]

        data_trip_A <- subset_triplet(data_central=all_data[[as.character(trip_A$Central)]],
                                      left=as.character(trip_A$Left),
                                      central=as.character(trip_A$Central),
                                      right=as.character(trip_A$Right))

        data_trip_B <- subset_triplet(data_central=all_data[[as.character(trip_B$Central)]],
                                      left=as.character(trip_B$Left),
                                      central=as.character(trip_B$Central),
                                      right=as.character(trip_B$Right))

        A_central <- to_torus(data_trip_A, 2)[sample(nrow(data_trip_A), n_sub), ]
        B_central <- to_torus(data_trip_B, 2)[sample(nrow(data_trip_B), n_sub), ]

        pvalues_control[trip_pair, c("pv", "pv_co")] <- test_both_modes(A_central, B_central)

        if (trip_pair %% 100 == 0) cat(sprintf("Control test: %d / %d pairs done\n", trip_pair, nrow(trip_pairs)))
    }

    saveRDS(pvalues_control, h0_file)

    # ================================================================
    # 4. H1 test: central residue vs. left/right neighbor, within the same tripeptide
    # ================================================================

    for (i_trip in seq_len(nrow(trip_list))) {

        central_name <- as.character(trip_list$Central[i_trip])
        data_trip <- subset_triplet(data_central=all_data[[central_name]],
                                    left=as.character(trip_list$Left[i_trip]),
                                    central=as.character(trip_list$Central[i_trip]),
                                    right=as.character(trip_list$Right[i_trip]))

        central <- to_torus(data_trip, 2)
        left    <- to_torus(data_trip, 1)
        right   <- to_torus(data_trip, 3)

        pvalues[i_trip, c("pv_left", "pv_left_co")]   <- test_both_modes(central, left)
        pvalues[i_trip, c("pv_right", "pv_right_co")] <- test_both_modes(central, right)

        if (i_trip %% 50 == 0) cat(sprintf("H1 test: %d / %d triplets done\n", i_trip, nrow(trip_list)))
    }

    saveRDS(pvalues, h1_file)
}

# ================================================================
# 5. Summary and plots: raw tangent-space test vs. center-outward test
# ================================================================

cat(sprintf("Rejection rate at alpha = %.2f:\n", alpha))
cat(sprintf("  left,  raw:             %.1f%%\n", mean(pvalues$pv_left    < alpha, na.rm = TRUE) * 100))
cat(sprintf("  left,  center-outward:  %.1f%%\n", mean(pvalues$pv_left_co < alpha, na.rm = TRUE) * 100))
cat(sprintf("  right, raw:             %.1f%%\n", mean(pvalues$pv_right    < alpha, na.rm = TRUE) * 100))
cat(sprintf("  right, center-outward:  %.1f%%\n", mean(pvalues$pv_right_co < alpha, na.rm = TRUE) * 100))

#' Stack a (raw, center-outward) p-value pair into a long data frame for plotting
to_long <- function(pv_raw, pv_co) {
    rbind(data.frame(pv = pv_raw, method = "raw"),
          data.frame(pv = pv_co,  method = "center-outward"))
}

ecdf_plot <- function(data_long, subtitle, title) {

    ggplot(data_long, aes(x = pv, color = method)) +
        stat_ecdf(linewidth = 1) +
        geom_abline(linetype = 'dashed', color = 'darkblue') +
        labs(x = 'p-value', y = 'ECDF', color = 'Test', subtitle = subtitle) +
        ggtitle(title) +
        theme_bw()
}

p1 <- ecdf_plot(to_long(pvalues_control$pv, pvalues_control$pv_co),
                'Testing independence between structures of different tripeptides',
                'p-value distribution under the null')

p2 <- ecdf_plot(to_long(pvalues$pv_left, pvalues$pv_left_co),
                'Testing independence between the structures of central and left amino-acids',
                'p-value distribution under a fixed alternative')

p3 <- ecdf_plot(to_long(pvalues$pv_right, pvalues$pv_right_co),
                'Testing independence between the structures of central and right amino-acids',
                'p-value distribution under a fixed alternative')

p_all <- ggarrange(p1, p2, p3, ncol = 3, labels = c('(a)', '(b)', '(c)'), common.legend = TRUE)
ggsave(file.path(results_dir, 'tripeptide_pvalues.pdf'), p_all, width = 12, height = 4)
