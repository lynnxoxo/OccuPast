# uncertainty.R
#
# Uncertainty helpers for OccuPast.
#
# Design:
# - within-replicate uncertainty is estimated from site-by-bin data
# - default estimator is site bootstrap
# - summary statistic per bin is the mean across sites
# - between-replicate uncertainty is pooled across temporal replicates
#
# Expected site_bin schema (from analyze_curve.R):
#   site_id, horizon_bin, bin_start, bin_end, value / value_norm, ...
#
# Expected replicate-estimate schema for pooling:
#   replicate_id, horizon_bin, bin_start, bin_end, estimate, var_hat

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   magrittr
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows group_by summarise left_join n_distinct
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort
#   @importFrom magrittr %>%

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.require_site_bin_uncertainty <- function(site_bin, value_col) {
  validate_required_fields(
    site_bin,
    c("site_id", "horizon_bin", "bin_start", "bin_end", value_col),
    "site_bin"
  )
  invisible(TRUE)
}

.require_replicate_estimates <- function(replicate_estimates) {
  validate_required_fields(
    replicate_estimates,
    c("replicate_id", "horizon_bin", "bin_start", "bin_end", "estimate", "var_hat"),
    "replicate_estimates"
  )
  invisible(TRUE)
}

.site_bin_to_matrix <- function(site_bin, value_col = "value_norm") {
  site_bin <- tibble::as_tibble(site_bin)
  if (exists(".standardize_bin_columns", mode = "function")) {
    site_bin <- .standardize_bin_columns(site_bin, arg = "site_bin")
  }
  .require_site_bin_uncertainty(site_bin, value_col)

  bin_tbl <- unique(site_bin[, c("horizon_bin", "bin_start", "bin_end"), drop = FALSE])
  bin_tbl <- bin_tbl[order(bin_tbl$horizon_bin), , drop = FALSE]

  site_ids <- sort(unique(site_bin$site_id))
  bin_ids <- bin_tbl$horizon_bin

  m <- matrix(
    0,
    nrow = length(site_ids),
    ncol = length(bin_ids),
    dimnames = list(site_ids, as.character(bin_ids))
  )

  row_idx <- match(site_bin$site_id, site_ids)
  col_idx <- match(site_bin$horizon_bin, bin_ids)

  vals <- site_bin[[value_col]]
  vals[is.na(vals)] <- 0

  for (i in seq_len(nrow(site_bin))) {
    m[row_idx[i], col_idx[i]] <- vals[i]
  }

  list(
    mat = m,
    site_ids = site_ids,
    bin_tbl = tibble::as_tibble(bin_tbl)
  )
}

.empty_uncertainty_tbl <- function() {
  tibble::tibble(
    horizon_bin = numeric(),
    bin_start = numeric(),
    bin_end = numeric(),
    estimate = numeric(),
    var_hat = numeric(),
    se = numeric()
  )
}

# -------------------------------------------------------------------------
# Within-replicate variance: site bootstrap
# -------------------------------------------------------------------------

#' Estimate within-replicate variance using a site bootstrap
#'
#' @param site_bin Site-by-bin table.
#' @param site_col Site identifier column. Currently must be "site_id".
#' @param bin_col Bin identifier column. Currently must be "horizon_bin".
#' @param value_col Value column to analyze, usually "value" or "value_norm".
#' @param B Number of bootstrap resamples.
#' @param seed Optional RNG seed.
#'
#' @return A list with:
#'   - data: one row per bin with estimate, var_hat, se
#'   - diagnostics: bootstrap diagnostics
within_variance_site_bootstrap <- function(site_bin,
                                           site_col = "site_id",
                                           bin_col = "horizon_bin",
                                           value_col = "value_norm",
                                           B = 200L,
                                           seed = NULL) {
  if (!identical(site_col, "site_id")) {
    rlang::abort("This implementation expects `site_col = \"site_id\"`.")
  }
  if (!identical(bin_col, "horizon_bin")) {
    rlang::abort("This implementation expects `bin_col = \"horizon_bin\"`.")
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 2) {
    rlang::abort("`B` must be a single integer >= 2.")
  }

  sbm <- .site_bin_to_matrix(site_bin, value_col = value_col)
  mat <- sbm$mat
  bin_tbl <- sbm$bin_tbl

  n_sites <- nrow(mat)
  n_bins <- ncol(mat)

  if (n_sites == 0 || n_bins == 0) {
    return(list(
      data = .empty_uncertainty_tbl(),
      diagnostics = list(
        n_sites = n_sites,
        n_bins = n_bins,
        B = as.integer(B)
      )
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Observed estimate: equal-site-weighted mean in each bin.
  estimate <- colMeans(mat)

  boot_means <- matrix(NA_real_, nrow = B, ncol = n_bins)

  for (b in seq_len(B)) {
    idx <- sample.int(n_sites, size = n_sites, replace = TRUE)
    boot_means[b, ] <- colMeans(mat[idx, , drop = FALSE])
  }

  # Bootstrap variance: Var_hat(theta) = var(theta*_1, ..., theta*_B).
  var_hat <- apply(boot_means, 2, stats::var)
  se <- sqrt(var_hat)

  out <- tibble::tibble(
    horizon_bin = bin_tbl$horizon_bin,
    bin_start = bin_tbl$bin_start,
    bin_end = bin_tbl$bin_end,
    estimate = as.numeric(estimate),
    var_hat = as.numeric(var_hat),
    se = as.numeric(se)
  )

  list(
    data = out,
    diagnostics = list(
      n_sites = n_sites,
      n_bins = n_bins,
      B = as.integer(B),
      summary_statistic = "mean_across_sites"
    )
  )
}

# -------------------------------------------------------------------------
# Within-replicate variance: site jackknife
# -------------------------------------------------------------------------

#' Estimate within-replicate variance using a leave-one-site-out jackknife
#'
#' @param site_bin Site-by-bin table.
#' @param site_col Site identifier column. Currently must be "site_id".
#' @param bin_col Bin identifier column. Currently must be "horizon_bin".
#' @param value_col Value column to analyze, usually "value" or "value_norm".
#'
#' @return A list with:
#'   - data: one row per bin with estimate, var_hat, se
#'   - diagnostics: jackknife diagnostics
within_variance_site_jackknife <- function(site_bin,
                                           site_col = "site_id",
                                           bin_col = "horizon_bin",
                                           value_col = "value_norm") {
  if (!identical(site_col, "site_id")) {
    rlang::abort("This implementation expects `site_col = \"site_id\"`.")
  }
  if (!identical(bin_col, "horizon_bin")) {
    rlang::abort("This implementation expects `bin_col = \"horizon_bin\"`.")
  }

  sbm <- .site_bin_to_matrix(site_bin, value_col = value_col)
  mat <- sbm$mat
  bin_tbl <- sbm$bin_tbl

  n_sites <- nrow(mat)
  n_bins <- ncol(mat)

  if (n_sites == 0 || n_bins == 0) {
    return(list(
      data = .empty_uncertainty_tbl(),
      diagnostics = list(
        n_sites = n_sites,
        n_bins = n_bins
      )
    ))
  }

  if (n_sites < 2) {
    out <- tibble::tibble(
      horizon_bin = bin_tbl$horizon_bin,
      bin_start = bin_tbl$bin_start,
      bin_end = bin_tbl$bin_end,
      estimate = as.numeric(colMeans(mat)),
      var_hat = NA_real_,
      se = NA_real_
    )

    return(list(
      data = out,
      diagnostics = list(
        n_sites = n_sites,
        n_bins = n_bins,
        note = "Jackknife variance undefined for fewer than 2 sites."
      )
    ))
  }

  estimate <- colMeans(mat)

  loo_means <- matrix(NA_real_, nrow = n_sites, ncol = n_bins)
  for (i in seq_len(n_sites)) {
    loo_means[i, ] <- colMeans(mat[-i, , drop = FALSE])
  }

  loo_bar <- colMeans(loo_means)
  # Delete-1 jackknife variance: ((n-1)/n) * sum((theta_-i - mean(theta_-i))^2).
  var_hat <- ((n_sites - 1) / n_sites) * colSums((loo_means - matrix(loo_bar, nrow = n_sites, ncol = n_bins, byrow = TRUE))^2)
  se <- sqrt(var_hat)

  out <- tibble::tibble(
    horizon_bin = bin_tbl$horizon_bin,
    bin_start = bin_tbl$bin_start,
    bin_end = bin_tbl$bin_end,
    estimate = as.numeric(estimate),
    var_hat = as.numeric(var_hat),
    se = as.numeric(se)
  )

  list(
    data = out,
    diagnostics = list(
      n_sites = n_sites,
      n_bins = n_bins,
      summary_statistic = "mean_across_sites"
    )
  )
}

# -------------------------------------------------------------------------
# Pooling across temporal replicates
# -------------------------------------------------------------------------

#' Pool within- and between-replicate uncertainty
#'
#' @param replicate_estimates Tibble with:
#'   replicate_id, horizon_bin, bin_start, bin_end, estimate, var_hat
#' @param conf_level Confidence level.
#' @param use_t If TRUE, use Rubin-style t critical values when possible.
#'
#' @return Tibble with pooled estimate and uncertainty components per bin.
pool_replicate_uncertainty <- function(replicate_estimates,
                                       conf_level = 0.95,
                                       use_t = TRUE) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort("`conf_level` must be a single number between 0 and 1.")
  }

  replicate_estimates <- tibble::as_tibble(replicate_estimates)
  if (exists(".standardize_bin_columns", mode = "function")) {
    replicate_estimates <- .standardize_bin_columns(replicate_estimates, arg = "replicate_estimates")
  }
  .require_replicate_estimates(replicate_estimates)

  pooled <- replicate_estimates %>%
    dplyr::group_by(horizon_bin, bin_start, bin_end) %>%
    dplyr::summarise(
      M = dplyr::n_distinct(replicate_id),
      estimate_mean = mean(estimate, na.rm = TRUE),
      W = mean(var_hat, na.rm = TRUE),
      B = if (dplyr::n() > 1) stats::var(estimate, na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    dplyr::rename(estimate = estimate_mean)

  pooled$W[!is.finite(pooled$W)] <- NA_real_
  pooled$B[!is.finite(pooled$B)] <- 0

  # Rubin total variance: T = W + (1 + 1/M)B.
  pooled$T <- pooled$W + (1 + 1 / pooled$M) * pooled$B
  pooled$se_within <- sqrt(pooled$W)
  pooled$se_between <- sqrt(pooled$B)
  pooled$se_total <- sqrt(pooled$T)

  if (use_t) {
    pooled$df <- Inf

    use_df <- pooled$M > 1 & is.finite(pooled$W) & is.finite(pooled$B) & pooled$B > 0
    pooled$df[use_df] <- (pooled$M[use_df] - 1) *
      (1 + pooled$W[use_df] / ((1 + 1 / pooled$M[use_df]) * pooled$B[use_df]))^2

    crit <- stats::qt((1 + conf_level) / 2, df = pooled$df)
    crit[!is.finite(crit)] <- stats::qnorm((1 + conf_level) / 2)
  } else {
    pooled$df <- Inf
    crit <- rep(stats::qnorm((1 + conf_level) / 2), nrow(pooled))
  }

  pooled$crit <- crit
  pooled$conf_level <- conf_level
  pooled$lower <- pooled$estimate - pooled$crit * pooled$se_total
  pooled$upper <- pooled$estimate + pooled$crit * pooled$se_total

  tibble::as_tibble(pooled)
}

# -------------------------------------------------------------------------
# Uncertainty decomposition helper
# -------------------------------------------------------------------------

#' Summarize uncertainty components in long form
#'
#' @param pooled_result Output from pool_replicate_uncertainty().
#'
#' @return Long-form tibble with within, between, and total variance/SE by bin.
summarize_uncertainty_components <- function(pooled_result) {
  pooled_result <- tibble::as_tibble(pooled_result)
  if (exists(".standardize_bin_columns", mode = "function")) {
    pooled_result <- .standardize_bin_columns(pooled_result, arg = "pooled_result")
  }

  validate_required_fields(
    pooled_result,
    c("horizon_bin", "bin_start", "bin_end", "W", "B", "T", "se_within", "se_between", "se_total"),
    "pooled_result"
  )

  out <- dplyr::bind_rows(
    tibble::tibble(
      horizon_bin = pooled_result$horizon_bin,
      bin_start = pooled_result$bin_start,
      bin_end = pooled_result$bin_end,
      component = "within",
      variance = pooled_result$W,
      se = pooled_result$se_within
    ),
    tibble::tibble(
      horizon_bin = pooled_result$horizon_bin,
      bin_start = pooled_result$bin_start,
      bin_end = pooled_result$bin_end,
      component = "between",
      variance = pooled_result$B,
      se = pooled_result$se_between
    ),
    tibble::tibble(
      horizon_bin = pooled_result$horizon_bin,
      bin_start = pooled_result$bin_start,
      bin_end = pooled_result$bin_end,
      component = "total",
      variance = pooled_result$T,
      se = pooled_result$se_total
    )
  )

  tibble::as_tibble(out)
}


# -------------------------------------------------------------------------
# Legacy aliases for earlier "bucket" terminology in argument names.
# -------------------------------------------------------------------------

within_variance_site_bootstrap_bucket <- within_variance_site_bootstrap
within_variance_site_jackknife_bucket <- within_variance_site_jackknife
