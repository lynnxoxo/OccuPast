# uncertainty.R
#
# Uncertainty helpers for the merged occupancy package.
#
# Design:
# - within-replicate uncertainty is estimated from site-by-bucket data
# - default estimator is site bootstrap
# - summary statistic per bucket is the mean across sites
# - between-replicate uncertainty is pooled across temporal replicates
#
# Expected site_bucket schema (from analyze_curve.R):
#   site_id, horizon_bucket, bucket_start, bucket_end, value / value_norm, ...
#
# Expected replicate-estimate schema for pooling:
#   replicate_id, horizon_bucket, bucket_start, bucket_end, estimate, var_hat

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows group_by summarise left_join n_distinct
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.require_site_bucket_uncertainty <- function(site_bucket, value_col) {
  validate_required_fields(
    site_bucket,
    c("site_id", "horizon_bucket", "bucket_start", "bucket_end", value_col),
    "site_bucket"
  )
  invisible(TRUE)
}

.require_replicate_estimates <- function(replicate_estimates) {
  validate_required_fields(
    replicate_estimates,
    c("replicate_id", "horizon_bucket", "bucket_start", "bucket_end", "estimate", "var_hat"),
    "replicate_estimates"
  )
  invisible(TRUE)
}

.site_bucket_to_matrix <- function(site_bucket, value_col = "value_norm") {
  site_bucket <- tibble::as_tibble(site_bucket)
  .require_site_bucket_uncertainty(site_bucket, value_col)

  bucket_tbl <- unique(site_bucket[, c("horizon_bucket", "bucket_start", "bucket_end"), drop = FALSE])
  bucket_tbl <- bucket_tbl[order(bucket_tbl$horizon_bucket), , drop = FALSE]

  site_ids <- sort(unique(site_bucket$site_id))
  bucket_ids <- bucket_tbl$horizon_bucket

  m <- matrix(
    0,
    nrow = length(site_ids),
    ncol = length(bucket_ids),
    dimnames = list(site_ids, as.character(bucket_ids))
  )

  row_idx <- match(site_bucket$site_id, site_ids)
  col_idx <- match(site_bucket$horizon_bucket, bucket_ids)

  vals <- site_bucket[[value_col]]
  vals[is.na(vals)] <- 0

  for (i in seq_len(nrow(site_bucket))) {
    m[row_idx[i], col_idx[i]] <- vals[i]
  }

  list(
    mat = m,
    site_ids = site_ids,
    bucket_tbl = tibble::as_tibble(bucket_tbl)
  )
}

.empty_uncertainty_tbl <- function() {
  tibble::tibble(
    horizon_bucket = numeric(),
    bucket_start = numeric(),
    bucket_end = numeric(),
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
#' @param site_bucket Site-by-bucket table.
#' @param site_col Site identifier column. Currently must be "site_id".
#' @param bucket_col Bucket identifier column. Currently must be "horizon_bucket".
#' @param value_col Value column to analyze, usually "value" or "value_norm".
#' @param B Number of bootstrap resamples.
#' @param seed Optional RNG seed.
#'
#' @return A list with:
#'   - data: one row per bucket with estimate, var_hat, se
#'   - diagnostics: bootstrap diagnostics
within_variance_site_bootstrap <- function(site_bucket,
                                           site_col = "site_id",
                                           bucket_col = "horizon_bucket",
                                           value_col = "value_norm",
                                           B = 200L,
                                           seed = NULL) {
  if (!identical(site_col, "site_id")) {
    rlang::abort("This implementation expects `site_col = \"site_id\"`.")
  }
  if (!identical(bucket_col, "horizon_bucket")) {
    rlang::abort("This implementation expects `bucket_col = \"horizon_bucket\"`.")
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 2) {
    rlang::abort("`B` must be a single integer >= 2.")
  }

  sbm <- .site_bucket_to_matrix(site_bucket, value_col = value_col)
  mat <- sbm$mat
  bucket_tbl <- sbm$bucket_tbl

  n_sites <- nrow(mat)
  n_buckets <- ncol(mat)

  if (n_sites == 0 || n_buckets == 0) {
    return(list(
      data = .empty_uncertainty_tbl(),
      diagnostics = list(
        n_sites = n_sites,
        n_buckets = n_buckets,
        B = as.integer(B)
      )
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Observed estimate: mean across sites
  estimate <- colMeans(mat)

  boot_means <- matrix(NA_real_, nrow = B, ncol = n_buckets)

  for (b in seq_len(B)) {
    idx <- sample.int(n_sites, size = n_sites, replace = TRUE)
    boot_means[b, ] <- colMeans(mat[idx, , drop = FALSE])
  }

  var_hat <- apply(boot_means, 2, stats::var)
  se <- sqrt(var_hat)

  out <- tibble::tibble(
    horizon_bucket = bucket_tbl$horizon_bucket,
    bucket_start = bucket_tbl$bucket_start,
    bucket_end = bucket_tbl$bucket_end,
    estimate = as.numeric(estimate),
    var_hat = as.numeric(var_hat),
    se = as.numeric(se)
  )

  list(
    data = out,
    diagnostics = list(
      n_sites = n_sites,
      n_buckets = n_buckets,
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
#' @param site_bucket Site-by-bucket table.
#' @param site_col Site identifier column. Currently must be "site_id".
#' @param bucket_col Bucket identifier column. Currently must be "horizon_bucket".
#' @param value_col Value column to analyze, usually "value" or "value_norm".
#'
#' @return A list with:
#'   - data: one row per bucket with estimate, var_hat, se
#'   - diagnostics: jackknife diagnostics
within_variance_site_jackknife <- function(site_bucket,
                                           site_col = "site_id",
                                           bucket_col = "horizon_bucket",
                                           value_col = "value_norm") {
  if (!identical(site_col, "site_id")) {
    rlang::abort("This implementation expects `site_col = \"site_id\"`.")
  }
  if (!identical(bucket_col, "horizon_bucket")) {
    rlang::abort("This implementation expects `bucket_col = \"horizon_bucket\"`.")
  }

  sbm <- .site_bucket_to_matrix(site_bucket, value_col = value_col)
  mat <- sbm$mat
  bucket_tbl <- sbm$bucket_tbl

  n_sites <- nrow(mat)
  n_buckets <- ncol(mat)

  if (n_sites == 0 || n_buckets == 0) {
    return(list(
      data = .empty_uncertainty_tbl(),
      diagnostics = list(
        n_sites = n_sites,
        n_buckets = n_buckets
      )
    ))
  }

  if (n_sites < 2) {
    out <- tibble::tibble(
      horizon_bucket = bucket_tbl$horizon_bucket,
      bucket_start = bucket_tbl$bucket_start,
      bucket_end = bucket_tbl$bucket_end,
      estimate = as.numeric(colMeans(mat)),
      var_hat = NA_real_,
      se = NA_real_
    )

    return(list(
      data = out,
      diagnostics = list(
        n_sites = n_sites,
        n_buckets = n_buckets,
        note = "Jackknife variance undefined for fewer than 2 sites."
      )
    ))
  }

  estimate <- colMeans(mat)

  loo_means <- matrix(NA_real_, nrow = n_sites, ncol = n_buckets)
  for (i in seq_len(n_sites)) {
    loo_means[i, ] <- colMeans(mat[-i, , drop = FALSE])
  }

  loo_bar <- colMeans(loo_means)
  var_hat <- ((n_sites - 1) / n_sites) * colSums((loo_means - matrix(loo_bar, nrow = n_sites, ncol = n_buckets, byrow = TRUE))^2)
  se <- sqrt(var_hat)

  out <- tibble::tibble(
    horizon_bucket = bucket_tbl$horizon_bucket,
    bucket_start = bucket_tbl$bucket_start,
    bucket_end = bucket_tbl$bucket_end,
    estimate = as.numeric(estimate),
    var_hat = as.numeric(var_hat),
    se = as.numeric(se)
  )

  list(
    data = out,
    diagnostics = list(
      n_sites = n_sites,
      n_buckets = n_buckets,
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
#'   replicate_id, horizon_bucket, bucket_start, bucket_end, estimate, var_hat
#' @param conf_level Confidence level.
#' @param use_t If TRUE, use Rubin-style t critical values when possible.
#'
#' @return Tibble with pooled estimate and uncertainty components per bucket.
pool_replicate_uncertainty <- function(replicate_estimates,
                                       conf_level = 0.95,
                                       use_t = TRUE) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort("`conf_level` must be a single number between 0 and 1.")
  }

  replicate_estimates <- tibble::as_tibble(replicate_estimates)
  .require_replicate_estimates(replicate_estimates)

  pooled <- replicate_estimates |>
    dplyr::group_by(horizon_bucket, bucket_start, bucket_end) |>
    dplyr::summarise(
      M = dplyr::n_distinct(replicate_id),
      estimate_mean = mean(estimate, na.rm = TRUE),
      W = mean(var_hat, na.rm = TRUE),
      B = if (dplyr::n() > 1) stats::var(estimate, na.rm = TRUE) else 0,
      .groups = "drop"
    ) |>
    dplyr::rename(estimate = estimate_mean)

  pooled$W[!is.finite(pooled$W)] <- NA_real_
  pooled$B[!is.finite(pooled$B)] <- 0

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
#' @return Long-form tibble with within, between, and total variance/SE by bucket.
summarize_uncertainty_components <- function(pooled_result) {
  pooled_result <- tibble::as_tibble(pooled_result)

  validate_required_fields(
    pooled_result,
    c("horizon_bucket", "bucket_start", "bucket_end", "W", "B", "T", "se_within", "se_between", "se_total"),
    "pooled_result"
  )

  out <- dplyr::bind_rows(
    tibble::tibble(
      horizon_bucket = pooled_result$horizon_bucket,
      bucket_start = pooled_result$bucket_start,
      bucket_end = pooled_result$bucket_end,
      component = "within",
      variance = pooled_result$W,
      se = pooled_result$se_within
    ),
    tibble::tibble(
      horizon_bucket = pooled_result$horizon_bucket,
      bucket_start = pooled_result$bucket_start,
      bucket_end = pooled_result$bucket_end,
      component = "between",
      variance = pooled_result$B,
      se = pooled_result$se_between
    ),
    tibble::tibble(
      horizon_bucket = pooled_result$horizon_bucket,
      bucket_start = pooled_result$bucket_start,
      bucket_end = pooled_result$bucket_end,
      component = "total",
      variance = pooled_result$T,
      se = pooled_result$se_total
    )
  )

  tibble::as_tibble(out)
}
