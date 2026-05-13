# significance.R
#
# Peak significance testing for the merged occupancy package.
#
# Supports significance tests at:
# - site level
# - region level
# - supra-regional / global level
#
# Default peak definition:
#   D_t = mu_t - (mu_{t-1} + mu_{t+1}) / 2
#
# Null hypothesis:
#   H0: D_t <= 0
#
# Optional baseline mode:
#   D_t = mu_t - mean(mu_baseline_window)
#
# Inference:
# - site bootstrap
# - optionally aggregated across multiple temporal replicates
#
# Expected site_bucket schema:
#   site_id, horizon_bucket, bucket_start, bucket_end,
#   value / value_norm, optionally site_region and other metadata

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang
#
# Suggested namespace usage:
#   @importFrom dplyr group_by summarise bind_rows n_distinct
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.require_site_bucket_significance <- function(site_bucket, value_col) {
  validate_required_fields(
    site_bucket,
    c("site_id", "horizon_bucket", "bucket_start", "bucket_end", value_col),
    "site_bucket"
  )
  invisible(TRUE)
}

.extract_site_bucket_tbl <- function(x) {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }

  if (is.list(x) && !inherits(x, "data.frame")) {
    if ("data" %in% names(x) && "site_bucket" %in% names(x$data)) {
      return(tibble::as_tibble(x$data$site_bucket))
    }
    if ("site_bucket" %in% names(x)) {
      return(tibble::as_tibble(x$site_bucket))
    }
  }

  NULL
}

.extract_site_bucket_replicates <- function(x) {
  # plain site_bucket tibble
  if (is.data.frame(x)) {
    return(list(
      replicates = list(tibble::as_tibble(x)),
      bucket_grid = NULL
    ))
  }

  if (is.list(x) && !inherits(x, "data.frame")) {
    # finalized ensemble object
    if ("replicate_data" %in% names(x) && "site_bucket" %in% names(x$replicate_data)) {
      bg <- NULL
      if ("pooled" %in% names(x) && "canonical_grid" %in% names(x$pooled)) {
        bg <- tibble::as_tibble(x$pooled$canonical_grid)
      } else if ("pooled" %in% names(x) && "canonical_grid" %in% names(x$pooled$estimates)) {
        bg <- tibble::as_tibble(x$pooled$estimates$canonical_grid)
      }
      return(list(
        replicates = lapply(x$replicate_data$site_bucket, tibble::as_tibble),
        bucket_grid = bg
      ))
    }

    # analysis object
    if ("data" %in% names(x) && "site_bucket" %in% names(x$data)) {
      return(list(
        replicates = list(tibble::as_tibble(x$data$site_bucket)),
        bucket_grid = NULL
      ))
    }

    # list of objects/tables
    maybe_list <- lapply(x, .extract_site_bucket_tbl)
    maybe_list <- Filter(Negate(is.null), maybe_list)
    if (length(maybe_list) > 0) {
      return(list(
        replicates = maybe_list,
        bucket_grid = NULL
      ))
    }
  }

  rlang::abort(
    "Could not extract site-bucket replicate tables. Supply a site_bucket table, an analysis object, a finalized ensemble object, or a list of such objects."
  )
}

.complete_curve_grid <- function(curve_tbl, bucket_grid) {
  curve_tbl <- tibble::as_tibble(curve_tbl)
  bucket_grid <- tibble::as_tibble(bucket_grid)

  validate_required_fields(
    bucket_grid,
    c("horizon_bucket", "bucket_start", "bucket_end"),
    "bucket_grid"
  )

  validate_required_fields(
    curve_tbl,
    c("horizon_bucket", "bucket_start", "bucket_end", "curve_value"),
    "curve_tbl"
  )

  out <- merge(
    bucket_grid,
    curve_tbl,
    by = c("horizon_bucket", "bucket_start", "bucket_end"),
    all.x = TRUE,
    sort = TRUE
  )

  if (!"n_sites" %in% names(out)) {
    out$n_sites <- NA_integer_
  }

  out$curve_value[is.na(out$curve_value)] <- 0
  out$n_sites[is.na(out$n_sites)] <- 0L

  tibble::as_tibble(out)
}

.get_bucket_order_tbl <- function(site_bucket) {
  unique(site_bucket[, c("horizon_bucket", "bucket_start", "bucket_end"), drop = FALSE]) |>
    tibble::as_tibble() |>
    (\(z) z[order(z$horizon_bucket), , drop = FALSE])()
}

.compute_reference_buckets <- function(bucket_tbl,
                                       target_buckets,
                                       mode = c("neighbors", "baseline"),
                                       baseline_buckets = NULL) {
  mode <- match.arg(mode)

  h <- bucket_tbl$horizon_bucket
  target_buckets <- sort(as.numeric(target_buckets))

  if (!all(target_buckets %in% h)) {
    rlang::abort("Some `target_buckets` were not found in the bucket table.")
  }

  if (mode == "neighbors") {
    idx <- match(target_buckets, h)

    left_idx <- min(idx) - 1L
    right_idx <- max(idx) + 1L

    if (left_idx < 1L || right_idx > length(h)) {
      rlang::abort(
        "Neighbor-based peak testing for a bucket range requires both a left and right neighboring bucket."
      )
    }

    return(c(h[left_idx], h[right_idx]))
  }

  if (is.null(baseline_buckets) || length(baseline_buckets) == 0L) {
    rlang::abort("For `mode = \"baseline\"`, `baseline_buckets` must be supplied.")
  }

  baseline_buckets <- as.numeric(baseline_buckets)
  if (!all(baseline_buckets %in% h)) {
    rlang::abort("Some `baseline_buckets` were not found in the bucket table.")
  }

  baseline_buckets
}

.compute_group_curve <- function(site_bucket,
                                 value_col,
                                 level = c("global", "region", "site"),
                                 level_id = NULL,
                                 region_col = "site_region") {
  level <- match.arg(level)
  site_bucket <- tibble::as_tibble(site_bucket)

  if (level == "global") {
    out <- site_bucket |>
      dplyr::group_by(horizon_bucket, bucket_start, bucket_end) |>
      dplyr::summarise(
        curve_value = mean(.data[[value_col]], na.rm = TRUE),
        n_sites = dplyr::n_distinct(site_id),
        .groups = "drop"
      )
    return(tibble::as_tibble(out))
  }

  if (is.null(level_id) || length(level_id) != 1L || is.na(level_id)) {
    rlang::abort("`level_id` must be supplied for `level = \"region\"` and `level = \"site\"`.")
  }

  if (level == "region") {
    validate_required_fields(site_bucket, region_col, "site_bucket")
    sb <- site_bucket[!is.na(site_bucket[[region_col]]) & site_bucket[[region_col]] == level_id, , drop = FALSE]

    out <- sb |>
      dplyr::group_by(horizon_bucket, bucket_start, bucket_end) |>
      dplyr::summarise(
        curve_value = mean(.data[[value_col]], na.rm = TRUE),
        n_sites = dplyr::n_distinct(site_id),
        .groups = "drop"
      )
    return(tibble::as_tibble(out))
  }

  sb <- site_bucket[!is.na(site_bucket$site_id) & site_bucket$site_id == level_id, , drop = FALSE]

  out <- sb |>
    dplyr::group_by(horizon_bucket, bucket_start, bucket_end) |>
    dplyr::summarise(
      curve_value = mean(.data[[value_col]], na.rm = TRUE),
      n_sites = dplyr::n_distinct(site_id),
      .groups = "drop"
    )

  tibble::as_tibble(out)
}

.compute_peak_contrast <- function(curve_tbl,
                                   target_buckets,
                                   mode = c("neighbors", "baseline"),
                                   baseline_buckets = NULL) {
  mode <- match.arg(mode)
  curve_tbl <- tibble::as_tibble(curve_tbl)

  validate_required_fields(
    curve_tbl,
    c("horizon_bucket", "bucket_start", "bucket_end", "curve_value"),
    "curve_tbl"
  )

  bucket_tbl <- .get_bucket_order_tbl(curve_tbl)

  target_buckets <- sort(as.numeric(target_buckets))
  ref_buckets <- .compute_reference_buckets(
    bucket_tbl = bucket_tbl,
    target_buckets = target_buckets,
    mode = mode,
    baseline_buckets = baseline_buckets
  )

  target_vals <- curve_tbl$curve_value[match(target_buckets, curve_tbl$horizon_bucket)]
  ref_vals <- curve_tbl$curve_value[match(ref_buckets, curve_tbl$horizon_bucket)]

  target_mean <- mean(target_vals, na.rm = TRUE)
  ref_mean <- mean(ref_vals, na.rm = TRUE)
  contrast <- target_mean - ref_mean

  list(
    target_buckets = target_buckets,
    reference_buckets = ref_buckets,
    target_value = target_mean,
    reference_value = ref_mean,
    contrast = contrast
  )
}

.bootstrap_peak_contrast_one <- function(site_bucket,
                                         value_col,
                                         target_buckets,
                                         mode = c("neighbors", "baseline"),
                                         baseline_buckets = NULL,
                                         level = c("global", "region", "site"),
                                         level_id = NULL,
                                         region_col = "site_region",
                                         B = 1000L,
                                         seed = NULL,
                                         bucket_grid = NULL) {
  mode <- match.arg(mode)
  level <- match.arg(level)

  site_bucket <- tibble::as_tibble(site_bucket)
  .require_site_bucket_significance(site_bucket, value_col)

  if (level == "region") {
    validate_required_fields(site_bucket, region_col, "site_bucket")
    site_bucket <- site_bucket[!is.na(site_bucket[[region_col]]) & site_bucket[[region_col]] == level_id, , drop = FALSE]
  }

  if (level == "site") {
    site_bucket <- site_bucket[!is.na(site_bucket$site_id) & site_bucket$site_id == level_id, , drop = FALSE]
  }

  site_ids <- sort(unique(site_bucket$site_id))
  n_sites <- length(site_ids)

  if (n_sites == 0L) {
    rlang::abort("No sites available for the requested test level/group.")
  }

  if (!is.null(seed)) set.seed(seed)

  observed_curve_raw <- .compute_group_curve(
    site_bucket = site_bucket,
    value_col = value_col,
    level = "global"
  )

  # If a canonical bucket grid is supplied, use it.
  # Otherwise derive from the observed subset.
  if (is.null(bucket_grid)) {
    bucket_grid <- .get_bucket_order_tbl(observed_curve_raw)
  } else {
    bucket_grid <- tibble::as_tibble(bucket_grid)
  }

  observed_curve <- .complete_curve_grid(
    curve_tbl = observed_curve_raw,
    bucket_grid = bucket_grid
  )

  observed <- .compute_peak_contrast(
    curve_tbl = observed_curve,
    target_buckets = target_buckets,
    mode = mode,
    baseline_buckets = baseline_buckets
  )

  boot_contrast <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n_sites, size = n_sites, replace = TRUE)
    sampled_sites <- site_ids[idx]

    sampled_tbl <- dplyr::bind_rows(lapply(sampled_sites, function(sid) {
      site_bucket[site_bucket$site_id == sid, , drop = FALSE]
    }))

    boot_curve_raw <- .compute_group_curve(
      site_bucket = sampled_tbl,
      value_col = value_col,
      level = "global"
    )

    boot_curve <- .complete_curve_grid(
      curve_tbl = boot_curve_raw,
      bucket_grid = bucket_grid
    )

    boot_contrast[b] <- .compute_peak_contrast(
      curve_tbl = boot_curve,
      target_buckets = target_buckets,
      mode = mode,
      baseline_buckets = baseline_buckets
    )$contrast
  }

  p_value <- mean(boot_contrast <= 0, na.rm = TRUE)
  ci <- stats::quantile(boot_contrast, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  list(
    observed = observed,
    bootstrap_distribution = boot_contrast,
    p_value = p_value,
    ci_lower = ci[1],
    ci_upper = ci[2],
    diagnostics = list(
      n_sites = n_sites,
      B = as.integer(B),
      level = level,
      level_id = level_id,
      mode = mode
    )
  )
}

# -------------------------------------------------------------------------
# Core multi-replicate test
# -------------------------------------------------------------------------

#' Test peak significance across one or more replicate site-bucket tables
#'
#' @param site_bucket_replicates Either:
#'   - one site_bucket tibble / analysis object, or
#'   - a list of site_bucket tibbles / analysis objects
#' @param target_bucket Bucket to test as a peak.
#' @param level One of "global", "region", "site".
#' @param level_id Required for region/site tests.
#' @param value_col Value column, usually "value_norm" or "value".
#' @param mode "neighbors" (default) or "baseline".
#' @param baseline_buckets Optional vector of baseline buckets.
#' @param region_col Region column name.
#' @param B Number of site-bootstrap resamples per replicate.
#' @param seed Optional seed.
#'
#' @return Structured list with observed contrast, CI, p-value, and replicate details.
test_peak_significance_core <- function(site_bucket_replicates,
                                        target_buckets,
                                        level = c("global", "region", "site"),
                                        level_id = NULL,
                                        value_col = "value_norm",
                                        mode = c("neighbors", "baseline"),
                                        baseline_buckets = NULL,
                                        region_col = "site_region",
                                        B = 1000L,
                                        seed = NULL) {
  level <- match.arg(level)
  mode <- match.arg(mode)

  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 2L) {
    rlang::abort("`B` must be a single integer >= 2.")
  }

  target_buckets <- sort(as.numeric(target_buckets))
  if (length(target_buckets) < 1L) {
    rlang::abort("`target_buckets` must contain at least one bucket.")
  }

  extracted <- .extract_site_bucket_replicates(site_bucket_replicates)
  site_bucket_replicates <- extracted$replicates
  bucket_grid <- extracted$bucket_grid

  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, length(site_bucket_replicates), replace = TRUE)

  rep_results <- vector("list", length(site_bucket_replicates))

  for (i in seq_along(site_bucket_replicates)) {
    sb <- tibble::as_tibble(site_bucket_replicates[[i]])

    rep_results[[i]] <- .bootstrap_peak_contrast_one(
      site_bucket = sb,
      value_col = value_col,
      target_buckets = target_buckets,
      mode = mode,
      baseline_buckets = baseline_buckets,
      level = level,
      level_id = level_id,
      region_col = region_col,
      B = B,
      seed = rep_seeds[i],
      bucket_grid = bucket_grid
    )
  }
  observed_contrasts <- vapply(rep_results, function(x) x$observed$contrast, numeric(1))
  observed_target_vals <- vapply(rep_results, function(x) x$observed$target_value, numeric(1))
  observed_ref_vals <- vapply(rep_results, function(x) x$observed$reference_value, numeric(1))

  pooled_boot <- unlist(lapply(rep_results, `[[`, "bootstrap_distribution"), use.names = FALSE)

  p_value <- mean(pooled_boot <= 0, na.rm = TRUE)
  ci <- stats::quantile(pooled_boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  reference_buckets <- rep_results[[1]]$observed$reference_buckets

  list(
    result = tibble::tibble(
      level = level,
      level_id = if (is.null(level_id)) NA_character_ else as.character(level_id),
      target_buckets = paste(target_buckets, collapse = ","),
      reference_type = mode,
      reference_buckets = paste(reference_buckets, collapse = ","),
      observed_target_value = mean(observed_target_vals, na.rm = TRUE),
      observed_reference_value = mean(observed_ref_vals, na.rm = TRUE),
      observed_contrast = mean(observed_contrasts, na.rm = TRUE),
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_value,
      significant = !is.na(p_value) && p_value < 0.05
    ),
    replicate_results = tibble::tibble(
      replicate_id = seq_along(rep_results),
      observed_contrast = observed_contrasts,
      observed_target_value = observed_target_vals,
      observed_reference_value = observed_ref_vals
    ),
    bootstrap_distribution = pooled_boot,
    diagnostics = list(
      n_replicates = length(rep_results),
      B_per_replicate = as.integer(B),
      level = level,
      level_id = level_id,
      mode = mode,
      value_col = value_col
    )
  )
}

# -------------------------------------------------------------------------
# Convenience wrappers
# -------------------------------------------------------------------------

#' Test peak significance for the supra-regional/global curve
test_peak_significance_global <- function(site_bucket_replicates,
                                          target_buckets,
                                          value_col = "value_norm",
                                          mode = c("neighbors", "baseline"),
                                          baseline_buckets = NULL,
                                          B = 1000L,
                                          seed = NULL) {
  test_peak_significance_core(
    site_bucket_replicates = site_bucket_replicates,
    target_buckets = target_buckets,
    level = "global",
    level_id = NULL,
    value_col = value_col,
    mode = mode,
    baseline_buckets = baseline_buckets,
    B = B,
    seed = seed
  )
}

#' Test peak significance for one regional curve
test_peak_significance_region <- function(site_bucket_replicates,
                                          region_id,
                                          target_buckets,
                                          value_col = "value_norm",
                                          mode = c("neighbors", "baseline"),
                                          baseline_buckets = NULL,
                                          region_col = "site_region",
                                          B = 1000L,
                                          seed = NULL) {
  test_peak_significance_core(
    site_bucket_replicates = site_bucket_replicates,
    target_buckets = target_buckets,
    level = "region",
    level_id = region_id,
    value_col = value_col,
    mode = mode,
    baseline_buckets = baseline_buckets,
    region_col = region_col,
    B = B,
    seed = seed
  )
}

#' Test peak significance for one site curve
test_peak_significance_site <- function(site_bucket_replicates,
                                        site_id,
                                        target_buckets,
                                        value_col = "value_norm",
                                        mode = c("neighbors", "baseline"),
                                        baseline_buckets = NULL,
                                        B = 1000L,
                                        seed = NULL) {
  test_peak_significance_core(
    site_bucket_replicates = site_bucket_replicates,
    target_buckets = target_buckets,
    level = "site",
    level_id = site_id,
    value_col = value_col,
    mode = mode,
    baseline_buckets = baseline_buckets,
    B = B,
    seed = seed
  )
}
