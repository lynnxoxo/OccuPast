# significance.R
#
# Peak significance testing for OccuPast.
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
# Expected site_bin schema:
#   site_id, horizon_bin, bin_start, bin_end,
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

.require_site_bin_significance <- function(site_bin, value_col) {
  validate_required_fields(
    site_bin,
    c("site_id", "horizon_bin", "bin_start", "bin_end", value_col),
    "site_bin"
  )
  invisible(TRUE)
}

.extract_site_bin_tbl <- function(x) {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }

  if (is.list(x) && !inherits(x, "data.frame")) {
    if ("data" %in% names(x) && "site_bin" %in% names(x$data)) {
      return(tibble::as_tibble(x$data$site_bin))
    }
    if ("site_bin" %in% names(x)) {
      return(tibble::as_tibble(x$site_bin))
    }
  }

  NULL
}

.extract_site_bin_replicates <- function(x) {
  # plain site_bin tibble
  if (is.data.frame(x)) {
    return(list(
      replicates = list(tibble::as_tibble(x)),
      bin_grid = NULL
    ))
  }

  if (is.list(x) && !inherits(x, "data.frame")) {
    # finalized ensemble object
    if ("replicate_data" %in% names(x) && "site_bin" %in% names(x$replicate_data)) {
      bg <- NULL
      if ("pooled" %in% names(x) && "canonical_grid" %in% names(x$pooled)) {
        bg <- tibble::as_tibble(x$pooled$canonical_grid)
      } else if ("pooled" %in% names(x) && "canonical_grid" %in% names(x$pooled$estimates)) {
        bg <- tibble::as_tibble(x$pooled$estimates$canonical_grid)
      }
      return(list(
        replicates = lapply(x$replicate_data$site_bin, tibble::as_tibble),
        bin_grid = bg
      ))
    }

    # analysis object
    if ("data" %in% names(x) && "site_bin" %in% names(x$data)) {
      return(list(
        replicates = list(tibble::as_tibble(x$data$site_bin)),
        bin_grid = NULL
      ))
    }

    # list of objects/tables
    maybe_list <- lapply(x, .extract_site_bin_tbl)
    maybe_list <- Filter(Negate(is.null), maybe_list)
    if (length(maybe_list) > 0) {
      return(list(
        replicates = maybe_list,
        bin_grid = NULL
      ))
    }
  }

  rlang::abort(
    "Could not extract site-bin replicate tables. Supply a site_bin table, an analysis object, a finalized ensemble object, or a list of such objects."
  )
}

.complete_curve_grid <- function(curve_tbl, bin_grid) {
  curve_tbl <- tibble::as_tibble(curve_tbl)
  bin_grid <- tibble::as_tibble(bin_grid)

  validate_required_fields(
    bin_grid,
    c("horizon_bin", "bin_start", "bin_end"),
    "bin_grid"
  )

  validate_required_fields(
    curve_tbl,
    c("horizon_bin", "bin_start", "bin_end", "curve_value"),
    "curve_tbl"
  )

  out <- merge(
    bin_grid,
    curve_tbl,
    by = c("horizon_bin", "bin_start", "bin_end"),
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

.get_bin_order_tbl <- function(site_bin) {
  unique(site_bin[, c("horizon_bin", "bin_start", "bin_end"), drop = FALSE]) |>
    tibble::as_tibble() |>
    (\(z) z[order(z$horizon_bin), , drop = FALSE])()
}

.compute_reference_bins <- function(bin_tbl,
                                    target_bins,
                                    mode = c("neighbors", "baseline"),
                                    baseline_bins = NULL) {
  mode <- match.arg(mode)

  h <- bin_tbl$horizon_bin
  target_bins <- sort(as.numeric(target_bins))

  if (!all(target_bins %in% h)) {
    rlang::abort("Some `target_bins` were not found in the bin table.")
  }

  if (mode == "neighbors") {
    idx <- match(target_bins, h)

    left_idx <- min(idx) - 1L
    right_idx <- max(idx) + 1L

    if (left_idx < 1L || right_idx > length(h)) {
      rlang::abort(
        "Neighbor-based peak testing for a bin range requires both a left and right neighboring bin."
      )
    }

    return(c(h[left_idx], h[right_idx]))
  }

  if (is.null(baseline_bins) || length(baseline_bins) == 0L) {
    rlang::abort("For `mode = \"baseline\"`, `baseline_bins` must be supplied.")
  }

  baseline_bins <- as.numeric(baseline_bins)
  if (!all(baseline_bins %in% h)) {
    rlang::abort("Some `baseline_bins` were not found in the bin table.")
  }

  baseline_bins
}

.compute_group_curve <- function(site_bin,
                                 value_col,
                                 level = c("global", "region", "site"),
                                 level_id = NULL,
                                 region_col = "site_region") {
  level <- match.arg(level)
  site_bin <- tibble::as_tibble(site_bin)

  if (level == "global") {
    out <- site_bin |>
      dplyr::group_by(horizon_bin, bin_start, bin_end) |>
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
    validate_required_fields(site_bin, region_col, "site_bin")
    sb <- site_bin[!is.na(site_bin[[region_col]]) & site_bin[[region_col]] == level_id, , drop = FALSE]

    out <- sb |>
      dplyr::group_by(horizon_bin, bin_start, bin_end) |>
      dplyr::summarise(
        curve_value = mean(.data[[value_col]], na.rm = TRUE),
        n_sites = dplyr::n_distinct(site_id),
        .groups = "drop"
      )
    return(tibble::as_tibble(out))
  }

  sb <- site_bin[!is.na(site_bin$site_id) & site_bin$site_id == level_id, , drop = FALSE]

  out <- sb |>
    dplyr::group_by(horizon_bin, bin_start, bin_end) |>
    dplyr::summarise(
      curve_value = mean(.data[[value_col]], na.rm = TRUE),
      n_sites = dplyr::n_distinct(site_id),
      .groups = "drop"
    )

  tibble::as_tibble(out)
}

# Peak contrast: D = mean(target bins) - mean(reference bins).
# Neighbor mode uses the bins immediately left and right of the target span.
.compute_peak_contrast <- function(curve_tbl,
                                   target_bins,
                                   mode = c("neighbors", "baseline"),
                                   baseline_bins = NULL) {
  mode <- match.arg(mode)
  curve_tbl <- tibble::as_tibble(curve_tbl)

  validate_required_fields(
    curve_tbl,
    c("horizon_bin", "bin_start", "bin_end", "curve_value"),
    "curve_tbl"
  )

  bin_tbl <- .get_bin_order_tbl(curve_tbl)

  target_bins <- sort(as.numeric(target_bins))
  ref_bins <- .compute_reference_bins(
    bin_tbl = bin_tbl,
    target_bins = target_bins,
    mode = mode,
    baseline_bins = baseline_bins
  )

  target_vals <- curve_tbl$curve_value[match(target_bins, curve_tbl$horizon_bin)]
  ref_vals <- curve_tbl$curve_value[match(ref_bins, curve_tbl$horizon_bin)]

  target_mean <- mean(target_vals, na.rm = TRUE)
  ref_mean <- mean(ref_vals, na.rm = TRUE)
  contrast <- target_mean - ref_mean

  list(
    target_bins = target_bins,
    reference_bins = ref_bins,
    target_value = target_mean,
    reference_value = ref_mean,
    contrast = contrast
  )
}

.bootstrap_peak_contrast_one <- function(site_bin,
                                         value_col,
                                         target_bins,
                                         mode = c("neighbors", "baseline"),
                                         baseline_bins = NULL,
                                         level = c("global", "region", "site"),
                                         level_id = NULL,
                                         region_col = "site_region",
                                         B = 1000L,
                                         seed = NULL,
                                         bin_grid = NULL) {
  mode <- match.arg(mode)
  level <- match.arg(level)

  site_bin <- tibble::as_tibble(site_bin)
  .require_site_bin_significance(site_bin, value_col)

  if (level == "region") {
    validate_required_fields(site_bin, region_col, "site_bin")
    site_bin <- site_bin[!is.na(site_bin[[region_col]]) & site_bin[[region_col]] == level_id, , drop = FALSE]
  }

  if (level == "site") {
    site_bin <- site_bin[!is.na(site_bin$site_id) & site_bin$site_id == level_id, , drop = FALSE]
  }

  site_ids <- sort(unique(site_bin$site_id))
  n_sites <- length(site_ids)

  if (n_sites == 0L) {
    rlang::abort("No sites available for the requested test level/group.")
  }

  if (!is.null(seed)) set.seed(seed)

  observed_curve_raw <- .compute_group_curve(
    site_bin = site_bin,
    value_col = value_col,
    level = "global"
  )

  # If a canonical bin grid is supplied, use it.
  # Otherwise derive from the observed subset.
  if (is.null(bin_grid)) {
    bin_grid <- .get_bin_order_tbl(observed_curve_raw)
  } else {
    bin_grid <- tibble::as_tibble(bin_grid)
  }

  observed_curve <- .complete_curve_grid(
    curve_tbl = observed_curve_raw,
    bin_grid = bin_grid
  )

  observed <- .compute_peak_contrast(
    curve_tbl = observed_curve,
    target_bins = target_bins,
    mode = mode,
    baseline_bins = baseline_bins
  )

  boot_contrast <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    idx <- sample.int(n_sites, size = n_sites, replace = TRUE)
    sampled_sites <- site_ids[idx]

    sampled_tbl <- dplyr::bind_rows(lapply(sampled_sites, function(sid) {
      site_bin[site_bin$site_id == sid, , drop = FALSE]
    }))

    boot_curve_raw <- .compute_group_curve(
      site_bin = sampled_tbl,
      value_col = value_col,
      level = "global"
    )

    boot_curve <- .complete_curve_grid(
      curve_tbl = boot_curve_raw,
      bin_grid = bin_grid
    )

    boot_contrast[b] <- .compute_peak_contrast(
      curve_tbl = boot_curve,
      target_bins = target_bins,
      mode = mode,
      baseline_bins = baseline_bins
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

#' Test peak significance across one or more replicate site-bin tables
#'
#' @param site_bin_replicates Either:
#'   - one site_bin tibble / analysis object, or
#'   - a list of site_bin tibbles / analysis objects
#' @param target_bin Bin to test as a peak.
#' @param level One of "global", "region", "site".
#' @param level_id Required for region/site tests.
#' @param value_col Value column, usually "value_norm" or "value".
#' @param mode "neighbors" (default) or "baseline".
#' @param baseline_bins Optional vector of baseline bins.
#' @param region_col Region column name.
#' @param B Number of site-bootstrap resamples per replicate.
#' @param seed Optional seed.
#'
#' @return Structured list with observed contrast, CI, p-value, and replicate details.
test_peak_significance_core <- function(site_bin_replicates,
                                        target_bins,
                                        level = c("global", "region", "site"),
                                        level_id = NULL,
                                        value_col = "value_norm",
                                        mode = c("neighbors", "baseline"),
                                        baseline_bins = NULL,
                                        region_col = "site_region",
                                        B = 1000L,
                                        seed = NULL) {
  level <- match.arg(level)
  mode <- match.arg(mode)

  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 2L) {
    rlang::abort("`B` must be a single integer >= 2.")
  }

  target_bins <- sort(as.numeric(target_bins))
  if (length(target_bins) < 1L) {
    rlang::abort("`target_bins` must contain at least one bin.")
  }

  extracted <- .extract_site_bin_replicates(site_bin_replicates)
  site_bin_replicates <- extracted$replicates
  bin_grid <- extracted$bin_grid

  if (!is.null(seed)) set.seed(seed)
  rep_seeds <- sample.int(.Machine$integer.max, length(site_bin_replicates), replace = TRUE)

  rep_results <- vector("list", length(site_bin_replicates))

  for (i in seq_along(site_bin_replicates)) {
    sb <- tibble::as_tibble(site_bin_replicates[[i]])

    rep_results[[i]] <- .bootstrap_peak_contrast_one(
      site_bin = sb,
      value_col = value_col,
      target_bins = target_bins,
      mode = mode,
      baseline_bins = baseline_bins,
      level = level,
      level_id = level_id,
      region_col = region_col,
      B = B,
      seed = rep_seeds[i],
      bin_grid = bin_grid
    )
  }
  observed_contrasts <- vapply(rep_results, function(x) x$observed$contrast, numeric(1))
  observed_target_vals <- vapply(rep_results, function(x) x$observed$target_value, numeric(1))
  observed_ref_vals <- vapply(rep_results, function(x) x$observed$reference_value, numeric(1))

  pooled_boot <- unlist(lapply(rep_results, `[[`, "bootstrap_distribution"), use.names = FALSE)

  p_value <- mean(pooled_boot <= 0, na.rm = TRUE)
  ci <- stats::quantile(pooled_boot, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  reference_bins <- rep_results[[1]]$observed$reference_bins

  list(
    result = tibble::tibble(
      level = level,
      level_id = if (is.null(level_id)) NA_character_ else as.character(level_id),
      target_bins = paste(target_bins, collapse = ","),
      reference_type = mode,
      reference_bins = paste(reference_bins, collapse = ","),
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
test_peak_significance_global <- function(site_bin_replicates,
                                          target_bins,
                                          value_col = "value_norm",
                                          mode = c("neighbors", "baseline"),
                                          baseline_bins = NULL,
                                          B = 1000L,
                                          seed = NULL) {
  test_peak_significance_core(
    site_bin_replicates = site_bin_replicates,
    target_bins = target_bins,
    level = "global",
    level_id = NULL,
    value_col = value_col,
    mode = mode,
    baseline_bins = baseline_bins,
    B = B,
    seed = seed
  )
}

#' Test peak significance for one regional curve
test_peak_significance_region <- function(site_bin_replicates,
                                          region_id,
                                          target_bins,
                                          value_col = "value_norm",
                                          mode = c("neighbors", "baseline"),
                                          baseline_bins = NULL,
                                          region_col = "site_region",
                                          B = 1000L,
                                          seed = NULL) {
  test_peak_significance_core(
    site_bin_replicates = site_bin_replicates,
    target_bins = target_bins,
    level = "region",
    level_id = region_id,
    value_col = value_col,
    mode = mode,
    baseline_bins = baseline_bins,
    region_col = region_col,
    B = B,
    seed = seed
  )
}

#' Test peak significance for one site curve
test_peak_significance_site <- function(site_bin_replicates,
                                        site_id,
                                        target_bins,
                                        value_col = "value_norm",
                                        mode = c("neighbors", "baseline"),
                                        baseline_bins = NULL,
                                        B = 1000L,
                                        seed = NULL) {
  test_peak_significance_core(
    site_bin_replicates = site_bin_replicates,
    target_bins = target_bins,
    level = "site",
    level_id = site_id,
    value_col = value_col,
    mode = mode,
    baseline_bins = baseline_bins,
    B = B,
    seed = seed
  )
}


# -------------------------------------------------------------------------
# Legacy aliases
# -------------------------------------------------------------------------
.extract_site_bucket_tbl <- .extract_site_bin_tbl
.extract_site_bucket_replicates <- .extract_site_bin_replicates
.require_site_bucket_significance <- .require_site_bin_significance
