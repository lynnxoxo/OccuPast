# ensemble.R
#
# Ensemble orchestration for OccuPast.
#
# This file coordinates:
# - temporal replicate generation
# - allocation per replicate
# - analytical summarization per replicate
# - within-replicate site bootstrap / jackknife
# - pooling across replicates
# - optional parallel execution at the replicate level
#
# Assumes these functions already exist:
#   from allocate_time.R:
#     harmonize_chronologies_merged()
#   from analyze_curve.R:
#     analyze_occupancy_curve()
#   from uncertainty.R:
#     within_variance_site_bootstrap()
#     within_variance_site_jackknife()
#     pool_replicate_uncertainty()
#     summarize_uncertainty_components()

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   parallel,
#   magrittr
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows left_join
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort
#   @importFrom magrittr %>%

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.assert_replicate_count <- function(M) {
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || M < 1) {
    rlang::abort("`M` must be a single integer >= 1.")
  }
}

.assert_offset_mode <- function(offset_mode) {
  allowed <- c("cycle", "random", "fixed")
  if (!offset_mode %in% allowed) {
    rlang::abort(
      paste0("`offset_mode` must be one of: ", paste(allowed, collapse = ", "))
    )
  }
}

.assert_parallel_args <- function(parallel, n_cores) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    rlang::abort("`parallel` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) || n_cores < 1) {
    rlang::abort("`n_cores` must be a single integer >= 1.")
  }
}

.derive_replicate_seeds <- function(M, seed = NULL) {
  if (is.null(seed)) {
    return(rep(NA_integer_, M))
  }
  set.seed(seed)
  sample.int(.Machine$integer.max, M, replace = TRUE)
}

# -------------------------------------------------------------------------
# Canonical regridding helpers
# -------------------------------------------------------------------------

#' Build a canonical bin grid from replicate estimates
#'
#' @param replicate_estimates Combined replicate estimates table.
#' @param bin_width Canonical bin width.
#' @param canonical_offset Offset for the canonical grid.
#'
#' @return Tibble with horizon_bin, bin_start, bin_end.
make_canonical_grid_from_replicates <- function(replicate_estimates,
                                                bin_width,
                                                canonical_offset = 0) {
  replicate_estimates <- tibble::as_tibble(replicate_estimates)
  if (exists(".standardize_bin_columns", mode = "function")) {
    replicate_estimates <- .standardize_bin_columns(replicate_estimates, arg = "replicate_estimates")
  }

  validate_required_fields(
    replicate_estimates,
    c("bin_start", "bin_end"),
    "replicate_estimates"
  )

  year_min <- min(replicate_estimates$bin_start, na.rm = TRUE)
  year_max <- max(replicate_estimates$bin_end, na.rm = TRUE)

  make_bin_grid(
    year_min = year_min,
    year_max = year_max,
    bin_width = bin_width,
    offset = canonical_offset
  )
}

#' Regrid one replicate estimate table onto a canonical grid by interval overlap
#'
#' @param replicate_tbl One replicate estimate table with:
#'   replicate_id, horizon_bin, bin_start, bin_end, estimate, var_hat, se
#' @param canonical_grid Canonical grid from make_bin_grid().
#'
#' @return Tibble on the canonical grid for that replicate.
regrid_replicate_estimates <- function(replicate_tbl, canonical_grid) {
  replicate_tbl <- tibble::as_tibble(replicate_tbl)
  canonical_grid <- tibble::as_tibble(canonical_grid)
  if (exists(".standardize_bin_columns", mode = "function")) {
    replicate_tbl <- .standardize_bin_columns(replicate_tbl, arg = "replicate_tbl")
    canonical_grid <- .standardize_bin_columns(canonical_grid, arg = "canonical_grid")
  }

  validate_required_fields(
    replicate_tbl,
    c("replicate_id", "bin_start", "bin_end", "estimate", "var_hat"),
    "replicate_tbl"
  )

  validate_required_fields(
    canonical_grid,
    c("horizon_bin", "bin_start", "bin_end"),
    "canonical_grid"
  )

  if (nrow(replicate_tbl) == 0 || nrow(canonical_grid) == 0) {
    return(tibble::tibble(
      replicate_id = numeric(),
      horizon_bin = numeric(),
      bin_start = numeric(),
      bin_end = numeric(),
      estimate = numeric(),
      var_hat = numeric(),
      se = numeric()
    ))
  }

  rep_id <- unique(replicate_tbl$replicate_id)
  if (length(rep_id) != 1L) {
    rlang::abort("`regrid_replicate_estimates()` expects exactly one replicate_id per input table.")
  }

  out_parts <- vector("list", nrow(canonical_grid))

  for (j in seq_len(nrow(canonical_grid))) {
    c_start <- canonical_grid$bin_start[j]
    c_end   <- canonical_grid$bin_end[j]
    c_width <- c_end - c_start

    overlaps <- pmax(
      0,
      pmin(replicate_tbl$bin_end, c_end) - pmax(replicate_tbl$bin_start, c_start)
    )

    keep <- overlaps > 0
    if (!any(keep)) {
      out_parts[[j]] <- tibble::tibble(
        replicate_id = rep_id,
        horizon_bin = canonical_grid$horizon_bin[j],
        bin_start = c_start,
        bin_end = c_end,
        estimate = 0,
        var_hat = 0,
        se = 0
      )
      next
    }

    # Regridding preserves mass by interval overlap: contribution = estimate * overlap / original_width.
    # Variance is scaled by the squared fraction, following Var(aX) = a^2 Var(X).
    # Fraction of each replicate interval contributing to this canonical interval
    rep_widths <- replicate_tbl$bin_end[keep] - replicate_tbl$bin_start[keep]
    frac <- overlaps[keep] / rep_widths

    # Redistribute estimate by overlap fraction
    est_j <- sum(replicate_tbl$estimate[keep] * frac, na.rm = TRUE)

    # Redistribute variance conservatively by squared fraction
    var_j <- sum(replicate_tbl$var_hat[keep] * (frac^2), na.rm = TRUE)

    out_parts[[j]] <- tibble::tibble(
      replicate_id = rep_id,
      horizon_bin = canonical_grid$horizon_bin[j],
      bin_start = c_start,
      bin_end = c_end,
      estimate = est_j,
      var_hat = var_j,
      se = sqrt(var_j)
    )
  }

  dplyr::bind_rows(out_parts)
}

#' Regrid one regional replicate estimate table onto a canonical grid by interval overlap
#'
#' @param region_replicate_tbl One replicate regional estimate table with:
#'   replicate_id, site_region, horizon_bin, bin_start, bin_end, estimate, var_hat, se
#' @param canonical_grid Canonical grid from make_bin_grid().
#'
#' @return Tibble on the canonical grid for that replicate and all regions.
regrid_region_replicate_estimates <- function(region_replicate_tbl, canonical_grid) {
  region_replicate_tbl <- tibble::as_tibble(region_replicate_tbl)
  canonical_grid <- tibble::as_tibble(canonical_grid)
  if (exists(".standardize_bin_columns", mode = "function")) {
    region_replicate_tbl <- .standardize_bin_columns(region_replicate_tbl, arg = "region_replicate_tbl")
    canonical_grid <- .standardize_bin_columns(canonical_grid, arg = "canonical_grid")
  }

  validate_required_fields(
    region_replicate_tbl,
    c("replicate_id", "site_region", "bin_start", "bin_end", "estimate", "var_hat"),
    "region_replicate_tbl"
  )

  validate_required_fields(
    canonical_grid,
    c("horizon_bin", "bin_start", "bin_end"),
    "canonical_grid"
  )

  if (nrow(region_replicate_tbl) == 0 || nrow(canonical_grid) == 0) {
    return(tibble::tibble(
      replicate_id = numeric(),
      site_region = character(),
      horizon_bin = numeric(),
      bin_start = numeric(),
      bin_end = numeric(),
      estimate = numeric(),
      var_hat = numeric(),
      se = numeric()
    ))
  }

  rep_id <- unique(region_replicate_tbl$replicate_id)
  if (length(rep_id) != 1L) {
    rlang::abort("`regrid_region_replicate_estimates()` expects exactly one replicate_id per input table.")
  }

  regs <- unique(region_replicate_tbl$site_region)
  regs <- regs[!is.na(regs)]

  out_all <- vector("list", length(regs))

  for (r in seq_along(regs)) {
    reg <- regs[r]
    reg_tbl <- region_replicate_tbl[region_replicate_tbl$site_region == reg, , drop = FALSE]

    out_parts <- vector("list", nrow(canonical_grid))

    for (j in seq_len(nrow(canonical_grid))) {
      c_start <- canonical_grid$bin_start[j]
      c_end   <- canonical_grid$bin_end[j]

      overlaps <- pmax(
        0,
        pmin(reg_tbl$bin_end, c_end) - pmax(reg_tbl$bin_start, c_start)
      )

      keep <- overlaps > 0
      if (!any(keep)) {
        out_parts[[j]] <- tibble::tibble(
          replicate_id = rep_id,
          site_region = reg,
          horizon_bin = canonical_grid$horizon_bin[j],
          bin_start = c_start,
          bin_end = c_end,
          estimate = 0,
          var_hat = 0,
          se = 0
        )
        next
      }

      # Same overlap-based regridding as the global curve, applied within region.
      rep_widths <- reg_tbl$bin_end[keep] - reg_tbl$bin_start[keep]
      frac <- overlaps[keep] / rep_widths

      est_j <- sum(reg_tbl$estimate[keep] * frac, na.rm = TRUE)
      var_j <- sum(reg_tbl$var_hat[keep] * (frac^2), na.rm = TRUE)

      out_parts[[j]] <- tibble::tibble(
        replicate_id = rep_id,
        site_region = reg,
        horizon_bin = canonical_grid$horizon_bin[j],
        bin_start = c_start,
        bin_end = c_end,
        estimate = est_j,
        var_hat = var_j,
        se = sqrt(var_j)
      )
    }

    out_all[[r]] <- dplyr::bind_rows(out_parts)
  }

  dplyr::bind_rows(out_all)
}

# -------------------------------------------------------------------------
# Offset generation
# -------------------------------------------------------------------------

#' Generate temporal offsets for replicate runs
#'
#' @param M Number of replicates.
#' @param bin_width Bin width.
#' @param mode One of "cycle", "random", "fixed".
#' @param fixed_offset Used when `mode = "fixed"`.
#' @param seed Optional seed.
#'
#' @return Numeric vector of length M.
generate_temporal_offsets <- function(M,
                                      bin_width,
                                      mode = c("cycle", "random", "fixed"),
                                      fixed_offset = 0,
                                      seed = NULL) {
  mode <- match.arg(mode)
  .assert_replicate_count(M)

  if (!is.numeric(bin_width) || length(bin_width) != 1L ||
      is.na(bin_width) || bin_width <= 0) {
    rlang::abort("`bin_width` must be a single positive number.")
  }

  if (mode == "fixed") {
    return(rep(as.numeric(fixed_offset), M))
  }

  if (mode == "cycle") {
    # Offset cycling tests whether the curve depends on arbitrary bin alignment.
    # Offsets are a nuisance discretisation choice, not an archaeological claim.
    if (M == 1L) return(0)
    return(seq(0, bin_width - 1, length.out = M))
  }

  # random
  if (!is.null(seed)) set.seed(seed)
  stats::runif(M, min = 0, max = bin_width)
}

# -------------------------------------------------------------------------
# Single replicate runner
# -------------------------------------------------------------------------

#' Run one temporal replicate
#'
#' @description
#' Runs one full replicate of the occupancy workflow:
#' chronology allocation, curve analysis, and within-replicate uncertainty.
#'
#' The default temporal allocation mode is \code{"stochastic"}, so each burial is
#' assigned to one time bin by random draw from its phase-level trapezoidal
#' profile. This makes each replicate a plausible temporal realization and gives
#' meaningful between-replicate variance when multiple replicates are pooled.
#'
#' @param prepared_inputs Output from \code{assemble_prepared_inputs()}.
#' @param replicate_id Integer replicate ID.
#' @param bin_width Bin width.
#' @param offset Temporal offset for this replicate.
#' @param allocation_mode One of \code{"stochastic"}, \code{"fractional"}, or
#'   \code{"deterministic_max"}. Default is \code{"stochastic"}.
#' @param profile Default profile for phase-system matched rows.
#' @param normalization One of \code{"none"}, \code{"site_size"},
#'   \code{"site_duration"}, or \code{"site_size_duration"}.
#' @param variance_method One of \code{"site_bootstrap"} or
#'   \code{"site_jackknife"}.
#' @param bootstrap_B Number of bootstrap resamples if bootstrap is used.
#' @param use_site_metadata_fallback Whether to rescue unmatched burials with site metadata bounds.
#' @param augment_to_site_size Whether to generate synthetic burials up to metadata site_size.
#' @param seed Optional replicate seed.
#'
#' @return Structured list for one replicate.
run_one_replicate <- function(prepared_inputs,
                              replicate_id,
                              bin_width = 25,
                              offset = 0,
                              allocation_mode = c("stochastic", "fractional", "deterministic_max"),
                              profile = c("trapezoid", "uniform"),
                              normalization = c("none", "site_size", "site_duration", "site_size_duration"),
                              variance_method = c("site_bootstrap", "site_jackknife"),
                              bootstrap_B = 200L,
                              use_site_metadata_fallback = FALSE,
                              augment_to_site_size = FALSE,
                              seed = NULL,
                              retain_burial_allocations = FALSE) {
  allocation_mode <- match.arg(allocation_mode)
  profile <- match.arg(profile)
  normalization <- match.arg(normalization)
  variance_method <- match.arg(variance_method)

  if (!is.list(prepared_inputs) || is.null(prepared_inputs$data)) {
    rlang::abort("`prepared_inputs` must be an assembled object with a `$data` element.")
  }
  if (is.null(prepared_inputs$data$mortuary) || is.null(prepared_inputs$data$chronology)) {
    rlang::abort("`prepared_inputs$data` must contain `mortuary` and `chronology`.")
  }

  alloc <- harmonize_chronologies_merged(
    prepared_inputs = prepared_inputs,
    bin_width = bin_width,
    offset = offset,
    allocation_mode = allocation_mode,
    profile = profile,
    seed = seed,
    use_site_metadata_fallback = use_site_metadata_fallback,
    augment_to_site_size = augment_to_site_size
  )

  analysis <- analyze_occupancy_curve(
    allocation_result = alloc,
    site_metadata = prepared_inputs$data$site_metadata,
    normalization = normalization
  )

  site_bin <- analysis$data$site_bin
  value_col <- analysis$settings$value_col

  if (variance_method == "site_bootstrap") {
    var_res <- within_variance_site_bootstrap(
      site_bin = site_bin,
      value_col = value_col,
      B = bootstrap_B,
      seed = seed
    )
  } else {
    var_res <- within_variance_site_jackknife(
      site_bin = site_bin,
      value_col = value_col
    )
  }

  rep_curve <- analysis$data$bin_curve %>%
    dplyr::select(horizon_bin, bin_start, bin_end, mean_value) %>%
    dplyr::rename(estimate = mean_value)

  rep_est <- dplyr::left_join(
    rep_curve,
    var_res$data[, c("horizon_bin", "bin_start", "bin_end", "var_hat", "se"), drop = FALSE],
    by = c("horizon_bin", "bin_start", "bin_end")
  )

  rep_est$replicate_id <- replicate_id
  rep_est <- rep_est[, c("replicate_id", "horizon_bin", "bin_start", "bin_end", "estimate", "var_hat", "se")]

  rep_region_est <- NULL
  if (!is.null(analysis$data$region_bin)) {
    reg_tbl <- analysis$data$region_bin

    if ("site_region" %in% names(reg_tbl) && "mean_value" %in% names(reg_tbl)) {
      reg_curve <- reg_tbl %>%
        dplyr::select(site_region, horizon_bin, bin_start, bin_end, mean_value) %>%
        dplyr::rename(estimate = mean_value)

      # Approximate within-replicate regional variance using site-bin subsets per region
      site_bin_tbl <- analysis$data$site_bin
      reg_var_parts <- lapply(unique(reg_curve$site_region), function(reg) {
        sb_reg <- site_bin_tbl[!is.na(site_bin_tbl$site_region) & site_bin_tbl$site_region == reg, , drop = FALSE]

        if (nrow(sb_reg) == 0) return(NULL)

        if (variance_method == "site_bootstrap") {
          vr <- within_variance_site_bootstrap(
            site_bin = sb_reg,
            value_col = value_col,
            B = bootstrap_B,
            seed = seed
          )
        } else {
          vr <- within_variance_site_jackknife(
            site_bin = sb_reg,
            value_col = value_col
          )
        }

        vr$data$site_region <- reg
        vr$data
      })

      reg_var_tbl <- dplyr::bind_rows(reg_var_parts)

      rep_region_est <- dplyr::left_join(
        reg_curve,
        reg_var_tbl[, c("site_region", "horizon_bin", "bin_start", "bin_end", "var_hat", "se"), drop = FALSE],
        by = c("site_region", "horizon_bin", "bin_start", "bin_end")
      )

      rep_region_est$replicate_id <- replicate_id
      rep_region_est <- rep_region_est[, c(
        "replicate_id", "site_region",
        "horizon_bin", "bin_start", "bin_end",
        "estimate", "var_hat", "se"
      )]
    }
  }
  list(
    replicate_id = replicate_id,
    offset = offset,
    allocation = alloc,
    analysis = analysis,
    variance = var_res,
    replicate_estimates = tibble::as_tibble(rep_est),
    region_replicate_estimates = if (is.null(rep_region_est)) NULL else tibble::as_tibble(rep_region_est),
    burial_allocations = if (isTRUE(retain_burial_allocations)) tibble::as_tibble(alloc$data) else NULL,
    diagnostics = list(
      allocation = alloc$diagnostics,
      analysis = analysis$diagnostics,
      variance = var_res$diagnostics
    ),
    settings = list(
      allocation_mode = allocation_mode,
      profile = profile,
      normalization = normalization,
      variance_method = variance_method,
      bin_width = bin_width,
      offset = offset,
      use_site_metadata_fallback = use_site_metadata_fallback,
      augment_to_site_size = augment_to_site_size,
      retain_burial_allocations = retain_burial_allocations
    )
  )
}
# -------------------------------------------------------------------------
# Ensemble runner
# -------------------------------------------------------------------------

#' Run a temporal ensemble
#'
#' @description
#' Runs multiple temporal replicates of the occupancy workflow and stores the
#' per-replicate results needed for later pooling, plotting, and significance
#' testing.
#'
#' By default, temporal allocation is \code{"stochastic"}, so each replicate is a
#' Monte Carlo realization in which each burial is assigned to one bin drawn
#' from its phase-level trapezoidal temporal profile. This makes the ensemble's
#' between-replicate component interpretable as temporal allocation uncertainty.
#'
#' @param prepared_inputs Output from \code{assemble_prepared_inputs()}.
#' @param M Number of replicates.
#' @param bin_width Bin width.
#' @param offset_mode One of \code{"cycle"}, \code{"random"}, or \code{"fixed"}.
#' @param fixed_offset Used when \code{offset_mode = "fixed"}.
#' @param allocation_mode One of \code{"stochastic"}, \code{"fractional"}, or
#'   \code{"deterministic_max"}. Default is \code{"stochastic"}.
#' @param profile Default profile for phase-system matched rows.
#' @param normalization One of \code{"none"}, \code{"site_size"},
#'   \code{"site_duration"}, or \code{"site_size_duration"}.
#' @param variance_method One of \code{"site_bootstrap"} or
#'   \code{"site_jackknife"}.
#' @param bootstrap_B Number of bootstrap resamples if bootstrap is used.
#' @param use_site_metadata_fallback Whether unmatched burials may be allocated
#'   from site metadata bounds.
#' @param augment_to_site_size Whether synthetic burials may be generated up to
#'   \code{site_size}.
#' @param parallel Whether to run replicates in parallel.
#' @param n_cores Number of cores when \code{parallel = TRUE}.
#' @param seed Optional master seed.
#'
#' @return Structured ensemble object.
run_temporal_ensemble <- function(prepared_inputs,
                                  M = 50L,
                                  bin_width = 25,
                                  offset_mode = c("cycle", "random", "fixed"),
                                  fixed_offset = 0,
                                  allocation_mode = c("stochastic", "fractional", "deterministic_max"),
                                  profile = c("trapezoid", "uniform"),
                                  normalization = c("none", "site_size", "site_duration", "site_size_duration"),
                                  variance_method = c("site_bootstrap", "site_jackknife"),
                                  bootstrap_B = 200L,
                                  use_site_metadata_fallback = FALSE,
                                  augment_to_site_size = FALSE,
                                  retain_burial_allocations = FALSE,
                                  parallel = FALSE,
                                  n_cores = 1L,
                                  seed = NULL) {
  offset_mode <- match.arg(offset_mode)
  allocation_mode <- match.arg(allocation_mode)
  profile <- match.arg(profile)
  normalization <- match.arg(normalization)
  variance_method <- match.arg(variance_method)

  .assert_replicate_count(M)
  .assert_offset_mode(offset_mode)
  .assert_parallel_args(parallel, n_cores)

  offsets <- generate_temporal_offsets(
    M = M,
    bin_width = bin_width,
    mode = offset_mode,
    fixed_offset = fixed_offset,
    seed = seed
  )

  rep_seeds <- .derive_replicate_seeds(M = M, seed = seed)

  run_rep <- function(i) {
    run_one_replicate(
      prepared_inputs = prepared_inputs,
      replicate_id = i,
      bin_width = bin_width,
      offset = offsets[i],
      allocation_mode = allocation_mode,
      profile = profile,
      normalization = normalization,
      variance_method = variance_method,
      bootstrap_B = bootstrap_B,
      use_site_metadata_fallback = use_site_metadata_fallback,
      augment_to_site_size = augment_to_site_size,
      seed = rep_seeds[i],
      retain_burial_allocations = retain_burial_allocations
    )
  }

  replicate_results <- NULL

  if (!parallel || M == 1L || n_cores == 1L) {
    replicate_results <- lapply(seq_len(M), run_rep)
  } else {
    # Unix/macOS: mclapply. Windows: PSOCK cluster fallback.
    if (.Platform$OS.type != "windows") {
      replicate_results <- parallel::mclapply(
        X = seq_len(M),
        FUN = run_rep,
        mc.cores = as.integer(n_cores)
      )
    } else {
      cl <- parallel::makeCluster(as.integer(n_cores))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
          library(magrittr)
          library(dplyr)
          library(tidyr)
          library(purrr)
          library(tibble)
          library(rlang)
        })
        NULL
      })

      parallel::clusterExport(
        cl = cl,
        varlist = c(
          "prepared_inputs",
          "bin_width",
          "offsets",
          "allocation_mode",
          "profile",
          "normalization",
          "variance_method",
          "bootstrap_B",
          "use_site_metadata_fallback",
          "augment_to_site_size",
          "rep_seeds",
          "run_one_replicate"
        ),
        envir = environment()
      )

      # Export everything currently sourced into the global environment,
      # including hidden helper functions starting with "."
      needed_funs <- ls(envir = .GlobalEnv, all.names = TRUE)
      parallel::clusterExport(
        cl = cl,
        varlist = needed_funs,
        envir = .GlobalEnv
      )

      replicate_results <- parallel::parLapply(
        cl = cl,
        X = seq_len(M),
        fun = run_rep
      )
    }
  }

  rep_estimates <- dplyr::bind_rows(lapply(replicate_results, `[[`, "replicate_estimates"))
  region_rep_estimates <- dplyr::bind_rows(
    lapply(replicate_results, function(x) x$region_replicate_estimates)
  )

  list(
    replicate_results = replicate_results,
    replicate_estimates = tibble::as_tibble(rep_estimates),
    region_replicate_estimates = tibble::as_tibble(region_rep_estimates),
    settings = list(
      M = M,
      bin_width = bin_width,
      offset_mode = offset_mode,
      fixed_offset = fixed_offset,
      allocation_mode = allocation_mode,
      profile = profile,
      normalization = normalization,
      variance_method = variance_method,
      bootstrap_B = bootstrap_B,
      use_site_metadata_fallback = use_site_metadata_fallback,
      augment_to_site_size = augment_to_site_size,
      parallel = parallel,
      n_cores = as.integer(n_cores),
      seed = seed,
      retain_burial_allocations = retain_burial_allocations
    ),
    diagnostics = list(
      offsets = offsets,
      n_replicates = length(replicate_results),
      n_estimate_rows = nrow(rep_estimates)
    )
  )
}

# -------------------------------------------------------------------------
# Final pooling
# -------------------------------------------------------------------------

#' Finalize a temporal ensemble by pooling replicate uncertainty
#'
#' @param ensemble_result Output from run_temporal_ensemble().
#' @param conf_level Confidence level for intervals.
#' @param use_t Whether to use t critical values in pooling.
#' @param canonical_offset Offset for the final pooled bin grid.
#'
#' @return Structured pooled result object.
finalize_ensemble <- function(ensemble_result,
                              conf_level = 0.95,
                              use_t = TRUE,
                              canonical_offset = 0) {
  if (!is.list(ensemble_result) || is.null(ensemble_result$replicate_estimates)) {
    rlang::abort("`ensemble_result` must be an object returned by `run_temporal_ensemble()`.")
  }

  rep_est_raw <- tibble::as_tibble(ensemble_result$replicate_estimates)

  if (nrow(rep_est_raw) == 0) {
    rlang::abort("`ensemble_result$replicate_estimates` is empty.")
  }

  bin_width <- ensemble_result$settings$bin_width
  if (is.null(bin_width) || !is.numeric(bin_width) || length(bin_width) != 1L || is.na(bin_width) || bin_width <= 0) {
    rlang::abort("Could not determine a valid `bin_width` from `ensemble_result$settings`.")
  }

  canonical_grid <- make_canonical_grid_from_replicates(
    replicate_estimates = rep_est_raw,
    bin_width = bin_width,
    canonical_offset = canonical_offset
  )

  split_reps <- split(rep_est_raw, rep_est_raw$replicate_id)

  rep_est_regridded <- dplyr::bind_rows(
    lapply(split_reps, regrid_replicate_estimates, canonical_grid = canonical_grid)
  )

  pooled <- pool_replicate_uncertainty(
    replicate_estimates = rep_est_regridded,
    conf_level = conf_level,
    use_t = use_t
  )

  components <- summarize_uncertainty_components(pooled)

  # Regional pooling
  pooled_region <- NULL
  region_rep_raw <- NULL
  region_rep_regridded <- NULL

  if (!is.null(ensemble_result$region_replicate_estimates)) {
    region_rep_raw <- tibble::as_tibble(ensemble_result$region_replicate_estimates)

    if (nrow(region_rep_raw) > 0) {
      split_region_reps <- split(region_rep_raw, region_rep_raw$replicate_id)

      region_rep_regridded <- dplyr::bind_rows(
        lapply(split_region_reps, regrid_region_replicate_estimates, canonical_grid = canonical_grid)
      )

      pooled_region <- region_rep_regridded %>%
        dplyr::group_by(site_region, horizon_bin, bin_start, bin_end) %>%
        dplyr::summarise(
          M = dplyr::n_distinct(replicate_id),
          estimate_mean = mean(estimate, na.rm = TRUE),
          W = mean(var_hat, na.rm = TRUE),
          B = if (dplyr::n() > 1) stats::var(estimate, na.rm = TRUE) else 0,
          .groups = "drop"
        ) %>%
        dplyr::rename(estimate = estimate_mean)

      pooled_region$W[!is.finite(pooled_region$W)] <- NA_real_
      pooled_region$B[!is.finite(pooled_region$B)] <- 0
      pooled_region$T <- pooled_region$W + (1 + 1 / pooled_region$M) * pooled_region$B
      pooled_region$se_within <- sqrt(pooled_region$W)
      pooled_region$se_between <- sqrt(pooled_region$B)
      pooled_region$se_total <- sqrt(pooled_region$T)

      if (use_t) {
        pooled_region$df <- Inf
        use_df <- pooled_region$M > 1 & is.finite(pooled_region$W) & is.finite(pooled_region$B) & pooled_region$B > 0
        pooled_region$df[use_df] <- (pooled_region$M[use_df] - 1) *
          (1 + pooled_region$W[use_df] / ((1 + 1 / pooled_region$M[use_df]) * pooled_region$B[use_df]))^2

        crit <- stats::qt((1 + conf_level) / 2, df = pooled_region$df)
        crit[!is.finite(crit)] <- stats::qnorm((1 + conf_level) / 2)
      } else {
        pooled_region$df <- Inf
        crit <- rep(stats::qnorm((1 + conf_level) / 2), nrow(pooled_region))
      }

      pooled_region$crit <- crit
      pooled_region$conf_level <- conf_level
      pooled_region$lower <- pooled_region$estimate - pooled_region$crit * pooled_region$se_total
      pooled_region$upper <- pooled_region$estimate + pooled_region$crit * pooled_region$se_total
    }
  }

  # Collect replicate-level data for later significance testing.
  site_bin_replicates <- lapply(
    ensemble_result$replicate_results,
    function(x) x$analysis$data$site_bin
  )

  bin_curve_replicates <- lapply(
    ensemble_result$replicate_results,
    function(x) x$analysis$data$bin_curve
  )

  region_bin_replicates <- lapply(
    ensemble_result$replicate_results,
    function(x) x$analysis$data$region_bin
  )

  burial_allocation_replicates <- lapply(
    ensemble_result$replicate_results,
    function(x) x$burial_allocations
  )

  if (all(vapply(burial_allocation_replicates, is.null, logical(1)))) {
    burial_allocation_replicates <- NULL
  }

  list(
    pooled = list(
      estimates = pooled,
      components = components,
      region_estimates = if (is.null(pooled_region)) NULL else tibble::as_tibble(pooled_region),
      canonical_grid = canonical_grid
    ),
    replicate_estimates = rep_est_regridded,
    raw_replicate_estimates = rep_est_raw,
    region_replicate_estimates = if (is.null(region_rep_regridded)) NULL else tibble::as_tibble(region_rep_regridded),
    raw_region_replicate_estimates = if (is.null(region_rep_raw)) NULL else tibble::as_tibble(region_rep_raw),
    replicate_results = ensemble_result$replicate_results,
    replicate_data = list(
      site_bin = site_bin_replicates,
      bin_curve = bin_curve_replicates,
      region_bin = region_bin_replicates,
      burial_allocations = burial_allocation_replicates
    ),
    diagnostics = list(
      ensemble = ensemble_result$diagnostics,
      pooling = list(
        n_replicates = dplyr::n_distinct(rep_est_regridded$replicate_id),
        n_pooled_bins = nrow(pooled),
        n_pooled_region_rows = if (is.null(pooled_region)) 0L else nrow(pooled_region),
        conf_level = conf_level,
        use_t = use_t,
        bin_width = bin_width,
        canonical_offset = canonical_offset
      )
    ),
    settings = ensemble_result$settings
  )
}
