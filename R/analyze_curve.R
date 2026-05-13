# analyze_curve.R
#
# Analytical summaries for OccuPast.
#
# Expected allocation input schema:
#   UID, input_type, record_id, burial_id, unwrap_index,
#   site_id, phase_id, system_name, phase_name,
#   sex_gender, age, chronology_source, is_synthetic,
#   horizon_bin, bin_start, bin_end,
#   overlap_len, raw_weight, alloc_weight,
#   allocation_mode, profile, chosen_bin, chosen_prob
#
# Optional site metadata schema:
#   site_id, site_name, coord_y, coord_x, site_country, site_region,
#   site_admin, site_start, site_end, site_size, site_dig_date

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   magrittr
#
# Suggested namespace usage:
#   @importFrom dplyr summarise group_by ungroup left_join mutate n_distinct
#   @importFrom tibble as_tibble
#   @importFrom rlang abort
#   @importFrom magrittr %>%

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.extract_allocation_tbl <- function(x) {
  if (is.list(x) && !inherits(x, "data.frame") && "data" %in% names(x)) {
    out <- tibble::as_tibble(x$data)
  } else {
    out <- tibble::as_tibble(x)
  }

  if (exists(".standardize_bin_columns", mode = "function")) {
    out <- .standardize_bin_columns(out, arg = "allocation_tbl")
  }

  out
}

.extract_site_metadata_tbl <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !inherits(x, "data.frame") && "data" %in% names(x)) {
    return(tibble::as_tibble(x$data))
  }
  tibble::as_tibble(x)
}

.require_allocation_table <- function(allocation_tbl) {
  validate_required_fields(
    allocation_tbl,
    c(
      "UID", "input_type", "record_id", "burial_id", "unwrap_index",
      "site_id", "phase_id", "system_name", "phase_name",
      "sex_gender", "age", "chronology_source", "is_synthetic",
      "horizon_bin", "bin_start", "bin_end",
      "overlap_len", "raw_weight", "alloc_weight",
      "allocation_mode", "profile", "chosen_bin", "chosen_prob"
    ),
    "allocation_tbl"
  )
  invisible(TRUE)
}

.require_site_metadata_analysis <- function(site_metadata) {
  validate_required_fields(
    site_metadata,
    c(
      "site_id", "site_name", "coord_y", "coord_x", "site_country",
      "site_region", "site_admin", "site_start", "site_end",
      "site_size", "site_dig_date"
    ),
    "site_metadata"
  )
  invisible(TRUE)
}

.safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

.compute_site_duration <- function(site_start, site_end) {
  s <- .safe_numeric(site_start)
  e <- .safe_numeric(site_end)
  out <- e - s
  out[!is.finite(out) | out <= 0] <- NA_real_
  out
}

# -------------------------------------------------------------------------
# Site-bin construction
# -------------------------------------------------------------------------

#' Build site-by-bin contributions
#' Useful for analysing the reliance of the dataset on one site
#'
#' @param allocation_tbl Allocation output from harmonize_chronologies_merged().
#' @param site_metadata Optional prepared site metadata.
#' @param attach_metadata If TRUE, join site metadata onto the site-bin table.
#'
#' @return Tibble with one row per site per bin.
build_site_bin <- function(allocation_tbl,
                           site_metadata = NULL,
                           attach_metadata = TRUE) {
  allocation_tbl <- tibble::as_tibble(allocation_tbl)
  .require_allocation_table(allocation_tbl)

  sb <- allocation_tbl %>%
    dplyr::group_by(site_id, horizon_bin, bin_start, bin_end) %>%
    dplyr::summarise(
      value = sum(alloc_weight, na.rm = TRUE),

      n_alloc_rows = dplyr::n(),
      n_UID = dplyr::n_distinct(UID),
      n_real_UID = dplyr::n_distinct(UID[!is_synthetic]),
      n_synthetic_UID = dplyr::n_distinct(UID[is_synthetic]),

      value_real = sum(alloc_weight[!is_synthetic], na.rm = TRUE),
      value_synthetic = sum(alloc_weight[is_synthetic], na.rm = TRUE),

      value_phase_system = sum(alloc_weight[chronology_source == "phase_system"], na.rm = TRUE),
      value_site_fallback = sum(alloc_weight[chronology_source == "site_metadata_fallback"], na.rm = TRUE),

      n_phase_system_UID = dplyr::n_distinct(UID[chronology_source == "phase_system"]),
      n_site_fallback_UID = dplyr::n_distinct(UID[chronology_source == "site_metadata_fallback"]),

      .groups = "drop"
    )

  if (!is.null(site_metadata) && isTRUE(attach_metadata)) {
    site_metadata <- tibble::as_tibble(site_metadata)
    .require_site_metadata_analysis(site_metadata)

    meta_small <- site_metadata %>%
      dplyr::select(
        site_id, site_name, coord_y, coord_x, site_country, site_region,
        site_admin, site_start, site_end, site_size, site_dig_date
      )

    sb <- dplyr::left_join(sb, meta_small, by = "site_id")
  }

  tibble::as_tibble(sb)
}

# -------------------------------------------------------------------------
# Normalization
# -------------------------------------------------------------------------

#' Normalize site-by-bin contributions
#'
#' @param site_bin Site-by-bin table from build_site_bin().
#' @param mode One of "none", "site_size", "site_duration", "site_size_duration".
#'
#' @return Site-by-bin tibble with normalization columns added.
normalize_site_bin <- function(site_bin,
                               mode = c("none", "site_size", "site_duration", "site_size_duration")) {
  mode <- match.arg(mode)
  site_bin <- tibble::as_tibble(site_bin)

  validate_required_fields(
    site_bin,
    c("site_id", "horizon_bin", "value"),
    "site_bin"
  )

  out <- site_bin

  if (!"site_size" %in% names(out)) {
    out$site_size <- NA
  }
  if (!"site_start" %in% names(out)) {
    out$site_start <- NA
  }
  if (!"site_end" %in% names(out)) {
    out$site_end <- NA
  }

  out$site_size_num <- .safe_numeric(out$site_size)
  out$site_size_num[!is.finite(out$site_size_num) | out$site_size_num <= 0] <- NA_real_

  out$site_duration <- .compute_site_duration(out$site_start, out$site_end)

  # Normalisation turns allocated counts into comparable intensities.
  # Depending on mode: Y*,sj = Y_sj, Y_sj/S_s, Y_sj/D_s, or Y_sj/(S_s D_s).
  denom <- switch(
    mode,
    none = rep(1, nrow(out)),
    site_size = out$site_size_num,
    site_duration = out$site_duration,
    site_size_duration = out$site_size_num * out$site_duration
  )

  if (mode == "none") {
    out$value_norm <- out$value
    out$value_real_norm <- if ("value_real" %in% names(out)) out$value_real else NA_real_
    out$value_synthetic_norm <- if ("value_synthetic" %in% names(out)) out$value_synthetic else NA_real_
    out$value_phase_system_norm <- if ("value_phase_system" %in% names(out)) out$value_phase_system else NA_real_
    out$value_site_fallback_norm <- if ("value_site_fallback" %in% names(out)) out$value_site_fallback else NA_real_
  } else {
    denom_bad <- !is.finite(denom) | denom <= 0
    denom[denom_bad] <- NA_real_

    out$value_norm <- out$value / denom
    out$value_real_norm <- if ("value_real" %in% names(out)) out$value_real / denom else NA_real_
    out$value_synthetic_norm <- if ("value_synthetic" %in% names(out)) out$value_synthetic / denom else NA_real_
    out$value_phase_system_norm <- if ("value_phase_system" %in% names(out)) out$value_phase_system / denom else NA_real_
    out$value_site_fallback_norm <- if ("value_site_fallback" %in% names(out)) out$value_site_fallback / denom else NA_real_
  }

  out$normalization_mode <- mode

  tibble::as_tibble(out)
}

# -------------------------------------------------------------------------
# Overall bin summaries
# -------------------------------------------------------------------------

#' Summarize the overall bin curve
#'
#' @param site_bin Normalized or unnormalized site-by-bin table.
#' @param value_col Contribution column to summarize.
#'
#' @return One row per bin.
summarize_bin_curve <- function(site_bin,
                                value_col = "value_norm") {
  site_bin <- tibble::as_tibble(site_bin)

  validate_required_fields(
    site_bin,
    c("site_id", "horizon_bin", "bin_start", "bin_end", value_col),
    "site_bin"
  )

  has_real <- "value_real_norm" %in% names(site_bin)
  has_synth <- "value_synthetic_norm" %in% names(site_bin)
  has_phase <- "value_phase_system_norm" %in% names(site_bin)
  has_fallback <- "value_site_fallback_norm" %in% names(site_bin)

  if (!has_real) site_bin$value_real_norm <- NA_real_
  if (!has_synth) site_bin$value_synthetic_norm <- NA_real_
  if (!has_phase) site_bin$value_phase_system_norm <- NA_real_
  if (!has_fallback) site_bin$value_site_fallback_norm <- NA_real_

  # The estimand here is an equal-site-weighted mean per temporal bin.
  out <- site_bin %>%
    dplyr::group_by(horizon_bin, bin_start, bin_end) %>%
    dplyr::summarise(
      total_value = sum(.data[[value_col]], na.rm = TRUE),
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      median_value = stats::median(.data[[value_col]], na.rm = TRUE),
      n_sites_contributing = sum(!is.na(.data[[value_col]]) & .data[[value_col]] > 0),
      n_sites_observed = dplyr::n_distinct(site_id),

      total_real = sum(value_real_norm, na.rm = TRUE),
      total_synthetic = sum(value_synthetic_norm, na.rm = TRUE),
      total_phase_system = sum(value_phase_system_norm, na.rm = TRUE),
      total_site_fallback = sum(value_site_fallback_norm, na.rm = TRUE),

      .groups = "drop"
    )

  tibble::as_tibble(out)
}

# -------------------------------------------------------------------------
# Regional summaries
# -------------------------------------------------------------------------

#' Summarize region-by-bin curves
#'
#' @param site_bin Normalized or unnormalized site-by-bin table.
#' @param value_col Contribution column to summarize.
#' @param region_col Region column name.
#'
#' @return One row per region per bin.
summarize_region_bin <- function(site_bin,
                                 value_col = "value_norm",
                                 region_col = "site_region") {
  site_bin <- tibble::as_tibble(site_bin)

  validate_required_fields(
    site_bin,
    c("site_id", "horizon_bin", "bin_start", "bin_end", value_col, region_col),
    "site_bin"
  )

  has_real <- "value_real_norm" %in% names(site_bin)
  has_synth <- "value_synthetic_norm" %in% names(site_bin)
  has_phase <- "value_phase_system_norm" %in% names(site_bin)
  has_fallback <- "value_site_fallback_norm" %in% names(site_bin)

  if (!has_real) site_bin$value_real_norm <- NA_real_
  if (!has_synth) site_bin$value_synthetic_norm <- NA_real_
  if (!has_phase) site_bin$value_phase_system_norm <- NA_real_
  if (!has_fallback) site_bin$value_site_fallback_norm <- NA_real_

  out <- site_bin %>%
    dplyr::filter(!is.na(.data[[region_col]]) & .data[[region_col]] != "") %>%
    dplyr::group_by(.data[[region_col]], horizon_bin, bin_start, bin_end) %>%
    dplyr::summarise(
      total_value = sum(.data[[value_col]], na.rm = TRUE),
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      median_value = stats::median(.data[[value_col]], na.rm = TRUE),
      n_sites_contributing = sum(!is.na(.data[[value_col]]) & .data[[value_col]] > 0),
      n_sites_observed = dplyr::n_distinct(site_id),

      total_real = sum(value_real_norm, na.rm = TRUE),
      total_synthetic = sum(value_synthetic_norm, na.rm = TRUE),
      total_phase_system = sum(value_phase_system_norm, na.rm = TRUE),
      total_site_fallback = sum(value_site_fallback_norm, na.rm = TRUE),

      .groups = "drop"
    )

  names(out)[1] <- region_col
  tibble::as_tibble(out)
}

# -------------------------------------------------------------------------
# Main analysis wrapper
# -------------------------------------------------------------------------

#' Analyze occupancy curves from temporal allocation output
#'
#' @param allocation_result Result object from harmonize_chronologies_merged(), or
#'   a raw allocation tibble.
#' @param site_metadata Optional prepared site metadata, or NULL.
#' @param normalization One of "none", "site_size", "site_duration", "site_size_duration".
#' @param region_col Region field to use from metadata/site_bin.
#'
#' @return Structured list with site-bin, overall curve, and regional curves.
analyze_occupancy_curve <- function(allocation_result,
                                    site_metadata = NULL,
                                    normalization = c("none", "site_size", "site_duration", "site_size_duration"),
                                    region_col = "site_region") {
  normalization <- match.arg(normalization)

  allocation_tbl <- .extract_allocation_tbl(allocation_result)
  .require_allocation_table(allocation_tbl)

  site_metadata <- .extract_site_metadata_tbl(site_metadata)
  if (!is.null(site_metadata)) {
    .require_site_metadata_analysis(site_metadata)
  }

  site_bin_raw <- build_site_bin(
    allocation_tbl = allocation_tbl,
    site_metadata = site_metadata,
    attach_metadata = TRUE
  )

  site_bin <- normalize_site_bin(
    site_bin = site_bin_raw,
    mode = normalization
  )

  value_col <- if (normalization == "none") "value" else "value_norm"

  bin_curve <- summarize_bin_curve(
    site_bin = site_bin,
    value_col = value_col
  )

  region_bin <- NULL
  if (region_col %in% names(site_bin)) {
    region_bin <- summarize_region_bin(
      site_bin = site_bin,
      value_col = value_col,
      region_col = region_col
    )
  }

  diagnostics <- list(
    n_allocation_rows = nrow(allocation_tbl),
    n_site_bin_rows = nrow(site_bin),
    n_bin_rows = nrow(bin_curve),
    n_region_bin_rows = if (is.null(region_bin)) 0L else nrow(region_bin),
    n_unique_sites = dplyr::n_distinct(site_bin$site_id),
    n_unique_bins = dplyr::n_distinct(site_bin$horizon_bin),
    n_unique_regions = if (!is.null(region_bin)) dplyr::n_distinct(region_bin[[region_col]]) else 0L,
    normalization = normalization,
    region_col = region_col,
    total_allocated_weight = sum(allocation_tbl$alloc_weight, na.rm = TRUE),
    total_site_value = sum(site_bin[[value_col]], na.rm = TRUE)
  )

  list(
    data = list(
      site_bin = site_bin,
      bin_curve = bin_curve,
      region_bin = region_bin
    ),
    diagnostics = diagnostics,
    settings = list(
      normalization = normalization,
      region_col = region_col,
      value_col = value_col
    )
  )
}


# -------------------------------------------------------------------------
# Legacy aliases for earlier "bucket" terminology.
# -------------------------------------------------------------------------

build_site_bucket <- build_site_bin
normalize_site_bucket <- normalize_site_bin
summarize_bucket_curve <- summarize_bin_curve
summarize_region_bucket <- summarize_region_bin
