# allocate_time.R
#
# Temporal allocation helpers for the merged occupancy package.
#
# Assumes preparation has already been completed with prepare_data.R, so that:
# - chronology is a strict phase-system table
# - mortuary data are already burial-level and include UID
# - chronology linkage happens by phase_id
# - site metadata may optionally provide fallback bounds and site_size
#
# Canonical input schemas expected here:
#
# chronology:
#   phase_id, system_name, phase_name,
#   horizon_start, horizon_end, fade_in_start, fade_out_end
#
# burial-level mortuary:
#   UID, input_type, record_id, burial_id, unwrap_index,
#   site_id, phase_id, sex_gender, age
#
# optional site metadata:
#   site_id, site_name, coord_y, coord_x, site_country, site_region,
#   site_admin, site_start, site_end, site_size, site_dig_date
#
# Main design choices:
# - buckets are half-open intervals [bucket_start, bucket_end)
# - default chronology profile is trapezoid:
#     fade_in_start -> horizon_start : linear rise
#     horizon_start -> horizon_end   : plateau
#     horizon_end   -> fade_out_end  : linear fall
# - site metadata fallback uses a simple uniform profile over [site_start, site_end)
# - allocation modes:
#     fractional        : keep normalized weights across all overlapping buckets
#     stochastic        : draw one bucket using normalized weights
#     deterministic_max : choose the single bucket with maximum normalized weight

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.require_allocation_inputs <- function(mortuary, chronology) {
  validate_required_fields(
    mortuary,
    c("UID", "input_type", "record_id", "burial_id", "unwrap_index",
      "site_id", "phase_id", "sex_gender", "age"),
    "mortuary"
  )

  validate_required_fields(
    chronology,
    c("phase_id", "system_name", "phase_name",
      "horizon_start", "horizon_end", "fade_in_start", "fade_out_end"),
    "chronology"
  )

  invisible(TRUE)
}

.require_site_metadata_inputs <- function(site_metadata) {
  validate_required_fields(
    site_metadata,
    c("site_id", "site_start", "site_end", "site_size"),
    "site_metadata"
  )
  invisible(TRUE)
}

.empty_allocation_table <- function() {
  tibble::tibble(
    UID = character(),
    input_type = character(),
    record_id = numeric(),
    burial_id = character(),
    unwrap_index = integer(),
    site_id = character(),
    phase_id = character(),
    system_name = character(),
    phase_name = character(),
    sex_gender = character(),
    age = character(),
    chronology_source = character(),
    is_synthetic = logical(),
    horizon_bucket = numeric(),
    bucket_start = numeric(),
    bucket_end = numeric(),
    overlap_len = numeric(),
    raw_weight = numeric(),
    alloc_weight = numeric(),
    allocation_mode = character(),
    profile = character(),
    chosen_bucket = logical(),
    chosen_prob = numeric()
  )
}

.assert_positive_scalar <- function(x, arg = "x") {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0) {
    rlang::abort(sprintf("`%s` must be a single positive number.", arg))
  }
}

.safe_max_numeric <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  max(x)
}

# -------------------------------------------------------------------------
# Bucket helpers
# -------------------------------------------------------------------------

make_bucket_grid <- function(year_min, year_max,
                             bucket_width = 25,
                             offset = 0) {
  .assert_positive_scalar(bucket_width, "bucket_width")

  if (!is.numeric(year_min) || !is.numeric(year_max) ||
      length(year_min) != 1L || length(year_max) != 1L ||
      is.na(year_min) || is.na(year_max)) {
    rlang::abort("`year_min` and `year_max` must be single finite numeric values.")
  }

  if (year_max <= year_min) {
    rlang::abort("`year_max` must be greater than `year_min`.")
  }

  start0 <- floor((year_min - offset) / bucket_width) * bucket_width + offset
  end0   <- ceiling((year_max - offset) / bucket_width) * bucket_width + offset

  starts <- seq(start0, end0 - bucket_width, by = bucket_width)

  tibble::tibble(
    horizon_bucket = starts,
    bucket_start = starts,
    bucket_end = starts + bucket_width
  )
}

bucket_seq_with_offset <- function(s, e, width, offset = 0, eps = 1e-9) {
  if (any(is.na(c(s, e))) || e <= s) return(numeric(0))
  b0 <- floor((s - offset) / width) * width + offset
  b1 <- floor(((e - eps) - offset) / width) * width + offset
  seq(b0, b1, by = width)
}

# -------------------------------------------------------------------------
# Profile helpers
# -------------------------------------------------------------------------

.profile_value <- function(x,
                           fade_in_start,
                           horizon_start,
                           horizon_end,
                           fade_out_end,
                           profile = c("trapezoid", "uniform")) {
  profile <- match.arg(profile)

  if (!is.finite(x)) return(0)

  fi <- fade_in_start
  hs <- horizon_start
  he <- horizon_end
  fo <- fade_out_end

  if (profile == "uniform") {
    return(ifelse(x >= fi && x < fo, 1, 0))
  }

  if (!is.finite(fi) || !is.finite(hs) || !is.finite(he) || !is.finite(fo)) return(0)
  if (x < fi || x >= fo) return(0)

  if (x < hs) {
    if (hs <= fi) return(1)
    return((x - fi) / (hs - fi))
  }

  if (x < he) {
    return(1)
  }

  if (fo <= he) return(1)
  (fo - x) / (fo - he)
}

.integrate_profile_interval <- function(a, b,
                                        fade_in_start,
                                        horizon_start,
                                        horizon_end,
                                        fade_out_end,
                                        profile = c("trapezoid", "uniform")) {
  profile <- match.arg(profile)

  if (!is.finite(a) || !is.finite(b) || b <= a) return(0)

  fi <- fade_in_start
  hs <- horizon_start
  he <- horizon_end
  fo <- fade_out_end

  if (profile == "uniform") {
    lo <- max(a, fi)
    hi <- min(b, fo)
    return(max(0, hi - lo))
  }

  cuts <- sort(unique(c(a, b, fi, hs, he, fo)))
  cuts <- cuts[cuts >= a & cuts <= b]

  if (length(cuts) < 2L) return(0)

  total <- 0
  for (i in seq_len(length(cuts) - 1L)) {
    x0 <- cuts[i]
    x1 <- cuts[i + 1L]
    if (x1 <= x0) next

    y0 <- .profile_value(x0, fi, hs, he, fo, profile = "trapezoid")
    y1 <- .profile_value(x1, fi, hs, he, fo, profile = "trapezoid")

    total <- total + 0.5 * (y0 + y1) * (x1 - x0)
  }

  total
}

# -------------------------------------------------------------------------
# Overlap / weighting
# -------------------------------------------------------------------------

weighted_overlaps <- function(fade_in_start,
                              horizon_start,
                              horizon_end,
                              fade_out_end,
                              bucket_width = 25,
                              offset = 0,
                              profile = c("trapezoid", "uniform")) {
  profile <- match.arg(profile)
  .assert_positive_scalar(bucket_width, "bucket_width")

  eff_start <- fade_in_start
  eff_end <- fade_out_end

  if (any(is.na(c(eff_start, eff_end))) || eff_end <= eff_start) {
    return(tibble::tibble(
      horizon_bucket = numeric(),
      bucket_start = numeric(),
      bucket_end = numeric(),
      overlap_len = numeric(),
      raw_weight = numeric(),
      alloc_weight = numeric()
    ))
  }

  bks <- bucket_seq_with_offset(eff_start, eff_end, bucket_width, offset)
  if (!length(bks)) {
    return(tibble::tibble(
      horizon_bucket = numeric(),
      bucket_start = numeric(),
      bucket_end = numeric(),
      overlap_len = numeric(),
      raw_weight = numeric(),
      alloc_weight = numeric()
    ))
  }

  bucket_start <- bks
  bucket_end <- bks + bucket_width

  overlap_len <- vapply(
    seq_along(bucket_start),
    function(i) {
      max(0, min(eff_end, bucket_end[i]) - max(eff_start, bucket_start[i]))
    },
    numeric(1)
  )

  raw_weight <- vapply(
    seq_along(bucket_start),
    function(i) {
      .integrate_profile_interval(
        a = bucket_start[i],
        b = bucket_end[i],
        fade_in_start = fade_in_start,
        horizon_start = horizon_start,
        horizon_end = horizon_end,
        fade_out_end = fade_out_end,
        profile = profile
      )
    },
    numeric(1)
  )

  keep <- overlap_len > 0 & raw_weight > 0
  if (!any(keep)) {
    return(tibble::tibble(
      horizon_bucket = numeric(),
      bucket_start = numeric(),
      bucket_end = numeric(),
      overlap_len = numeric(),
      raw_weight = numeric(),
      alloc_weight = numeric()
    ))
  }

  out <- tibble::tibble(
    horizon_bucket = bucket_start[keep],
    bucket_start = bucket_start[keep],
    bucket_end = bucket_end[keep],
    overlap_len = overlap_len[keep],
    raw_weight = raw_weight[keep]
  )

  total_w <- sum(out$raw_weight)
  if (!is.finite(total_w) || total_w <= 0) {
    out$alloc_weight <- 1 / nrow(out)
  } else {
    out$alloc_weight <- out$raw_weight / total_w
  }

  out
}

# -------------------------------------------------------------------------
# Allocation helpers
# -------------------------------------------------------------------------

choose_one_bucket <- function(weights_tbl,
                              seed = NULL,
                              method = c("stochastic", "deterministic_max")) {
  method <- match.arg(method)

  if (nrow(weights_tbl) == 0) {
    return(tibble::tibble(
      horizon_bucket = numeric(),
      bucket_start = numeric(),
      bucket_end = numeric(),
      overlap_len = numeric(),
      raw_weight = numeric(),
      alloc_weight = numeric(),
      chosen_bucket = logical(),
      chosen_prob = numeric()
    ))
  }

  if (method == "deterministic_max") {
    idx <- which.max(weights_tbl$alloc_weight)
  } else {
    if (!is.null(seed)) set.seed(seed)
    idx <- sample.int(
      n = nrow(weights_tbl),
      size = 1L,
      prob = weights_tbl$alloc_weight
    )
  }

  out <- weights_tbl[idx, , drop = FALSE]
  out$chosen_bucket <- TRUE
  out$chosen_prob <- out$alloc_weight
  out$alloc_weight <- 1
  out
}

#' Allocate one burial-level row across time buckets
#'
#' @description
#' Allocates a single burial-level record to one or more temporal buckets using
#' the chronology bounds attached to its phase or fallback metadata interval.
#'
#' The temporal weighting profile is normally trapezoidal at the phase level:
#' \itemize{
#'   \item \code{fade_in_start -> horizon_start}: linear rise
#'   \item \code{horizon_start -> horizon_end}: plateau
#'   \item \code{horizon_end -> fade_out_end}: linear fall
#' }
#'
#' Under \code{allocation_mode = "stochastic"} (the default), the burial is
#' assigned to exactly one bucket by random draw from the normalized bucket
#' weights implied by the temporal profile. This makes each replicate a plausible
#' temporal realization rather than a deterministic expected-value curve.
#'
#' Under \code{allocation_mode = "fractional"}, the burial contributes
#' proportionally across all overlapping buckets. Under
#' \code{allocation_mode = "deterministic_max"}, it is assigned to the bucket
#' with maximal weight.
#'
#' @param row_one Single joined mortuary+chronology row as a one-row data.frame/tibble.
#' @param bucket_width Width of the temporal buckets.
#' @param offset Bucket grid offset.
#' @param allocation_mode One of \code{"stochastic"}, \code{"fractional"}, or
#'   \code{"deterministic_max"}. Default is \code{"stochastic"}.
#' @param profile One of \code{"trapezoid"} or \code{"uniform"}.
#' @param seed Optional seed for stochastic allocation.
#'
#' @return Tibble with one or more allocated bucket rows.
allocate_one_row <- function(row_one,
                             bucket_width = 25,
                             offset = 0,
                             allocation_mode = c("stochastic", "fractional", "deterministic_max"),
                             profile = c("trapezoid", "uniform"),
                             seed = NULL) {
  allocation_mode <- match.arg(allocation_mode)
  profile <- match.arg(profile)

  weights_tbl <- weighted_overlaps(
    fade_in_start = row_one$fade_in_start[[1]],
    horizon_start = row_one$horizon_start[[1]],
    horizon_end = row_one$horizon_end[[1]],
    fade_out_end = row_one$fade_out_end[[1]],
    bucket_width = bucket_width,
    offset = offset,
    profile = profile
  )

  if (nrow(weights_tbl) == 0) {
    return(.empty_allocation_table())
  }

  if (allocation_mode == "fractional") {
    weights_tbl$chosen_bucket <- FALSE
    weights_tbl$chosen_prob <- NA_real_
    alloc_tbl <- weights_tbl
  } else {
    alloc_tbl <- choose_one_bucket(
      weights_tbl = weights_tbl,
      seed = seed,
      method = allocation_mode
    )
  }

  core <- tibble::tibble(
    UID = row_one$UID[[1]],
    input_type = row_one$input_type[[1]],
    record_id = row_one$record_id[[1]],
    burial_id = row_one$burial_id[[1]],
    unwrap_index = row_one$unwrap_index[[1]],
    site_id = row_one$site_id[[1]],
    phase_id = row_one$phase_id[[1]],
    system_name = row_one$system_name[[1]],
    phase_name = row_one$phase_name[[1]],
    sex_gender = row_one$sex_gender[[1]],
    age = row_one$age[[1]],
    chronology_source = row_one$chronology_source[[1]],
    is_synthetic = row_one$is_synthetic[[1]]
  )

  core[rep(1, nrow(alloc_tbl)), , drop = FALSE] |>
    tibble::as_tibble() |>
    dplyr::bind_cols(
      alloc_tbl,
      tibble::tibble(
        allocation_mode = allocation_mode,
        profile = profile
      )
    ) |>
    dplyr::select(
      UID, input_type, record_id, burial_id, unwrap_index,
      site_id, phase_id, system_name, phase_name,
      sex_gender, age, chronology_source, is_synthetic,
      horizon_bucket, bucket_start, bucket_end,
      overlap_len, raw_weight, alloc_weight,
      allocation_mode, profile,
      chosen_bucket, chosen_prob
    )
}

# -------------------------------------------------------------------------
# Fallback / augmentation resolver
# -------------------------------------------------------------------------

resolve_allocation_bounds <- function(mortuary,
                                      chronology,
                                      site_metadata = NULL,
                                      use_site_metadata_fallback = FALSE,
                                      augment_to_site_size = FALSE,
                                      synthetic_prefix = "SYN") {
  mortuary <- tibble::as_tibble(mortuary)
  chronology <- tibble::as_tibble(chronology)
  .require_allocation_inputs(mortuary, chronology)

  if (!is.null(site_metadata)) {
    site_metadata <- tibble::as_tibble(site_metadata)
    .require_site_metadata_inputs(site_metadata)
  }

  phase_dups <- chronology$phase_id[duplicated(chronology$phase_id)]
  if (length(phase_dups) > 0) {
    rlang::abort(
      paste0(
        "Chronology contains duplicate `phase_id` values: ",
        paste(unique(phase_dups), collapse = ", ")
      )
    )
  }

  log_tbl <- log_exclusions_init()

  joined <- dplyr::left_join(mortuary, chronology, by = "phase_id")
  joined$chronology_source <- "phase_system"
  joined$is_synthetic <- FALSE

  incomplete_idx <- which(
    is.na(joined$system_name) |
      is.na(joined$phase_name) |
      is.na(joined$horizon_start) |
      is.na(joined$horizon_end) |
      is.na(joined$fade_in_start) |
      is.na(joined$fade_out_end)
  )

  if (length(incomplete_idx) > 0 && use_site_metadata_fallback) {
    if (is.null(site_metadata)) {
      warning(
        "Site metadata fallback was requested but no `site_metadata` was supplied. Unmatched burials will be dropped.",
        call. = FALSE
      )
    } else {
      meta_small <- site_metadata[, c("site_id", "site_start", "site_end"), drop = FALSE]
      meta_small <- meta_small[!duplicated(meta_small$site_id), , drop = FALSE]

      joined <- dplyr::left_join(
        joined,
        meta_small,
        by = "site_id"
      )

      fallback_ok <- !is.na(joined$site_start) & !is.na(joined$site_end)
      rescue_idx <- incomplete_idx[fallback_ok[incomplete_idx]]

      if (length(rescue_idx) > 0) {
        joined$system_name[rescue_idx] <- "site_metadata_fallback"
        joined$phase_name[rescue_idx] <- "site_metadata_fallback"
        joined$fade_in_start[rescue_idx] <- suppressWarnings(as.numeric(joined$site_start[rescue_idx]))
        joined$horizon_start[rescue_idx] <- suppressWarnings(as.numeric(joined$site_start[rescue_idx]))
        joined$horizon_end[rescue_idx] <- suppressWarnings(as.numeric(joined$site_end[rescue_idx]))
        joined$fade_out_end[rescue_idx] <- suppressWarnings(as.numeric(joined$site_end[rescue_idx]))
        joined$chronology_source[rescue_idx] <- "site_metadata_fallback"
      }

      for (pid in unique(joined$phase_id[rescue_idx])) {
        n_pid <- sum(joined$phase_id[rescue_idx] == pid, na.rm = TRUE)
        log_tbl <- log_exclusions_add(
          log = log_tbl,
          stage = "resolve_allocation_bounds",
          entity_type = "phase",
          entity_id = as.character(pid),
          n = n_pid,
          reason = "rescued_with_site_metadata_fallback",
          detail = "Used site_start/site_end because full chronology match was unavailable."
        )
      }
    }
  }

  still_bad_idx <- which(
    is.na(joined$system_name) |
      is.na(joined$phase_name) |
      is.na(joined$horizon_start) |
      is.na(joined$horizon_end) |
      is.na(joined$fade_in_start) |
      is.na(joined$fade_out_end)
  )

  if (length(still_bad_idx) > 0) {
    bad_ids <- unique(joined$phase_id[still_bad_idx])

    for (pid in bad_ids) {
      n_pid <- sum(joined$phase_id[still_bad_idx] == pid, na.rm = TRUE)
      log_tbl <- log_exclusions_add(
        log = log_tbl,
        stage = "resolve_allocation_bounds",
        entity_type = "phase",
        entity_id = as.character(pid),
        n = n_pid,
        reason = "missing_or_incomplete_chronology_match",
        detail = "Dropped because no complete chronology match was available and no usable site metadata fallback was found."
      )
    }

    warning(
      paste0(
        "Some mortuary rows could not be matched to a complete chronology definition and were dropped. ",
        "Affected phase_id values: ",
        paste(stats::na.omit(bad_ids), collapse = ", ")
      ),
      call. = FALSE
    )

    joined <- joined[-still_bad_idx, , drop = FALSE]
  }

  synthetic_tbl <- NULL

  if (augment_to_site_size) {
    if (is.null(site_metadata)) {
      warning(
        "Site-size augmentation was requested but no `site_metadata` was supplied. No synthetic burials were generated.",
        call. = FALSE
      )
    } else {
      meta_aug <- site_metadata[, c("site_id", "site_start", "site_end", "site_size"), drop = FALSE]
      meta_aug <- meta_aug[!duplicated(meta_aug$site_id), , drop = FALSE]

      obs_n <- as.data.frame(table(mortuary$site_id), stringsAsFactors = FALSE)
      names(obs_n) <- c("site_id", "observed_n")
      obs_n$observed_n <- as.integer(obs_n$observed_n)

      aug_df <- dplyr::left_join(meta_aug, obs_n, by = "site_id")
      aug_df$observed_n[is.na(aug_df$observed_n)] <- 0L
      aug_df$site_size_num <- suppressWarnings(as.numeric(aug_df$site_size))
      aug_df$site_start_num <- suppressWarnings(as.numeric(aug_df$site_start))
      aug_df$site_end_num <- suppressWarnings(as.numeric(aug_df$site_end))

      can_aug <- !is.na(aug_df$site_size_num) &
        !is.na(aug_df$site_start_num) &
        !is.na(aug_df$site_end_num) &
        aug_df$site_size_num > aug_df$observed_n

      aug_df$to_generate <- ifelse(
        can_aug,
        pmax(as.integer(round(aug_df$site_size_num)) - aug_df$observed_n, 0L),
        0L
      )

      synth_parts <- list()
      next_rec <- .safe_max_numeric(mortuary$record_id, default = 0)

      if (nrow(joined) > 0) {
        uid_n <- nrow(joined)
      } else {
        uid_n <- 0
      }

      synth_counter <- 0L

      for (i in seq_len(nrow(aug_df))) {
        n_gen <- aug_df$to_generate[i]
        if (!isTRUE(n_gen > 0)) next

        rec_ids <- seq.int(from = next_rec + 1, length.out = n_gen)
        next_rec <- max(rec_ids)

        uid_vals <- sprintf("%s%07d", synthetic_prefix, seq.int(from = synth_counter + 1L, length.out = n_gen))
        synth_counter <- synth_counter + n_gen

        synth_parts[[length(synth_parts) + 1L]] <- tibble::tibble(
          UID = uid_vals,
          input_type = rep("synthetic", n_gen),
          record_id = rec_ids,
          burial_id = rep(NA_character_, n_gen),
          unwrap_index = rep(NA_integer_, n_gen),
          site_id = rep(aug_df$site_id[i], n_gen),
          phase_id = rep(NA, n_gen),
          sex_gender = rep(NA_character_, n_gen),
          age = rep(NA_character_, n_gen),
          system_name = rep("site_metadata_fallback", n_gen),
          phase_name = rep("site_metadata_fallback", n_gen),
          horizon_start = rep(aug_df$site_start_num[i], n_gen),
          horizon_end = rep(aug_df$site_end_num[i], n_gen),
          fade_in_start = rep(aug_df$site_start_num[i], n_gen),
          fade_out_end = rep(aug_df$site_end_num[i], n_gen),
          chronology_source = rep("site_metadata_fallback", n_gen),
          is_synthetic = rep(TRUE, n_gen)
        )

        log_tbl <- log_exclusions_add(
          log = log_tbl,
          stage = "resolve_allocation_bounds",
          entity_type = "site",
          entity_id = as.character(aug_df$site_id[i]),
          n = n_gen,
          reason = "synthetic_burials_generated",
          detail = paste0(
            "Observed burials=", aug_df$observed_n[i],
            "; site_size=", aug_df$site_size_num[i],
            "; generated=", n_gen,
            ". Used site_start/site_end for fallback allocation."
          )
        )
      }

      if (length(synth_parts) > 0) {
        synthetic_tbl <- dplyr::bind_rows(synth_parts)
      }
    }
  }

  joined$core_profile <- ifelse(joined$chronology_source == "site_metadata_fallback", "uniform", "trapezoid")

  resolved_real <- joined[, c(
    "UID", "input_type", "record_id", "burial_id", "unwrap_index",
    "site_id", "phase_id", "sex_gender", "age",
    "system_name", "phase_name",
    "horizon_start", "horizon_end", "fade_in_start", "fade_out_end",
    "chronology_source", "is_synthetic", "core_profile"
  ), drop = FALSE]

  if (!is.null(synthetic_tbl)) {
    synthetic_tbl$core_profile <- "uniform"
    resolved <- dplyr::bind_rows(resolved_real, synthetic_tbl)
  } else {
    resolved <- resolved_real
  }

  list(
    data = tibble::as_tibble(resolved),
    logs = list(exclusions = log_tbl),
    diagnostics = list(
      n_input_rows = nrow(mortuary),
      n_resolved_real_rows = nrow(resolved_real),
      n_synthetic_rows = if (is.null(synthetic_tbl)) 0L else nrow(synthetic_tbl),
      n_output_rows = nrow(resolved),
      n_dropped_missing_chronology = sum(log_tbl$n[log_tbl$reason == "missing_or_incomplete_chronology_match"], na.rm = TRUE),
      n_rescued_site_metadata = sum(log_tbl$n[log_tbl$reason == "rescued_with_site_metadata_fallback"], na.rm = TRUE),
      n_generated_synthetic = sum(log_tbl$n[log_tbl$reason == "synthetic_burials_generated"], na.rm = TRUE)
    )
  )
}

# -------------------------------------------------------------------------
# Main harmonization / allocation
# -------------------------------------------------------------------------

#' Harmonize chronology and allocate burials across temporal buckets
#'
#' @description
#' Main temporal allocation engine for the merged workflow.
#'
#' Burials are first matched to phase-level chronology definitions via
#' \code{phase_id}. Each phase defines a trapezoidal temporal profile from
#' \code{fade_in_start}, \code{horizon_start}, \code{horizon_end}, and
#' \code{fade_out_end}. These phase-level profiles are then used to assign
#' burials to temporal buckets.
#'
#' The default mode is \code{"stochastic"}, meaning that each burial is allocated
#' to exactly one bucket per replicate by random draw from the normalized bucket
#' weights implied by its phase-level trapezoidal profile. This makes between-
#' replicate variation interpretable as temporal allocation uncertainty.
#'
#' If \code{use_site_metadata_fallback = TRUE}, burials lacking a full phase
#' chronology match may instead be allocated from \code{site_start/site_end}
#' in metadata, using a uniform site-level interval. If
#' \code{augment_to_site_size = TRUE}, synthetic burials may also be generated
#' up to \code{site_size} and allocated from the same site-level bounds.
#'
#' @param prepared_inputs Optional assembled object from \code{assemble_prepared_inputs()}.
#' @param mortuary Optional burial-level mortuary tibble.
#' @param chronology Optional chronology tibble.
#' @param site_metadata Optional site metadata tibble.
#' @param bucket_width Width of the temporal buckets.
#' @param offset Bucket grid offset.
#' @param allocation_mode One of \code{"stochastic"}, \code{"fractional"}, or
#'   \code{"deterministic_max"}. Default is \code{"stochastic"}.
#' @param profile Default temporal profile for phase-based allocation:
#'   \code{"trapezoid"} or \code{"uniform"}.
#' @param seed Optional seed for stochastic allocation.
#' @param use_site_metadata_fallback Whether unmatched burials may fall back to
#'   \code{site_start/site_end}.
#' @param augment_to_site_size Whether synthetic burials may be generated up to
#'   \code{site_size} using site-level bounds.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{data}: bucket-level allocation table
#'   \item \code{logs}: exclusion and fallback logs
#'   \item \code{diagnostics}: allocation diagnostics
#' }
harmonize_chronologies_merged <- function(prepared_inputs = NULL,
                                          mortuary = NULL,
                                          chronology = NULL,
                                          site_metadata = NULL,
                                          bucket_width = 25,
                                          offset = 0,
                                          allocation_mode = c("stochastic", "fractional", "deterministic_max"),
                                          profile = c("trapezoid", "uniform"),
                                          seed = NULL,
                                          use_site_metadata_fallback = FALSE,
                                          augment_to_site_size = FALSE) {
  allocation_mode <- match.arg(allocation_mode)
  profile <- match.arg(profile)
  .assert_positive_scalar(bucket_width, "bucket_width")

  if (!is.null(prepared_inputs)) {
    if (is.null(prepared_inputs$data$mortuary) || is.null(prepared_inputs$data$chronology)) {
      rlang::abort("`prepared_inputs` must contain `$data$mortuary` and `$data$chronology`.")
    }
    mortuary <- prepared_inputs$data$mortuary
    chronology <- prepared_inputs$data$chronology
    if (!is.null(prepared_inputs$data$site_metadata)) {
      site_metadata <- prepared_inputs$data$site_metadata
    }
  }

  if (is.null(mortuary) || is.null(chronology)) {
    rlang::abort("Supply either `prepared_inputs` or both `mortuary` and `chronology`.")
  }

  mortuary <- tibble::as_tibble(mortuary)
  chronology <- tibble::as_tibble(chronology)
  if (!is.null(site_metadata)) {
    site_metadata <- tibble::as_tibble(site_metadata)
  }

  resolved <- resolve_allocation_bounds(
    mortuary = mortuary,
    chronology = chronology,
    site_metadata = site_metadata,
    use_site_metadata_fallback = use_site_metadata_fallback,
    augment_to_site_size = augment_to_site_size
  )

  resolved_data <- resolved$data
  allocation_log <- resolved$logs$exclusions

  if (nrow(resolved_data) == 0) {
    rlang::abort(
      "No allocatable mortuary rows remain after chronology resolution and fallback handling."
    )
  }

  if (!is.null(seed) && allocation_mode == "stochastic") {
    set.seed(seed)
    row_seeds <- sample.int(.Machine$integer.max, nrow(resolved_data), replace = TRUE)
  } else {
    row_seeds <- rep(NA_integer_, nrow(resolved_data))
  }

  out_list <- lapply(seq_len(nrow(resolved_data)), function(i) {
    row_i <- resolved_data[i, , drop = FALSE]

    row_profile <- if (!is.null(row_i$core_profile) && !is.na(row_i$core_profile[[1]])) {
      row_i$core_profile[[1]]
    } else {
      profile
    }

    allocate_one_row(
      row_one = row_i,
      bucket_width = bucket_width,
      offset = offset,
      allocation_mode = allocation_mode,
      profile = row_profile,
      seed = row_seeds[i]
    )
  })

  out <- dplyr::bind_rows(out_list)

  diagnostics <- list(
    n_input_rows = nrow(mortuary),
    n_rows_after_resolution = nrow(resolved_data),
    n_dropped_missing_chronology = sum(allocation_log$n[allocation_log$reason == "missing_or_incomplete_chronology_match"], na.rm = TRUE),
    n_rescued_site_metadata = sum(allocation_log$n[allocation_log$reason == "rescued_with_site_metadata_fallback"], na.rm = TRUE),
    n_generated_synthetic = sum(allocation_log$n[allocation_log$reason == "synthetic_burials_generated"], na.rm = TRUE),
    n_output_rows = nrow(out),
    n_unique_UID = dplyr::n_distinct(out$UID),
    n_unique_sites = dplyr::n_distinct(out$site_id),
    n_unique_phases = dplyr::n_distinct(out$phase_id),
    allocation_mode = allocation_mode,
    bucket_width = bucket_width,
    offset = offset,
    mean_buckets_per_UID = if (nrow(out) == 0) 0 else mean(as.numeric(table(out$UID)))
  )

  list(
    data = out,
    logs = list(
      exclusions = allocation_log
    ),
    diagnostics = diagnostics
  )
}
