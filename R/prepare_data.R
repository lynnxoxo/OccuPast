# prepare_data.R
#
# Domain-level data preparation functions for OccuPast.
#
# Core design:
# - chronology is a strict phase-system table
# - chronology joins happen by phase_id
# - site_id is the analytical key for aggregation/bootstrap
# - site_name comes from site metadata, not mortuary inputs
# - individual mortuary input is already burial-level: one row = one burial
# - aggregated mortuary input is unwrapped to burial-level pseudo-records
# - burial_count exists only in prepared aggregated data, and is dropped after unwrapping
# - UID is assigned only after burial-level combination and deterministic sorting
# - sex_gender and age are preserved as pass-through attributes only

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   stringr
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows n_distinct arrange
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort

# -------------------------------------------------------------------------
# Validation helpers
# -------------------------------------------------------------------------

#' Validate required fields in a table
#'
#' @param data A data.frame or tibble.
#' @param required Character vector of required column names.
#' @param table_name Friendly name for error messages.
#'
#' @return Invisibly returns TRUE if validation passes.
validate_required_fields <- function(data, required, table_name = deparse(substitute(data))) {
  .require_data_frame(data, table_name)

  if (!is.character(required) || length(required) == 0) {
    rlang::abort("`required` must be a non-empty character vector.")
  }

  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    rlang::abort(
      paste0(
        "`", table_name, "` is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}

.empty_coercion_log <- function() {
  tibble::tibble(
    field = character(),
    row = integer(),
    original_value = character(),
    coerced_value = character(),
    issue = character()
  )
}

# -------------------------------------------------------------------------
# Chronology preparation
# -------------------------------------------------------------------------

#' Build a strict chronology table
#'
#' @param data A chronology source table.
#' @param colmap Named character vector mapping canonical names to source names.
#'   Allowed canonical names:
#'   - phase_id
#'   - system_name
#'   - phase_name
#'   - horizon_start
#'   - horizon_end
#'   - fade_in_start
#'   - fade_out_end
#' @param drop_invalid If TRUE, drop rows failing chronology-order checks and log them.
#'
#' @return A list with:
#'   - data: chronology tibble
#'   - logs: list(rename, coercions, exclusions)
#'   - diagnostics: list
build_chronology_table <- function(data, colmap, drop_invalid = TRUE) {
  .require_data_frame(data, "data")
  .require_named_character(colmap, "colmap")

  allowed_targets <- c(
    "phase_id",
    "system_name",
    "phase_name",
    "horizon_start",
    "horizon_end",
    "fade_in_start",
    "fade_out_end"
  )

  bad_targets <- setdiff(names(colmap), allowed_targets)
  if (length(bad_targets) > 0) {
    rlang::abort(
      paste0(
        "`colmap` contains unsupported canonical names: ",
        paste(bad_targets, collapse = ", "),
        ". Allowed names are: ",
        paste(allowed_targets, collapse = ", ")
      )
    )
  }

  exclusions <- log_exclusions_init()

  ren <- rename_columns_flex(data, colmap, strict = FALSE)
  dat <- ren$data

  dat <- .add_missing_columns(dat, allowed_targets, fill = NA)
  dat <- tibble::as_tibble(dat[, allowed_targets, drop = FALSE])

  for (nm in c("phase_id", "system_name", "phase_name")) {
    dat[[nm]] <- .null_if_empty_string(dat[[nm]])
    if (is.character(dat[[nm]])) {
      dat[[nm]] <- stringr::str_trim(dat[[nm]])
    }
  }

  # If no explicit phase_id is supplied, derive a stable join key from
  # system_name and phase_name. Chronological joins downstream depend on this.
  need_phase_id <- is.na(dat$phase_id) | dat$phase_id == ""
  if (any(need_phase_id)) {
    dat$phase_id[need_phase_id] <- ifelse(
      !is.na(dat$system_name[need_phase_id]) & !is.na(dat$phase_name[need_phase_id]),
      paste(dat$system_name[need_phase_id], dat$phase_name[need_phase_id], sep = "::"),
      NA_character_
    )
  }

  # Coerce chronology bounds with an audit trail; invalid dates are logged
  # rather than silently discarded.
  coe <- coerce_numeric_fields(
    dat,
    fields = c("horizon_start", "horizon_end", "fade_in_start", "fade_out_end"),
    na_strings = c("", "NA", "N/A", "unknown")
  )
  dat <- coe$data

  # Missing fade bounds collapse to the core horizon, giving a rectangular
  # profile instead of dropping otherwise usable phases.
  dat$fade_in_start[is.na(dat$fade_in_start) & !is.na(dat$horizon_start)] <-
    dat$horizon_start[is.na(dat$fade_in_start) & !is.na(dat$horizon_start)]

  dat$fade_out_end[is.na(dat$fade_out_end) & !is.na(dat$horizon_end)] <-
    dat$horizon_end[is.na(dat$fade_out_end) & !is.na(dat$horizon_end)]

  missing_core <- is.na(dat$phase_id) |
    is.na(dat$system_name) |
    is.na(dat$phase_name) |
    is.na(dat$horizon_start) |
    is.na(dat$horizon_end) |
    is.na(dat$fade_in_start) |
    is.na(dat$fade_out_end)

  bad_fade_in <- !is.na(dat$fade_in_start) &
    !is.na(dat$horizon_start) &
    (dat$fade_in_start > dat$horizon_start)

  bad_horizon <- !is.na(dat$horizon_start) &
    !is.na(dat$horizon_end) &
    (dat$horizon_start > dat$horizon_end)

  bad_fade_out <- !is.na(dat$horizon_end) &
    !is.na(dat$fade_out_end) &
    (dat$horizon_end > dat$fade_out_end)

  invalid_idx <- which(missing_core | bad_fade_in | bad_horizon | bad_fade_out)

  if (drop_invalid && length(invalid_idx) > 0) {
    for (i in invalid_idx) {
      reason_i <- if (missing_core[i]) {
        "missing_required_chronology_value"
      } else if (bad_fade_in[i]) {
        "fade_in_after_horizon_start"
      } else if (bad_horizon[i]) {
        "horizon_start_after_horizon_end"
      } else if (bad_fade_out[i]) {
        "horizon_end_after_fade_out_end"
      } else {
        "invalid_chronology_row"
      }

      exclusions <- log_exclusions_add(
        log = exclusions,
        stage = "build_chronology_table",
        entity_type = "phase",
        entity_id = dat$phase_id[i],
        n = 1L,
        reason = reason_i,
        detail = sprintf("row=%s", i)
      )
    }
    dat <- dat[-invalid_idx, , drop = FALSE]
  }

  dat <- tibble::as_tibble(dat[, allowed_targets, drop = FALSE])

  diagnostics <- list(
    n_input_rows = nrow(data),
    n_output_rows = nrow(dat),
    n_dropped = nrow(exclusions),
    rename_log = ren$rename_log,
    missing_mapped_columns = ren$missing_mapped,
    coercion_log = coe$log,
    invalid_rows_detected = length(invalid_idx),
    invalid_rows_dropped = if (drop_invalid) length(invalid_idx) else 0L,
    chronology_rules = c(
      "fade_in_start <= horizon_start",
      "horizon_start <= horizon_end",
      "horizon_end <= fade_out_end"
    )
  )

  list(
    data = dat,
    logs = list(
      rename = ren$rename_log,
      coercions = coe$log,
      exclusions = exclusions
    ),
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Mortuary preparation: individual
# -------------------------------------------------------------------------

#' Prepare individual mortuary records
#'
#' @param data Raw individual-level mortuary table.
#' @param colmap Named character vector mapping canonical names to source names.
#'   Supported canonical names:
#'   - record_id
#'   - burial_id
#'   - site_id
#'   - phase_id
#'   - sex_gender
#'   - age
#'
#' @return A list with:
#'   - data: canonical individual burial-level tibble
#'   - logs: list(rename, coercions, exclusions)
#'   - diagnostics: list
prepare_mortuary_individual <- function(data, colmap) {
  .require_data_frame(data, "data")
  .require_named_character(colmap, "colmap")

  allowed_targets <- c(
    "record_id",
    "burial_id",
    "site_id",
    "phase_id",
    "sex_gender",
    "age"
  )

  bad_targets <- setdiff(names(colmap), allowed_targets)
  if (length(bad_targets) > 0) {
    rlang::abort(
      paste0(
        "`colmap` contains unsupported canonical names: ",
        paste(bad_targets, collapse = ", "),
        ". Allowed names are: ",
        paste(allowed_targets, collapse = ", ")
      )
    )
  }

  exclusions <- log_exclusions_init()

  ren <- rename_columns_flex(data, colmap, strict = FALSE)
  dat <- ren$data
  dat <- .add_missing_columns(dat, allowed_targets, fill = NA)

  for (nm in c("burial_id", "site_id", "phase_id", "sex_gender", "age")) {
    dat[[nm]] <- .null_if_empty_string(dat[[nm]])
    if (is.character(dat[[nm]])) {
      dat[[nm]] <- stringr::str_trim(dat[[nm]])
    }
  }

  dat$record_id <- .null_if_empty_string(dat$record_id)
  if (is.character(dat$record_id)) {
    dat$record_id <- stringr::str_trim(dat$record_id)
  }

  coe <- coerce_numeric_fields(
    dat,
    fields = "record_id",
    na_strings = c("", "NA", "N/A", "unknown")
  )
  dat <- coe$data
  coercion_log <- coe$log

  # record_id is an ingest key only. Fill missing IDs deterministically so
  # combined burial-level tables can later receive stable UIDs.
  if (all(is.na(dat$record_id))) {
    dat$record_id <- seq_len(nrow(dat))
  } else {
    need_record_id <- is.na(dat$record_id)
    if (any(need_record_id)) {
      next_id <- suppressWarnings(max(dat$record_id, na.rm = TRUE))
      if (!is.finite(next_id)) next_id <- 0
      dat$record_id[need_record_id] <- seq.int(from = next_id + 1, length.out = sum(need_record_id))
    }
  }

  need_burial_id <- is.na(dat$burial_id) | dat$burial_id == ""
  if (any(need_burial_id)) {
    dat$burial_id[need_burial_id] <- paste0("burial_", dat$record_id[need_burial_id])
  }

  bad_record <- is.na(dat$record_id)
  bad_site <- is.na(dat$site_id)
  bad_phase <- is.na(dat$phase_id)

  invalid_idx <- which(bad_record | bad_site | bad_phase)

  if (length(invalid_idx) > 0) {
    for (i in invalid_idx) {
      reason_i <- if (bad_record[i]) {
        "missing_record_id"
      } else if (bad_site[i]) {
        "missing_site_id"
      } else {
        "missing_phase_id"
      }

      exclusions <- log_exclusions_add(
        log = exclusions,
        stage = "prepare_mortuary_individual",
        entity_type = "burial",
        entity_id = dat$burial_id[i],
        n = 1L,
        reason = reason_i,
        detail = sprintf("row=%s", i)
      )
    }
    dat <- dat[-invalid_idx, , drop = FALSE]
  }

  out <- tibble::as_tibble(
    dat[, c(
      "record_id",
      "burial_id",
      "site_id",
      "phase_id",
      "sex_gender",
      "age"
    ), drop = FALSE]
  )

  diagnostics <- list(
    n_input_rows = nrow(data),
    n_output_rows = nrow(out),
    n_dropped = nrow(exclusions),
    n_unique_sites = dplyr::n_distinct(out$site_id),
    n_unique_phases = dplyr::n_distinct(out$phase_id),
    n_unique_records = dplyr::n_distinct(out$record_id),
    n_unique_burials = dplyr::n_distinct(out$burial_id)
  )

  list(
    data = out,
    logs = list(
      rename = ren$rename_log,
      coercions = coercion_log,
      exclusions = exclusions
    ),
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Mortuary preparation: aggregated
# -------------------------------------------------------------------------

#' Prepare aggregated mortuary records
#'
#' @param data Raw aggregated mortuary table.
#' @param colmap Named character vector mapping canonical names to source names.
#'   Supported canonical names:
#'   - record_id
#'   - site_id
#'   - phase_id
#'   - burial_count
#'   - sex_gender
#'   - age
#' @param generate_record_id If TRUE, generate numeric record IDs when absent.
#'
#' @return A list with:
#'   - data: canonical aggregated tibble
#'   - logs: list(rename, coercions, exclusions)
#'   - diagnostics: list
prepare_mortuary_aggregated <- function(data, colmap,
                                        generate_record_id = TRUE) {
  .require_data_frame(data, "data")
  .require_named_character(colmap, "colmap")

  allowed_targets <- c(
    "record_id",
    "site_id",
    "phase_id",
    "burial_count",
    "sex_gender",
    "age"
  )

  bad_targets <- setdiff(names(colmap), allowed_targets)
  if (length(bad_targets) > 0) {
    rlang::abort(
      paste0(
        "`colmap` contains unsupported canonical names: ",
        paste(bad_targets, collapse = ", "),
        ". Allowed names are: ",
        paste(allowed_targets, collapse = ", ")
      )
    )
  }

  exclusions <- log_exclusions_init()

  ren <- rename_columns_flex(data, colmap, strict = FALSE)
  dat <- ren$data
  dat <- .add_missing_columns(dat, allowed_targets, fill = NA)

  for (nm in c("site_id", "phase_id", "sex_gender", "age")) {
    dat[[nm]] <- .null_if_empty_string(dat[[nm]])
    if (is.character(dat[[nm]])) {
      dat[[nm]] <- stringr::str_trim(dat[[nm]])
    }
  }

  dat$record_id <- .null_if_empty_string(dat$record_id)
  if (is.character(dat$record_id)) {
    dat$record_id <- stringr::str_trim(dat$record_id)
  }

  if (generate_record_id) {
    # Aggregated rows may not have record IDs; deterministic IDs make later
    # unwrapping reproducible.
    need_id <- is.na(dat$record_id) | dat$record_id == ""
    if (any(need_id)) {
      dat$record_id[need_id] <- seq_len(sum(need_id))
    }
  }

  coe <- coerce_numeric_fields(
    dat,
    fields = c("record_id", "burial_count"),
    na_strings = c("", "NA", "N/A", "unknown")
  )
  dat <- coe$data

  bad_id <- is.na(dat$record_id)
  bad_site <- is.na(dat$site_id)
  bad_phase <- is.na(dat$phase_id)
  bad_count <- is.na(dat$burial_count) | dat$burial_count <= 0

  invalid_idx <- which(bad_id | bad_site | bad_phase | bad_count)

  if (length(invalid_idx) > 0) {
    for (i in invalid_idx) {
      reason_i <- if (bad_id[i]) {
        "missing_record_id"
      } else if (bad_site[i]) {
        "missing_site_id"
      } else if (bad_phase[i]) {
        "missing_phase_id"
      } else {
        "invalid_burial_count"
      }

      exclusions <- log_exclusions_add(
        log = exclusions,
        stage = "prepare_mortuary_aggregated",
        entity_type = "record",
        entity_id = as.character(dat$record_id[i]),
        n = 1L,
        reason = reason_i,
        detail = sprintf("row=%s", i)
      )
    }
    dat <- dat[-invalid_idx, , drop = FALSE]
  }

  out <- tibble::as_tibble(
    dat[, c(
      "record_id",
      "site_id",
      "phase_id",
      "burial_count",
      "sex_gender",
      "age"
    ), drop = FALSE]
  )

  diagnostics <- list(
    n_input_rows = nrow(data),
    n_output_rows = nrow(out),
    n_dropped = nrow(exclusions),
    n_unique_sites = dplyr::n_distinct(out$site_id),
    n_unique_phases = dplyr::n_distinct(out$phase_id),
    n_unique_records = dplyr::n_distinct(out$record_id)
  )

  list(
    data = out,
    logs = list(
      rename = ren$rename_log,
      coercions = coe$log,
      exclusions = exclusions
    ),
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Site metadata preparation
# -------------------------------------------------------------------------

#' Prepare site metadata for visualization and lookup
#'
#' @param data Raw site metadata table.
#' @param colmap Named character vector mapping canonical names to source names.
#'   Supported canonical names:
#'   - site_id
#'   - site_name
#'   - coord_y
#'   - coord_x
#'   - site_country
#'   - site_region
#'   - site_admin
#'   - site_start
#'   - site_end
#'   - site_size
#'   - site_dig_date
#'
#' @return A list with:
#'   - data: canonical site metadata tibble
#'   - logs: list(rename, coercions, exclusions)
#'   - diagnostics: list
prepare_site_metadata <- function(data, colmap) {
  .require_data_frame(data, "data")
  .require_named_character(colmap, "colmap")

  allowed_targets <- c(
    "site_id",
    "site_name",
    "coord_y",
    "coord_x",
    "site_country",
    "site_region",
    "site_admin",
    "site_start",
    "site_end",
    "site_size",
    "site_dig_date"
  )

  bad_targets <- setdiff(names(colmap), allowed_targets)
  if (length(bad_targets) > 0) {
    rlang::abort(
      paste0(
        "`colmap` contains unsupported canonical names: ",
        paste(bad_targets, collapse = ", "),
        ". Allowed names are: ",
        paste(allowed_targets, collapse = ", ")
      )
    )
  }

  exclusions <- log_exclusions_init()

  ren <- rename_columns_flex(data, colmap, strict = FALSE)
  dat <- ren$data
  dat <- .add_missing_columns(dat, allowed_targets, fill = NA)

  text_fields <- c(
    "site_id",
    "site_name",
    "site_country",
    "site_region",
    "site_admin",
    "site_start",
    "site_end",
    "site_dig_date"
  )

  for (nm in text_fields) {
    dat[[nm]] <- .null_if_empty_string(dat[[nm]])
    if (is.character(dat[[nm]])) {
      dat[[nm]] <- stringr::str_trim(dat[[nm]])
    }
  }

  coe <- coerce_numeric_fields(
    dat,
    fields = c("coord_y", "coord_x", "site_size"),
    na_strings = c("", "NA", "N/A", "unknown")
  )
  dat <- coe$data

  bad_site <- is.na(dat$site_id)
  invalid_idx <- which(bad_site)

  if (length(invalid_idx) > 0) {
    for (i in invalid_idx) {
      exclusions <- log_exclusions_add(
        log = exclusions,
        stage = "prepare_site_metadata",
        entity_type = "site",
        entity_id = NA_character_,
        n = 1L,
        reason = "missing_site_id",
        detail = sprintf("row=%s", i)
      )
    }
    dat <- dat[-invalid_idx, , drop = FALSE]
  }

  dup_sites <- duplicated(dat$site_id)
  dup_idx <- which(dup_sites)

  if (length(dup_idx) > 0) {
    for (i in dup_idx) {
      exclusions <- log_exclusions_add(
        log = exclusions,
        stage = "prepare_site_metadata",
        entity_type = "site",
        entity_id = dat$site_id[i],
        n = 1L,
        reason = "duplicate_site_id_dropped",
        detail = sprintf("row=%s", i)
      )
    }
    dat <- dat[!dup_sites, , drop = FALSE]
  }

  out <- tibble::as_tibble(
    dat[, c(
      "site_id",
      "site_name",
      "coord_y",
      "coord_x",
      "site_country",
      "site_region",
      "site_admin",
      "site_start",
      "site_end",
      "site_size",
      "site_dig_date"
    ), drop = FALSE]
  )

  diagnostics <- list(
    n_input_rows = nrow(data),
    n_output_rows = nrow(out),
    n_dropped = nrow(exclusions),
    n_unique_sites = dplyr::n_distinct(out$site_id),
    n_with_coordinates = sum(!is.na(out$coord_y) & !is.na(out$coord_x)),
    n_with_site_dates = sum(!is.na(out$site_start) & !is.na(out$site_end)),
    n_with_site_size = sum(!is.na(out$site_size))
  )

  list(
    data = out,
    logs = list(
      rename = ren$rename_log,
      coercions = coe$log,
      exclusions = exclusions
    ),
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Phase-link validation
# -------------------------------------------------------------------------

#' Validate phase links between mortuary and chronology tables
#'
#' @param mortuary A burial-level mortuary table containing phase_id.
#' @param chronology A chronology table containing phase_id.
#'
#' @return A list with join diagnostics.
validate_phase_links <- function(mortuary, chronology) {
  validate_required_fields(mortuary, "phase_id", "mortuary")
  validate_required_fields(chronology, "phase_id", "chronology")

  mort_phases <- sort(unique(stats::na.omit(mortuary$phase_id)))
  chr_phases <- chronology$phase_id

  duplicate_chronology_phase_ids <- sort(unique(chr_phases[duplicated(chr_phases)]))
  unmatched_phase_ids <- sort(setdiff(mort_phases, unique(chr_phases)))
  matched_phase_ids <- sort(intersect(mort_phases, unique(chr_phases)))

  list(
    n_mortuary_rows = nrow(mortuary),
    n_distinct_mortuary_phases = length(mort_phases),
    n_distinct_chronology_phases = dplyr::n_distinct(chr_phases),
    n_matched_phase_ids = length(matched_phase_ids),
    n_unmatched_phase_ids = length(unmatched_phase_ids),
    unmatched_phase_ids = unmatched_phase_ids,
    duplicate_chronology_phase_ids = duplicate_chronology_phase_ids,
    ok = length(unmatched_phase_ids) == 0 && length(duplicate_chronology_phase_ids) == 0
  )
}

# -------------------------------------------------------------------------
# Aggregate unwrapping
# -------------------------------------------------------------------------

#' Unwrap aggregated mortuary rows to burial-level pseudo-records
#'
#' @param aggregated Prepared output from prepare_mortuary_aggregated(),
#'   or a canonical aggregated data.frame/tibble.
#'
#' @return A list with:
#'   - data: burial-level pseudo-record tibble
#'   - diagnostics: list
unwrap_mortuary_aggregated <- function(aggregated) {
  agg_data <- if (is.list(aggregated) && !is.null(aggregated$data)) aggregated$data else aggregated

  validate_required_fields(
    agg_data,
    c("record_id", "site_id", "phase_id", "burial_count", "sex_gender", "age"),
    "aggregated"
  )

  agg_data <- tibble::as_tibble(agg_data)

  if (nrow(agg_data) == 0) {
    out <- tibble::tibble(
      record_id = numeric(),
      burial_id = character(),
      unwrap_index = integer(),
      site_id = character(),
      phase_id = character(),
      sex_gender = character(),
      age = character(),
      input_type = character()
    )

    return(list(
      data = out,
      diagnostics = list(
        n_input_rows = 0L,
        n_output_rows = 0L,
        total_unwrapped = 0L
      )
    ))
  }

  # Unwrapping converts count rows into pseudo-burial rows. The unwrap_index
  # preserves which pseudo-records came from the same aggregate input row.
  expanded_list <- vector("list", nrow(agg_data))

  for (i in seq_len(nrow(agg_data))) {
    n_i <- as.integer(agg_data$burial_count[i])

    expanded_list[[i]] <- tibble::tibble(
      record_id = rep(agg_data$record_id[i], n_i),
      burial_id = rep(NA_character_, n_i),
      unwrap_index = seq_len(n_i),
      site_id = rep(agg_data$site_id[i], n_i),
      phase_id = rep(agg_data$phase_id[i], n_i),
      sex_gender = rep(agg_data$sex_gender[i], n_i),
      age = rep(agg_data$age[i], n_i),
      input_type = rep("aggregated", n_i)
    )
  }

  out <- dplyr::bind_rows(expanded_list)

  diagnostics <- list(
    n_input_rows = nrow(agg_data),
    n_output_rows = nrow(out),
    total_unwrapped = nrow(out),
    min_unwrap_count = min(agg_data$burial_count),
    max_unwrap_count = max(agg_data$burial_count)
  )

  list(
    data = out,
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Burial-level mortuary combination
# -------------------------------------------------------------------------

#' Combine burial-level individual rows and unwrapped aggregate rows
#'
#' @param individual Prepared output from prepare_mortuary_individual(),
#'   or a canonical individual data.frame/tibble.
#' @param aggregated_unwrapped Output from unwrap_mortuary_aggregated(),
#'   or a canonical burial-level aggregated data.frame/tibble.
#'
#' @return A list with:
#'   - data: combined burial-level mortuary tibble (without UID)
#'   - diagnostics: list
combine_mortuary_burial_level <- function(individual = NULL, aggregated_unwrapped = NULL) {
  if (is.null(individual) && is.null(aggregated_unwrapped)) {
    rlang::abort("At least one of `individual` or `aggregated_unwrapped` must be supplied.")
  }

  parts <- list()

  if (!is.null(individual)) {
    ind_data <- if (is.list(individual) && !is.null(individual$data)) individual$data else individual

    validate_required_fields(
      ind_data,
      c("record_id", "burial_id", "site_id", "phase_id", "sex_gender", "age"),
      "individual"
    )

    ind_tbl <- tibble::as_tibble(ind_data)
    ind_tbl$unwrap_index <- NA_integer_
    ind_tbl$input_type <- "individual"

    ind_tbl <- ind_tbl[, c(
      "record_id",
      "burial_id",
      "unwrap_index",
      "site_id",
      "phase_id",
      "sex_gender",
      "age",
      "input_type"
    ), drop = FALSE]

    parts <- c(parts, list(ind_tbl))
  }

  if (!is.null(aggregated_unwrapped)) {
    agg_data <- if (is.list(aggregated_unwrapped) && !is.null(aggregated_unwrapped$data)) {
      aggregated_unwrapped$data
    } else {
      aggregated_unwrapped
    }

    validate_required_fields(
      agg_data,
      c("record_id", "burial_id", "unwrap_index", "site_id", "phase_id", "sex_gender", "age", "input_type"),
      "aggregated_unwrapped"
    )

    agg_tbl <- tibble::as_tibble(agg_data)
    agg_tbl <- agg_tbl[, c(
      "record_id",
      "burial_id",
      "unwrap_index",
      "site_id",
      "phase_id",
      "sex_gender",
      "age",
      "input_type"
    ), drop = FALSE]

    parts <- c(parts, list(agg_tbl))
  }

  combined <- dplyr::bind_rows(parts)

  ind_sites <- if (is.null(individual)) character() else {
    ind_data <- if (is.list(individual) && !is.null(individual$data)) individual$data else individual
    sort(unique(stats::na.omit(ind_data$site_id)))
  }

  agg_sites <- if (is.null(aggregated_unwrapped)) character() else {
    agg_data <- if (is.list(aggregated_unwrapped) && !is.null(aggregated_unwrapped$data)) aggregated_unwrapped$data else aggregated_unwrapped
    sort(unique(stats::na.omit(agg_data$site_id)))
  }

  ind_phases <- if (is.null(individual)) character() else {
    ind_data <- if (is.list(individual) && !is.null(individual$data)) individual$data else individual
    sort(unique(stats::na.omit(ind_data$phase_id)))
  }

  agg_phases <- if (is.null(aggregated_unwrapped)) character() else {
    agg_data <- if (is.list(aggregated_unwrapped) && !is.null(aggregated_unwrapped$data)) aggregated_unwrapped$data else aggregated_unwrapped
    sort(unique(stats::na.omit(agg_data$phase_id)))
  }

  diagnostics <- list(
    n_combined_rows = nrow(combined),
    n_individual_rows = if (is.null(individual)) 0L else nrow(if (is.list(individual) && !is.null(individual$data)) individual$data else individual),
    n_aggregated_unwrapped_rows = if (is.null(aggregated_unwrapped)) 0L else nrow(if (is.list(aggregated_unwrapped) && !is.null(aggregated_unwrapped$data)) aggregated_unwrapped$data else aggregated_unwrapped),
    n_unique_sites = dplyr::n_distinct(combined$site_id),
    n_unique_phases = dplyr::n_distinct(combined$phase_id),
    overlapping_sites = intersect(ind_sites, agg_sites),
    n_overlapping_sites = length(intersect(ind_sites, agg_sites)),
    overlapping_phases = intersect(ind_phases, agg_phases),
    n_overlapping_phases = length(intersect(ind_phases, agg_phases))
  )

  list(
    data = combined,
    diagnostics = diagnostics
  )
}

# -------------------------------------------------------------------------
# Final assembly
# -------------------------------------------------------------------------

#' Assemble prepared chronology, burial-level mortuary, and optional site metadata
#'
#' @param chronology Prepared output from build_chronology_table(), or chronology tibble.
#' @param mortuary_burial_level Output from combine_mortuary_burial_level(),
#'   or a combined burial-level mortuary tibble.
#' @param site_metadata Optional prepared output from prepare_site_metadata(),
#'   or site metadata tibble.
#'
#' @return A structured list with assembled data, logs, diagnostics, and settings.
assemble_prepared_inputs <- function(chronology,
                                     mortuary_burial_level,
                                     site_metadata = NULL) {
  chr_data <- if (is.list(chronology) && !is.null(chronology$data)) chronology$data else chronology
  mort_data <- if (is.list(mortuary_burial_level) && !is.null(mortuary_burial_level$data)) mortuary_burial_level$data else mortuary_burial_level
  meta_data <- NULL
  if (!is.null(site_metadata)) {
    meta_data <- if (is.list(site_metadata) && !is.null(site_metadata$data)) site_metadata$data else site_metadata
  }

  validate_required_fields(
    chr_data,
    c("phase_id", "system_name", "phase_name", "horizon_start", "horizon_end", "fade_in_start", "fade_out_end"),
    "chronology"
  )

  validate_required_fields(
    mort_data,
    c("record_id", "burial_id", "unwrap_index", "site_id", "phase_id", "sex_gender", "age", "input_type"),
    "mortuary_burial_level"
  )

  mort_tbl <- tibble::as_tibble(mort_data)

  # Deterministic ordering before UID assignment. UIDs are analysis keys, so
  # they must not depend on the original row order of input files. For
  # individual rows unwrap_index is NA; replace with 0 for stable sorting.
  sort_unwrap <- ifelse(is.na(mort_tbl$unwrap_index), 0L, as.integer(mort_tbl$unwrap_index))

  ord <- order(
    mort_tbl$input_type,
    mort_tbl$record_id,
    sort_unwrap,
    na.last = TRUE
  )

  mort_tbl <- mort_tbl[ord, , drop = FALSE]
  mort_tbl$UID <- sprintf("UID%07d", seq_len(nrow(mort_tbl)))

  mort_tbl <- mort_tbl[, c(
    "UID",
    "input_type",
    "record_id",
    "burial_id",
    "unwrap_index",
    "site_id",
    "phase_id",
    "sex_gender",
    "age"
  ), drop = FALSE]

  phase_diag <- validate_phase_links(mort_tbl, chr_data)

  if (!phase_diag$ok) {
    msgs <- character()

    if (length(phase_diag$duplicate_chronology_phase_ids) > 0) {
      msgs <- c(
        msgs,
        paste0(
          "Duplicate chronology `phase_id` values: ",
          paste(phase_diag$duplicate_chronology_phase_ids, collapse = ", ")
        )
      )
    }

    if (length(phase_diag$unmatched_phase_ids) > 0) {
      show_n <- min(10L, length(phase_diag$unmatched_phase_ids))
      msgs <- c(
        msgs,
        paste0(
          "Mortuary `phase_id` values not found in chronology (showing first ",
          show_n, "): ",
          paste(utils::head(phase_diag$unmatched_phase_ids, show_n), collapse = ", ")
        )
      )
    }

    rlang::abort(paste(msgs, collapse = " | "))
  }

  metadata_join_diag <- NULL
  if (!is.null(meta_data)) {
    validate_required_fields(meta_data, "site_id", "site_metadata")

    mort_sites <- sort(unique(stats::na.omit(mort_tbl$site_id)))
    meta_sites <- sort(unique(stats::na.omit(meta_data$site_id)))

    metadata_join_diag <- list(
      n_mortuary_sites = length(mort_sites),
      n_metadata_sites = length(meta_sites),
      n_matched_sites = length(intersect(mort_sites, meta_sites)),
      n_unmatched_mortuary_sites = length(setdiff(mort_sites, meta_sites)),
      unmatched_mortuary_sites = sort(setdiff(mort_sites, meta_sites))
    )
  }

  list(
    data = list(
      chronology = tibble::as_tibble(chr_data),
      mortuary = tibble::as_tibble(mort_tbl),
      site_metadata = if (is.null(meta_data)) NULL else tibble::as_tibble(meta_data)
    ),
    diagnostics = list(
      chronology = if (is.list(chronology) && !is.null(chronology$diagnostics)) chronology$diagnostics else NULL,
      mortuary_burial_level = if (is.list(mortuary_burial_level) && !is.null(mortuary_burial_level$diagnostics)) mortuary_burial_level$diagnostics else NULL,
      site_metadata = if (!is.null(site_metadata) && is.list(site_metadata) && !is.null(site_metadata$diagnostics)) site_metadata$diagnostics else NULL,
      phase_links = phase_diag,
      metadata_links = metadata_join_diag
    ),
    logs = list(
      chronology = if (is.list(chronology) && !is.null(chronology$logs)) chronology$logs else NULL,
      site_metadata = if (!is.null(site_metadata) && is.list(site_metadata) && !is.null(site_metadata$logs)) site_metadata$logs else NULL
    ),
    settings = list(
      uid_sort = c("input_type", "record_id", "unwrap_index")
    )
  )
}
