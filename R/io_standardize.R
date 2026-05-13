# io_standardize.R
#
# Shared ingest / standardization helpers for OccuPast.
#
#
# Public helpers exposed for reuse:
#   - rename_columns_flex()
#   - coerce_numeric_fields()
#   - log_exclusions_init()
#   - log_exclusions_add()
#
# Internal helpers used across preparation modules:
#   - .require_data_frame()
#   - .require_named_character()
#   - .as_character_na_safe()
#   - .null_if_empty_string()
#   - .add_missing_columns()
#   - .standardize_bin_columns()

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   stringr
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows rename
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort
#   @importFrom stringr str_trim

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.require_data_frame <- function(data, arg = "data") {
  if (!is.data.frame(data)) {
    rlang::abort(sprintf("`%s` must be a data.frame or tibble.", arg))
  }
}

.require_named_character <- function(x, arg = "colmap") {
  if (!is.character(x) || is.null(names(x)) || any(names(x) == "")) {
    rlang::abort(sprintf("`%s` must be a named character vector.", arg))
  }
}

.as_character_na_safe <- function(x) {
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  out
}

.assert_unique_names <- function(nms, arg = "data") {
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0) {
    rlang::abort(
      paste0(
        "`", arg, "` has duplicated column names: ",
        paste(dup, collapse = ", ")
      )
    )
  }
}

.null_if_empty_string <- function(x) {
  x <- .as_character_na_safe(x)
  x[stringr::str_trim(x) == ""] <- NA_character_
  x
}

.add_missing_columns <- function(data, cols, fill = NA) {
  for (nm in cols) {
    if (!nm %in% names(data)) {
      data[[nm]] <- fill
    }
  }
  data
}

# Standardize temporal-bin column names while keeping legacy inputs readable.
# Bins used to be called buckets. Now, they aren't anymore.
.standardize_bin_columns <- function(data, arg = "data") {
  .require_data_frame(data, arg)
  data <- tibble::as_tibble(data)

  legacy_map <- c(
    horizon_bucket = "horizon_bin",
    bucket_start = "bin_start",
    bucket_end = "bin_end",
    chosen_bucket = "chosen_bin"
  )

  for (old in names(legacy_map)) {
    new <- unname(legacy_map[[old]])
    has_old <- old %in% names(data)
    has_new <- new %in% names(data)

    if (has_old && has_new) {
      rlang::abort(
        paste0(
          "`", arg, "` contains both legacy column `", old,
          "` and canonical column `", new, "`. Please keep only one."
        )
      )
    }

    if (has_old && !has_new) {
      names(data)[names(data) == old] <- new
    }
  }

  data
}

.standard_exclusion_log <- function() {
  tibble::tibble(
    stage = character(),
    entity_type = character(),
    entity_id = character(),
    n = integer(),
    reason = character(),
    detail = character()
  )
}

.standard_coercion_log <- function() {
  tibble::tibble(
    field = character(),
    row = integer(),
    original_value = character(),
    coerced_value = character(),
    issue = character()
  )
}

# -------------------------------------------------------------------------
# Exclusion log helpers
# The exclusion log documents dropped entries for troubleshooting
# -------------------------------------------------------------------------

#' Initialize an exclusion log
#'
#' @return Empty tibble with standard exclusion-log columns.
log_exclusions_init <- function() {
  .standard_exclusion_log()
}

#' Append an exclusion event to a log
#'
#' @param log Existing exclusion log tibble.
#' @param stage Processing stage name.
#' @param entity_type Type of excluded entity, e.g. "phase", "burial", "site".
#' @param entity_id Optional entity identifier.
#' @param n Integer count represented by the exclusion event.
#' @param reason Short exclusion reason.
#' @param detail Optional longer explanation.
#'
#' @return Updated exclusion log tibble.
log_exclusions_add <- function(log, stage, entity_type, entity_id = NA_character_,
                               n = 1L, reason, detail = NULL) {
  if (missing(log) || is.null(log)) {
    log <- .standard_exclusion_log()
  }

  required_cols <- c("stage", "entity_type", "entity_id", "n", "reason", "detail")
  missing_cols <- setdiff(required_cols, names(log))
  if (length(missing_cols) > 0) {
    rlang::abort(
      paste0(
        "`log` is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (!is.character(stage) || length(stage) != 1L || is.na(stage) || stage == "") {
    rlang::abort("`stage` must be a single non-empty character string.")
  }
  if (!is.character(entity_type) || length(entity_type) != 1L || is.na(entity_type) || entity_type == "") {
    rlang::abort("`entity_type` must be a single non-empty character string.")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0) {
    rlang::abort("`n` must be a single non-negative number.")
  }
  if (!is.character(reason) || length(reason) != 1L || is.na(reason) || reason == "") {
    rlang::abort("`reason` must be a single non-empty character string.")
  }

  dplyr::bind_rows(
    tibble::as_tibble(log),
    tibble::tibble(
      stage = stage,
      entity_type = entity_type,
      entity_id = if (length(entity_id) == 0) NA_character_ else as.character(entity_id[[1]]),
      n = as.integer(n),
      reason = reason,
      detail = if (is.null(detail)) NA_character_ else as.character(detail[[1]])
    )
  )
}

# -------------------------------------------------------------------------
# Flexible renaming
# -------------------------------------------------------------------------

#' Rename columns flexibly using a canonical-to-source mapping
#'
#' @param data A data.frame or tibble.
#' @param colmap Named character vector. Names are canonical output names;
#'   values are source column names in `data`.
#' @param strict If TRUE, abort when a mapped source column is missing.
#'
#' @return A list with:
#'   - data: renamed tibble
#'   - rename_log: tibble describing rename actions
#'   - missing_mapped: character vector of missing mapped source columns
rename_columns_flex <- function(data, colmap, strict = FALSE) {
  .require_data_frame(data, "data")
  .require_named_character(colmap, "colmap")
  .assert_unique_names(names(data), "data")

  data <- tibble::as_tibble(data)
  original_names <- names(data)

  target_names <- names(colmap)
  source_names <- unname(colmap)

  if (any(is.na(target_names)) || any(target_names == "")) {
    rlang::abort("`colmap` names must all be non-empty canonical target names.")
  }
  if (any(is.na(source_names)) || any(source_names == "")) {
    rlang::abort("`colmap` values must all be non-empty source column names.")
  }
  if (anyDuplicated(target_names) > 0) {
    dup <- unique(target_names[duplicated(target_names)])
    rlang::abort(
      paste0(
        "Canonical target names are duplicated in `colmap`: ",
        paste(dup, collapse = ", ")
      )
    )
  }

  present <- source_names %in% original_names
  missing_sources <- source_names[!present]

  if (strict && length(missing_sources) > 0) {
    rlang::abort(
      paste0(
        "Mapped source columns are missing from `data`: ",
        paste(missing_sources, collapse = ", ")
      )
    )
  }

  rename_log <- tibble::tibble(
    target_name = target_names,
    source_name = source_names,
    source_present = present,
    action = ifelse(
      !present,
      "missing_source",
      ifelse(target_names == source_names, "kept", "renamed")
    )
  )

  if (any(present)) {
    rename_pairs <- stats::setNames(source_names[present], target_names[present])

    collisions <- names(rename_pairs)[
      names(rename_pairs) %in% names(data) &
        names(rename_pairs) != unname(rename_pairs)
    ]

    if (length(collisions) > 0) {
      rlang::abort(
        paste0(
          "Renaming would overwrite existing columns: ",
          paste(collisions, collapse = ", "),
          ". Rename or remove these columns first."
        )
      )
    }

    data <- dplyr::rename(data, !!!rename_pairs)
  }

  .assert_unique_names(names(data), "renamed data")

  list(
    data = data,
    rename_log = rename_log,
    missing_mapped = missing_sources
  )
}

# -------------------------------------------------------------------------
# Numeric coercion with audit trail
# -------------------------------------------------------------------------

#' Coerce selected fields to numeric with explicit logging
#'
#' @param data A data.frame or tibble.
#' @param fields Character vector of column names to coerce.
#' @param na_strings Character vector of strings to treat as NA before coercion.
#'
#' @return A list with:
#'   - data: tibble with coerced columns
#'   - log: tibble of coercion issues
#'
#' @details
#' For character/factor fields:
#' - values in `na_strings` become NA
#' - commas are converted to dots
#' - surrounding whitespace is trimmed
#' - non-convertible values are turned into NA and logged
coerce_numeric_fields <- function(data, fields, na_strings = c("", "NA", "N/A")) {
  .require_data_frame(data, "data")
  data <- tibble::as_tibble(data)

  if (!is.character(fields) || length(fields) == 0) {
    rlang::abort("`fields` must be a non-empty character vector.")
  }

  missing_fields <- setdiff(fields, names(data))
  if (length(missing_fields) > 0) {
    rlang::abort(
      paste0(
        "These fields are not present in `data`: ",
        paste(missing_fields, collapse = ", ")
      )
    )
  }

  log_tbl <- .standard_coercion_log()

  for (field in fields) {
    x <- data[[field]]
    original <- x

    if (is.numeric(x)) {
      next
    }

    chr <- .as_character_na_safe(x)
    chr <- stringr::str_trim(chr)
    chr[chr %in% na_strings] <- NA_character_
    chr <- gsub(",", ".", chr, fixed = TRUE)

    suppressWarnings(num <- as.numeric(chr))

    bad <- which(!is.na(chr) & is.na(num))
    if (length(bad) > 0) {
      log_tbl <- dplyr::bind_rows(
        log_tbl,
        tibble::tibble(
          field = field,
          row = bad,
          original_value = chr[bad],
          coerced_value = NA_character_,
          issue = "non_numeric_to_na"
        )
      )
    }

    changed_ok <- which(!is.na(chr) & !is.na(num))
    if (length(changed_ok) > 0 && !is.numeric(original)) {
      original_chr <- .as_character_na_safe(original)
      original_chr <- stringr::str_trim(original_chr)
      changed_repr <- changed_ok[as.character(num[changed_ok]) != original_chr[changed_ok]]

      if (length(changed_repr) > 0) {
        log_tbl <- dplyr::bind_rows(
          log_tbl,
          tibble::tibble(
            field = field,
            row = changed_repr,
            original_value = original_chr[changed_repr],
            coerced_value = as.character(num[changed_repr]),
            issue = "string_to_numeric"
          )
        )
      }
    }

    data[[field]] <- num
  }

  list(
    data = data,
    log = log_tbl
  )
}
