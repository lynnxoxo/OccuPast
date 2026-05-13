# export.R
#
# Export helpers for OccuPast.
#
# Design:
# - flexible: accepts any mix of prepared / allocation / analysis / ensemble / final / significance objects
# - CSV-first for portability
# - optional plot export through plotting.R helpers
# - manifest records what was written and what was missing
#
# Expected object types:
# - prepared: output from assemble_prepared_inputs()
# - alloc: output from harmonize_chronologies_merged()
# - analysis: output from analyze_occupancy_curve()
# - ens: output from run_temporal_ensemble()
# - final: output from finalize_ensemble()
# - significance: one result object or a named list of result objects from significance.R

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows
#   @importFrom tibble tibble as_tibble
#   @importFrom rlang abort

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

.dir_create_safe <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.write_csv_safe <- function(tbl, path, na = "") {
  if (is.null(tbl)) return(FALSE)
  tbl <- tibble::as_tibble(tbl)
  .dir_create_safe(dirname(path))
  utils::write.csv(tbl, file = path, row.names = FALSE, na = na)
  TRUE
}

.extract_tbl_export <- function(x, path = NULL) {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) return(tibble::as_tibble(x))

  if (is.list(x) && !inherits(x, "data.frame") && !is.null(path)) {
    obj <- x
    for (nm in path) {
      if (is.null(obj[[nm]])) return(NULL)
      obj <- obj[[nm]]
    }
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
  }

  NULL
}

.flatten_list_for_manifest <- function(x, prefix = NULL) {
  if (is.null(x)) return(list())

  if (!is.list(x) || inherits(x, "data.frame")) {
    nm <- if (is.null(prefix)) "value" else prefix

    val <- if (length(x) == 0) {
      NA_character_
    } else if (length(x) == 1) {
      as.character(x)
    } else {
      paste(as.character(x), collapse = " | ")
    }

    return(stats::setNames(list(val), nm))
  }

  out <- list()
  nms <- names(x)
  if (is.null(nms)) nms <- seq_along(x)

  for (i in seq_along(x)) {
    nm <- as.character(nms[i])
    key <- if (is.null(prefix)) nm else paste(prefix, nm, sep = ".")
    out <- c(out, .flatten_list_for_manifest(x[[i]], prefix = key))
  }

  out
}

.record_export <- function(log_tbl, category, name, path = NA_character_, status = "written", detail = NA_character_) {
  dplyr::bind_rows(
    log_tbl,
    tibble::tibble(
      category = as.character(category),
      name = as.character(name),
      path = as.character(path),
      status = as.character(status),
      detail = as.character(detail)
    )
  )
}

.empty_export_log <- function() {
  tibble::tibble(
    category = character(),
    name = character(),
    path = character(),
    status = character(),
    detail = character()
  )
}

# -------------------------------------------------------------------------
# Manifest
# -------------------------------------------------------------------------

#' Make a run manifest table
#'
#' @param run_name Name of the run.
#' @param prepared Optional prepared object.
#' @param alloc Optional allocation object.
#' @param analysis Optional analysis object.
#' @param ens Optional ensemble object.
#' @param final Optional finalized ensemble object.
#' @param significance Optional significance object or named list of such objects.
#'
#' @return One-row tibble manifest.
make_run_manifest <- function(run_name = "run",
                              prepared = NULL,
                              alloc = NULL,
                              analysis = NULL,
                              ens = NULL,
                              final = NULL,
                              significance = NULL) {
  ts <- as.character(Sys.time())

  out <- list(
    run_name = run_name,
    timestamp = ts,
    has_prepared = !is.null(prepared),
    has_allocation = !is.null(alloc),
    has_analysis = !is.null(analysis),
    has_ensemble = !is.null(ens),
    has_final = !is.null(final),
    has_significance = !is.null(significance)
  )

  if (!is.null(prepared) && !is.null(prepared$settings)) {
    out <- c(out, .flatten_list_for_manifest(prepared$settings, prefix = "prepared"))
  }

  if (!is.null(alloc) && !is.null(alloc$diagnostics)) {
    out <- c(out, .flatten_list_for_manifest(alloc$diagnostics, prefix = "allocation"))
  }

  if (!is.null(analysis) && !is.null(analysis$settings)) {
    out <- c(out, .flatten_list_for_manifest(analysis$settings, prefix = "analysis"))
  }

  if (!is.null(ens) && !is.null(ens$settings)) {
    out <- c(out, .flatten_list_for_manifest(ens$settings, prefix = "ensemble"))
  }

  if (!is.null(final) && !is.null(final$settings)) {
    out <- c(out, .flatten_list_for_manifest(final$settings, prefix = "final"))
  }

  tibble::as_tibble(out)
}

# -------------------------------------------------------------------------
# Table exporters
# -------------------------------------------------------------------------

#' Export a single table to CSV
#'
#' @param tbl Table to write.
#' @param path Destination path.
#' @param na String to use for missing values.
#'
#' @return Path invisibly.
export_table_csv <- function(tbl, path, na = "") {
  if (is.null(tbl)) {
    rlang::abort("`tbl` is NULL.")
  }
  .write_csv_safe(tbl, path, na = na)
  invisible(path)
}

#' Export core result tables
#'
#' @param out_dir Output directory.
#' @param prepared Optional prepared object.
#' @param alloc Optional allocation object.
#' @param analysis Optional analysis object.
#' @param ens Optional ensemble object.
#' @param final Optional finalized ensemble object.
#'
#' @return Export log tibble.
export_results_csv <- function(out_dir,
                               prepared = NULL,
                               alloc = NULL,
                               analysis = NULL,
                               ens = NULL,
                               final = NULL) {
  tables_dir <- file.path(out_dir, "tables")
  .dir_create_safe(tables_dir)

  log_tbl <- .empty_export_log()

  # prepared
  prepared_tables <- list(
    chronology = .extract_tbl_export(prepared, c("data", "chronology")),
    mortuary = .extract_tbl_export(prepared, c("data", "mortuary")),
    site_metadata = .extract_tbl_export(prepared, c("data", "site_metadata"))
  )

  for (nm in names(prepared_tables)) {
    tbl <- prepared_tables[[nm]]
    path <- file.path(tables_dir, paste0("prepared_", nm, ".csv"))
    if (!is.null(tbl)) {
      .write_csv_safe(tbl, path)
      log_tbl <- .record_export(log_tbl, "tables", paste0("prepared_", nm), path, "written")
    } else {
      log_tbl <- .record_export(log_tbl, "tables", paste0("prepared_", nm), NA, "skipped", "not available")
    }
  }

  # allocation
  alloc_tbl <- .extract_tbl_export(alloc, c("data"))
  if (!is.null(alloc_tbl)) {
    path <- file.path(tables_dir, "allocation.csv")
    .write_csv_safe(alloc_tbl, path)
    log_tbl <- .record_export(log_tbl, "tables", "allocation", path, "written")
  } else {
    log_tbl <- .record_export(log_tbl, "tables", "allocation", NA, "skipped", "not available")
  }

  # analysis
  analysis_tables <- list(
    site_bin = .extract_tbl_export(analysis, c("data", "site_bin")),
    bin_curve = .extract_tbl_export(analysis, c("data", "bin_curve")),
    region_bin = .extract_tbl_export(analysis, c("data", "region_bin"))
  )

  for (nm in names(analysis_tables)) {
    tbl <- analysis_tables[[nm]]
    path <- file.path(tables_dir, paste0("analysis_", nm, ".csv"))
    if (!is.null(tbl)) {
      .write_csv_safe(tbl, path)
      log_tbl <- .record_export(log_tbl, "tables", paste0("analysis_", nm), path, "written")
    } else {
      log_tbl <- .record_export(log_tbl, "tables", paste0("analysis_", nm), NA, "skipped", "not available")
    }
  }

  # ensemble
  ens_tbl <- .extract_tbl_export(ens, c("replicate_estimates"))
  if (!is.null(ens_tbl)) {
    path <- file.path(tables_dir, "ensemble_replicate_estimates.csv")
    .write_csv_safe(ens_tbl, path)
    log_tbl <- .record_export(log_tbl, "tables", "ensemble_replicate_estimates", path, "written")
  } else {
    log_tbl <- .record_export(log_tbl, "tables", "ensemble_replicate_estimates", NA, "skipped", "not available")
  }

  # final
  final_tables <- list(
    pooled_estimates = .extract_tbl_export(final, c("pooled", "estimates")),
    uncertainty_components = .extract_tbl_export(final, c("pooled", "components")),
    replicate_estimates = .extract_tbl_export(final, c("replicate_estimates"))
  )

  for (nm in names(final_tables)) {
    tbl <- final_tables[[nm]]
    path <- file.path(tables_dir, paste0("final_", nm, ".csv"))
    if (!is.null(tbl)) {
      .write_csv_safe(tbl, path)
      log_tbl <- .record_export(log_tbl, "tables", paste0("final_", nm), path, "written")
    } else {
      log_tbl <- .record_export(log_tbl, "tables", paste0("final_", nm), NA, "skipped", "not available")
    }
  }

  # fallback / uncertainty diagnostic
  if (!is.null(final)) {
    diag_fallback <- try(
      flag_uncertain_fallback_bins(
        final = final,
        provenance_source = final,
        fallback_threshold = 0.25,
        ratio_threshold = 1,
        normalized = TRUE
      ),
      silent = TRUE
    )

    if (!inherits(diag_fallback, "try-error") && !is.null(diag_fallback)) {
      path <- file.path(tables_dir, "final_fallback_uncertainty_diagnostic.csv")
      .write_csv_safe(diag_fallback, path)
      log_tbl <- .record_export(
        log_tbl,
        "tables",
        "final_fallback_uncertainty_diagnostic",
        path,
        "written"
      )
    } else {
      log_tbl <- .record_export(
        log_tbl,
        "tables",
        "final_fallback_uncertainty_diagnostic",
        NA,
        "skipped",
        "diagnostic could not be generated"
      )
    }
  } else {
    log_tbl <- .record_export(
      log_tbl,
      "tables",
      "final_fallback_uncertainty_diagnostic",
      NA,
      "skipped",
      "final not available"
    )
  }

  log_tbl
}

#' Export diagnostics and logs
#'
#' @param out_dir Output directory.
#' @param prepared Optional prepared object.
#' @param alloc Optional allocation object.
#' @param analysis Optional analysis object.
#' @param ens Optional ensemble object.
#' @param final Optional finalized ensemble object.
#'
#' @return Export log tibble.
export_diagnostics_csv <- function(out_dir,
                                   prepared = NULL,
                                   alloc = NULL,
                                   analysis = NULL,
                                   ens = NULL,
                                   final = NULL) {
  diag_dir <- file.path(out_dir, "diagnostics")
  .dir_create_safe(diag_dir)

  log_tbl <- .empty_export_log()

  write_diag_list <- function(x, prefix) {
    local_log <- .empty_export_log()
    if (is.null(x)) {
      local_log <- .record_export(local_log, "diagnostics", prefix, NA, "skipped", "not available")
      return(local_log)
    }

    if (!is.list(x)) {
      tbl <- tibble::tibble(value = as.character(x))
      path <- file.path(diag_dir, paste0(prefix, ".csv"))
      .write_csv_safe(tbl, path)
      local_log <- .record_export(local_log, "diagnostics", prefix, path, "written")
      return(local_log)
    }

    for (nm in names(x)) {
      obj <- x[[nm]]
      path <- file.path(diag_dir, paste0(prefix, "_", nm, ".csv"))

      if (is.data.frame(obj)) {
        .write_csv_safe(obj, path)
        local_log <- .record_export(local_log, "diagnostics", paste0(prefix, "_", nm), path, "written")
      } else if (is.list(obj)) {
        flat <- .flatten_list_for_manifest(obj)
        tbl <- tibble::tibble(
          key = names(flat),
          value = unname(unlist(flat))
        )
        .write_csv_safe(tbl, path)
        local_log <- .record_export(local_log, "diagnostics", paste0(prefix, "_", nm), path, "written")
      } else {
        tbl <- tibble::tibble(value = as.character(obj))
        .write_csv_safe(tbl, path)
        local_log <- .record_export(local_log, "diagnostics", paste0(prefix, "_", nm), path, "written")
      }
    }

    local_log
  }

  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(prepared$diagnostics, "prepared_diagnostics"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(prepared$logs, "prepared_logs"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(alloc$diagnostics, "allocation_diagnostics"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(alloc$logs, "allocation_logs"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(analysis$diagnostics, "analysis_diagnostics"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(ens$diagnostics, "ensemble_diagnostics"))
  log_tbl <- dplyr::bind_rows(log_tbl, write_diag_list(final$diagnostics, "final_diagnostics"))

  log_tbl
}

#' Export significance results
#'
#' @param out_dir Output directory.
#' @param significance A significance result object or named list of them.
#'
#' @return Export log tibble.
export_significance_csv <- function(out_dir,
                                    significance = NULL) {
  sig_dir <- file.path(out_dir, "tables", "significance")
  .dir_create_safe(sig_dir)

  log_tbl <- .empty_export_log()

  if (is.null(significance)) {
    return(.record_export(log_tbl, "tables", "significance", NA, "skipped", "not available"))
  }

  if (!is.list(significance) || ("result" %in% names(significance))) {
    significance <- list(significance = significance)
  }

  for (nm in names(significance)) {
    obj <- significance[[nm]]

    if (is.null(obj)) {
      log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm), NA, "skipped", "NULL")
      next
    }

    if (!is.null(obj$result)) {
      path1 <- file.path(sig_dir, paste0(nm, "_summary.csv"))
      .write_csv_safe(obj$result, path1)
      log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm, "_summary"), path1, "written")

      if (!is.null(obj$replicate_results)) {
        path2 <- file.path(sig_dir, paste0(nm, "_replicate_results.csv"))
        .write_csv_safe(obj$replicate_results, path2)
        log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm, "_replicate_results"), path2, "written")
      }

      if (!is.null(obj$bootstrap_distribution)) {
        path3 <- file.path(sig_dir, paste0(nm, "_bootstrap_distribution.csv"))
        .write_csv_safe(tibble::tibble(value = obj$bootstrap_distribution), path3)
        log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm, "_bootstrap_distribution"), path3, "written")
      }
    } else if (is.data.frame(obj)) {
      path <- file.path(sig_dir, paste0(nm, ".csv"))
      .write_csv_safe(obj, path)
      log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm), path, "written")
    } else {
      log_tbl <- .record_export(log_tbl, "tables", paste0("significance_", nm), NA, "skipped", "unsupported object")
    }
  }

  log_tbl
}

# -------------------------------------------------------------------------
# Plot export bundle
# -------------------------------------------------------------------------

#' Export standard plot bundle
#'
#' @param out_dir Output directory.
#' @param prepared Optional prepared object.
#' @param analysis Optional analysis object.
#' @param final Optional finalized ensemble object.
#' @param significance Optional significance object or list.
#' @param do_spatial Whether to export spatial plots.
#' @param spatial_bins Optional vector of bins for spatial slices.
#' @param spatial_bbox Optional bbox c(xmin, xmax, ymin, ymax).
#' @param region_col Region column name.
#'
#' @return Export log tibble.
export_plot_bundle <- function(out_dir,
                               prepared = NULL,
                               analysis = NULL,
                               final = NULL,
                               ens = NULL,
                               significance = NULL,
                               do_spatial = TRUE,
                               spatial_bins = NULL,
                               spatial_bbox = NULL,
                               region_col = "site_region",
                               spatial_buckets = NULL) {
  if (is.null(spatial_bins) && !is.null(spatial_buckets)) {
    spatial_bins <- spatial_buckets
  }
  plots_dir <- file.path(out_dir, "plots")
  .dir_create_safe(plots_dir)

  log_tbl <- .empty_export_log()

  # -----------------------------------------------------------------------
  # Harmonization
  # -----------------------------------------------------------------------
  if (!is.null(prepared)) {
    try({
      plot_harmonization_table_merged(
        chronology = prepared,
        save_path = file.path(plots_dir, "harmonization_table.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "harmonization_table",
        file.path(plots_dir, "harmonization_table.png"),
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "harmonization_table",
      NA, "skipped", "prepared not available"
    )
  }

  # -----------------------------------------------------------------------
  # Global pooled occupancy curve + uncertainty
  # -----------------------------------------------------------------------
  if (!is.null(final)) {
    try({
      plot_occupancy_curve(
        pooled_estimates = final,
        save_path = file.path(plots_dir, "occupancy_curve.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "occupancy_curve",
        file.path(plots_dir, "occupancy_curve.png"),
        "written"
      )
    }, silent = TRUE)

    try({
      plot_uncertainty_components(
        components = final,
        save_path = file.path(plots_dir, "uncertainty_components.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "uncertainty_components",
        file.path(plots_dir, "uncertainty_components.png"),
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "occupancy_curve",
      NA, "skipped", "final not available"
    )
    log_tbl <- .record_export(
      log_tbl, "plots", "uncertainty_components",
      NA, "skipped", "final not available"
    )
  }

  # -----------------------------------------------------------------------
  # Regional curves
  # Prefer pooled regional results from final; otherwise use analysis
  # -----------------------------------------------------------------------
  regional_source <- NULL
  if (!is.null(final) && !is.null(final$pooled) && !is.null(final$pooled$region_estimates)) {
    regional_source <- final
  } else if (!is.null(analysis)) {
    regional_source <- analysis
  }

  if (!is.null(regional_source)) {
    try({
      plot_regional_curves(
        region_bin = regional_source,
        region_col = region_col,
        facet = TRUE,
        save_path = file.path(plots_dir, "regional_curves_faceted.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "regional_curves_faceted",
        file.path(plots_dir, "regional_curves_faceted.png"),
        "written"
      )
    }, silent = TRUE)

    try({
      plot_regional_curves(
        region_bin = regional_source,
        region_col = region_col,
        facet = FALSE,
        save_dir = file.path(plots_dir, "regions")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "regional_curves_by_region",
        file.path(plots_dir, "regions"),
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "regional_curves_faceted",
      NA, "skipped", "neither final pooled regional results nor analysis available"
    )
    log_tbl <- .record_export(
      log_tbl, "plots", "regional_curves_by_region",
      NA, "skipped", "neither final pooled regional results nor analysis available"
    )
  }

  # -----------------------------------------------------------------------
  # Provenance contributions
  # Prefer final so contributions are averaged across replicates
  # -----------------------------------------------------------------------
  provenance_source <- NULL
  if (!is.null(final) && !is.null(final$replicate_data) && !is.null(final$replicate_data$site_bin)) {
    provenance_source <- final
  } else if (!is.null(ens) && !is.null(ens$replicate_results)) {
    provenance_source <- ens
  } else if (!is.null(analysis)) {
    provenance_source <- analysis
  }

  if (!is.null(provenance_source)) {
    try({
      plot_provenance_contributions(
        site_bin = provenance_source,
        normalized = TRUE,
        save_path = file.path(plots_dir, "provenance_contributions.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "provenance_contributions",
        file.path(plots_dir, "provenance_contributions.png"),
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "provenance_contributions",
      NA, "skipped", "neither final replicate site_bin, ensemble, nor analysis available"
    )
  }

  # -----------------------------------------------------------------------
  # Fallback / uncertainty diagnostic
  # -----------------------------------------------------------------------
  if (!is.null(final)) {
    try({
      diag_fallback <- flag_uncertain_fallback_bins(
        final = final,
        provenance_source = final,
        fallback_threshold = 0.25,
        ratio_threshold = 1,
        normalized = TRUE
      )

      plot_fallback_uncertainty_diagnostic(
        diag_tbl = diag_fallback,
        save_path = file.path(plots_dir, "fallback_uncertainty_diagnostic.png")
      )

      log_tbl <- .record_export(
        log_tbl, "plots", "fallback_uncertainty_diagnostic",
        file.path(plots_dir, "fallback_uncertainty_diagnostic.png"),
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "fallback_uncertainty_diagnostic",
      NA, "skipped", "final not available"
    )
  }

  # -----------------------------------------------------------------------
  # Spatial outputs
  # Prefer final, then ens, then analysis
  # -----------------------------------------------------------------------
  spatial_source <- NULL
  if (!is.null(final) && !is.null(final$replicate_data) && !is.null(final$replicate_data$site_bin)) {
    spatial_source <- final
  } else if (!is.null(ens) && !is.null(ens$replicate_results)) {
    spatial_source <- ens
  } else if (!is.null(analysis)) {
    spatial_source <- analysis
  }

  if (isTRUE(do_spatial) && !is.null(spatial_source) && !is.null(spatial_bins)) {
    try({
      plot_spatial_slices(
        site_bin = spatial_source,
        selected_bins = spatial_bins,
        bbox = spatial_bbox,
        facet = TRUE,
        save_path = file.path(plots_dir, "spatial_slices_faceted.png")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "spatial_slices_faceted",
        file.path(plots_dir, "spatial_slices_faceted.png"),
        "written"
      )
    }, silent = TRUE)

    try({
      plot_spatial_slices(
        site_bin = spatial_source,
        selected_bins = spatial_bins,
        bbox = spatial_bbox,
        facet = FALSE,
        save_dir = file.path(plots_dir, "spatial")
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "spatial_slices_by_bin",
        file.path(plots_dir, "spatial"),
        "written"
      )
    }, silent = TRUE)

    try({
      anim_path <- file.path(plots_dir, "spatial_density.gif")
      animate_spatial_density(
        site_bin = spatial_source,
        bbox = spatial_bbox,
        save_path = anim_path
      )
      log_tbl <- .record_export(
        log_tbl, "plots", "spatial_density_animation",
        anim_path,
        "written"
      )
    }, silent = TRUE)
  } else {
    log_tbl <- .record_export(
      log_tbl, "plots", "spatial_slices_faceted",
      NA, "skipped", "spatial source or spatial_bins not available"
    )
    log_tbl <- .record_export(
      log_tbl, "plots", "spatial_slices_by_bin",
      NA, "skipped", "spatial source or spatial_bins not available"
    )
    log_tbl <- .record_export(
      log_tbl, "plots", "spatial_density_animation",
      NA, "skipped", "spatial source or spatial_bins not available"
    )
  }

  log_tbl
}

# -------------------------------------------------------------------------
# Full run bundle
# -------------------------------------------------------------------------

#' Export a full run bundle
#'
#' @param out_dir Root output directory.
#' @param run_name Name of the run.
#' @param prepared Optional prepared object.
#' @param alloc Optional allocation object.
#' @param analysis Optional analysis object.
#' @param ens Optional ensemble object.
#' @param final Optional finalized ensemble object.
#' @param significance Optional significance result object or named list.
#' @param export_plots Whether to export plots.
#' @param do_spatial Whether to export spatial plots.
#' @param spatial_bins Optional vector of bins for spatial slices.
#' @param spatial_bbox Optional bbox c(xmin, xmax, ymin, ymax).
#' @param region_col Region column name.
#'
#' @return Structured list with export paths, manifest, and export log.
export_run_bundle <- function(out_dir,
                              run_name = "run",
                              prepared = NULL,
                              alloc = NULL,
                              analysis = NULL,
                              ens = NULL,
                              final = NULL,
                              significance = NULL,
                              export_plots = TRUE,
                              do_spatial = TRUE,
                              spatial_bins = NULL,
                              spatial_bbox = NULL,
                              region_col = "site_region",
                              spatial_buckets = NULL) {
  if (is.null(spatial_bins) && !is.null(spatial_buckets)) {
    spatial_bins <- spatial_buckets
  }
  run_dir <- file.path(out_dir, run_name)
  .dir_create_safe(run_dir)
  .dir_create_safe(file.path(run_dir, "tables"))
  .dir_create_safe(file.path(run_dir, "diagnostics"))
  .dir_create_safe(file.path(run_dir, "plots"))
  .dir_create_safe(file.path(run_dir, "manifests"))

  export_log <- .empty_export_log()

  # Manifest
  manifest <- make_run_manifest(
    run_name = run_name,
    prepared = prepared,
    alloc = alloc,
    analysis = analysis,
    ens = ens,
    final = final,
    significance = significance
  )

  manifest_path <- file.path(run_dir, "manifests", "run_manifest.csv")
  .write_csv_safe(manifest, manifest_path)
  export_log <- .record_export(export_log, "manifest", "run_manifest", manifest_path, "written")

  # Tables
  export_log <- dplyr::bind_rows(
    export_log,
    export_results_csv(
      out_dir = run_dir,
      prepared = prepared,
      alloc = alloc,
      analysis = analysis,
      ens = ens,
      final = final
    )
  )

  # Diagnostics
  export_log <- dplyr::bind_rows(
    export_log,
    export_diagnostics_csv(
      out_dir = run_dir,
      prepared = prepared,
      alloc = alloc,
      analysis = analysis,
      ens = ens,
      final = final
    )
  )

  # Significance
  export_log <- dplyr::bind_rows(
    export_log,
    export_significance_csv(
      out_dir = run_dir,
      significance = significance
    )
  )

  # Plots
  if (isTRUE(export_plots)) {
    export_log <- dplyr::bind_rows(
      export_log,
      export_plot_bundle(
        out_dir = run_dir,
        prepared = prepared,
        analysis = analysis,
        final = final,
        ens = ens,
        significance = significance,
        do_spatial = do_spatial,
        spatial_bins = spatial_bins,
        spatial_bbox = spatial_bbox,
        region_col = region_col
      )
    )
  } else {
    export_log <- .record_export(export_log, "plots", "plot_bundle", NA, "skipped", "export_plots = FALSE")
  }

  # Export log itself
  log_path <- file.path(run_dir, "manifests", "export_log.csv")
  .write_csv_safe(export_log, log_path)

  list(
    run_dir = run_dir,
    manifest = manifest,
    export_log = export_log,
    paths = list(
      manifest = manifest_path,
      export_log = log_path
    )
  )
}
