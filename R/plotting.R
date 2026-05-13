# plotting.R
#
# Plotting helpers for OccuPast.
#
# Expected inputs come from:
# - prepare_data.R
# - analyze_curve.R
# - uncertainty.R
# - ensemble.R
#
# Design:
# - ggplot2-first
# - OccuPast-like palette and minimal style
# - global, regional, diagnostic, and spatial plots separated
# - optional direct saving via ggplot2::ggsave()
# - animation via gganimate if installed

# Suggested imports in DESCRIPTION:
# Imports:
#   dplyr,
#   tibble,
#   rlang,
#   ggplot2
#
# Optional Suggests:
#   gganimate,
#   sf,
#   maps
#
# Suggested namespace usage:
#   @importFrom dplyr bind_rows mutate select filter group_by summarise
#   @importFrom tibble as_tibble
#   @importFrom rlang abort

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#OccuPast palette
.ocp_palette <- list(
  coral = "#E57050",
  teal = "#86B0A0",
  pale = "#CCDADE",
  dark = "#4A4A4A",
  red = "#d62728",
  orange = "#ff7f0e",
  green = "#2ca02c",
  grey = "#7f7f7f",
  point_red = "red3"
)

#OccuPast ggplot2 theme
.ocp_theme <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

.ocp_occupancy_curve_style <- function(p,
                                       title,
                                       subtitle = NULL,
                                       y_lab = "Normalized occupancy") {
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Pseudo-horizon bin",
      y = y_lab
    ) +
    .ocp_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "grey85", fill = NA),
      legend.position = "none"
    )
}

.ocp_spatial_map_style <- function(use_normalized = FALSE,
                                   point_colour = .ocp_palette$point_red,
                                   alpha_range = c(0.2, 0.7)) {

  list(
    ggplot2::scale_fill_viridis_c(option = "C"),

    ggplot2::scale_alpha(range = alpha_range, guide = "none"),

    ggplot2::scale_size_continuous(
      name = if (use_normalized) "Normalized occupancy" else "Expected burials",
      guide = ggplot2::guide_legend(
        override.aes = list(
          shape = 19,
          colour = point_colour,
          fill = "white",
          alpha = 1,
          stroke = 0
        )
      )
    ),

    ggplot2::labs(
      x = "Longitude",
      y = "Latitude",
      fill = "Density level"
    ),

    .ocp_theme(),

    ggplot2::theme(
      legend.position = "right",
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    ),

    ggplot2::guides(
      alpha = "none",
      fill = ggplot2::guide_colorbar(order = 1)
    )
  )
}

.get_world_basemap <- function(bbox = NULL) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    rlang::abort("Packages `rnaturalearth` and `sf` are required for the world basemap.")
  }

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  if (!is.null(bbox) && length(bbox) == 4L) {
    bb <- sf::st_bbox(c(
      xmin = bbox[1],
      xmax = bbox[2],
      ymin = bbox[3],
      ymax = bbox[4]
    ), crs = sf::st_crs(world))

    world <- sf::st_crop(world, bb)
  }

  world
}

.save_plot_optional <- function(plot_obj,
                                save_path = NULL,
                                width = 10,
                                height = 6,
                                dpi = 300,
                                bg = "white") {
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = save_path,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = bg
    )
  }
  invisible(plot_obj)
}

#' Flag bins dominated by between-run variance and site fallback contribution
#'
#' @description
#' Joins pooled uncertainty components to provenance contributions by bin and
#' identifies bins where between-replicate variance exceeds within-replicate
#' variance and site-fallback contribution is high.
#'
#' This is useful for diagnosing intervals where uncertainty is driven mainly by
#' weak chronological resolution, especially allocations based on
#' `site_start/site_end` fallback rather than phase-level chronology.
#'
#' @param final Finalized ensemble object from `finalize_ensemble()`, or a pooled
#'   uncertainty object containing `pooled$estimates`.
#' @param provenance_source Object usable by `plot_provenance_contributions()`,
#'   typically `final`, `ens`, or `analysis`.
#' @param fallback_threshold Minimum fallback share required to flag a bin.
#'   Default is 0.25.
#' @param ratio_threshold Minimum `B/W` ratio required to flag a bin.
#'   Default is 1.
#' @param normalized If TRUE, use normalized provenance columns where available.
#'
#' @return Tibble with bin-level uncertainty/provenance diagnostics and a
#'   logical `flagged` column.
flag_uncertain_fallback_bins <- function(final,
                                         provenance_source = final,
                                         fallback_threshold = 0.25,
                                         ratio_threshold = 1,
                                         normalized = TRUE) {
  pooled <- .extract_tbl(final, path = c("pooled", "estimates"))
  if (is.null(pooled)) {
    pooled <- .extract_tbl(final)
  }

  if (is.null(pooled)) {
    rlang::abort("Could not extract pooled estimates from `final`.")
  }

  .require_cols_plot(
    pooled,
    c("horizon_bin", "W", "B", "T"),
    "pooled_estimates"
  )

  # Use the same replicate-aware site-bin extraction used for spatial/provenance plotting
  sb <- .extract_site_bin_for_spatial_plot(provenance_source)

  .require_cols_plot(
    sb,
    c("horizon_bin"),
    "site_bin"
  )

  choose_cols <- function(tbl, normalized = TRUE) {
    if (normalized &&
        all(c("value_phase_system_norm",
              "value_site_fallback_norm",
              "value_synthetic_norm") %in% names(tbl))) {
      return(list(
        phase = "value_phase_system_norm",
        fallback = "value_site_fallback_norm",
        synthetic = "value_synthetic_norm"
      ))
    }

    if (all(c("value_phase_system",
              "value_site_fallback",
              "value_synthetic") %in% names(tbl))) {
      return(list(
        phase = "value_phase_system",
        fallback = "value_site_fallback",
        synthetic = "value_synthetic"
      ))
    }

    rlang::abort("Required provenance contribution columns were not found.")
  }

  cols <- choose_cols(sb, normalized = normalized)

  prov <- sb |>
    dplyr::group_by(horizon_bin) |>
    dplyr::summarise(
      phase_system = sum(.data[[cols$phase]], na.rm = TRUE),
      site_fallback = sum(.data[[cols$fallback]], na.rm = TRUE),
      synthetic = sum(.data[[cols$synthetic]], na.rm = TRUE),
      .groups = "drop"
    )

  prov$total_contribution <- prov$phase_system + prov$site_fallback + prov$synthetic
  prov$fallback_share <- ifelse(
    prov$total_contribution > 0,
    prov$site_fallback / prov$total_contribution,
    NA_real_
  )

  out <- dplyr::left_join(
    pooled[, c("horizon_bin", "W", "B", "T"), drop = FALSE],
    prov[, c("horizon_bin", "phase_system", "site_fallback", "synthetic", "total_contribution", "fallback_share"), drop = FALSE],
    by = "horizon_bin"
  )

  out$BW_ratio <- ifelse(is.finite(out$W) & out$W > 0, out$B / out$W, NA_real_)
  out$BT_share <- ifelse(is.finite(out$T) & out$T > 0, out$B / out$T, NA_real_)

  out$flag_B_gt_W <- is.finite(out$B) & is.finite(out$W) & out$B > out$W
  out$flag_fallback_high <- is.finite(out$fallback_share) & out$fallback_share >= fallback_threshold
  out$flag_ratio_high <- is.finite(out$BW_ratio) & out$BW_ratio >= ratio_threshold

  out$flagged <- out$flag_B_gt_W & out$flag_fallback_high & out$flag_ratio_high

  out |>
    dplyr::arrange(dplyr::desc(flagged), dplyr::desc(BW_ratio), dplyr::desc(fallback_share)) |>
    tibble::as_tibble()
}

.extract_tbl <- function(x, path = NULL) {
  if (is.null(x)) return(NULL)

  if (is.data.frame(x)) return(tibble::as_tibble(x))

  if (is.list(x) && !is.null(path)) {
    obj <- x
    for (nm in path) {
      if (is.null(obj[[nm]])) {
        return(NULL)
      }
      obj <- obj[[nm]]
    }
    if (is.data.frame(obj)) return(tibble::as_tibble(obj))
  }

  if (is.list(x) && !inherits(x, "data.frame") && "data" %in% names(x)) {
    if (is.data.frame(x$data)) return(tibble::as_tibble(x$data))
  }

  NULL
}

.extract_site_bin_replicates_for_plot <- function(x) {
  # finalized ensemble
  if (is.list(x) && !inherits(x, "data.frame") &&
      "replicate_data" %in% names(x) &&
      "site_bin" %in% names(x$replicate_data)) {
    return(lapply(x$replicate_data$site_bin, tibble::as_tibble))
  }

  # raw ensemble
  if (is.list(x) && !inherits(x, "data.frame") &&
      "replicate_results" %in% names(x)) {
    reps <- lapply(x$replicate_results, function(rr) {
      if (!is.null(rr$analysis) &&
          !is.null(rr$analysis$data) &&
          !is.null(rr$analysis$data$site_bin)) {
        tibble::as_tibble(rr$analysis$data$site_bin)
      } else {
        NULL
      }
    })
    reps <- Filter(Negate(is.null), reps)
    if (length(reps) > 0) return(reps)
  }

  # single analysis object
  sb <- .extract_tbl(x, path = c("data", "site_bin"))
  if (!is.null(sb)) return(list(tibble::as_tibble(sb)))

  # plain table
  if (is.data.frame(x)) return(list(tibble::as_tibble(x)))

  NULL
}

site_bin_to_canonical_grid <- function(site_bin_tbl, canonical_grid) {
  site_bin_tbl <- tibble::as_tibble(site_bin_tbl)
  canonical_grid <- tibble::as_tibble(canonical_grid)

  .require_cols_plot(
    site_bin_tbl,
    c("site_id", "bin_start", "bin_end"),
    "site_bin_tbl"
  )

  .require_cols_plot(
    canonical_grid,
    c("horizon_bin", "bin_start", "bin_end"),
    "canonical_grid"
  )

  numeric_cols <- names(site_bin_tbl)[vapply(site_bin_tbl, is.numeric, logical(1))]
  exclude_numeric <- c("bin_start", "bin_end", "horizon_bin", "replicate_id")
  value_cols <- setdiff(numeric_cols, exclude_numeric)

  context_cols <- intersect(
    c(
      "site_name", "coord_y", "coord_x", "site_country", "site_region",
      "site_admin", "site_start", "site_end", "site_size", "site_dig_date",
      "normalization_mode"
    ),
    names(site_bin_tbl)
  )

  first_non_missing <- function(x) {
    idx <- which(!is.na(x) & x != "")
    if (length(idx) == 0) {
      if (length(x) == 0) return(NA)
      return(x[1])
    }
    x[idx[1]]
  }

  sites <- unique(site_bin_tbl$site_id)
  out_all <- vector("list", length(sites))

  for (s in seq_along(sites)) {
    sid <- sites[s]
    stbl <- site_bin_tbl[site_bin_tbl$site_id == sid, , drop = FALSE]

    out_parts <- vector("list", nrow(canonical_grid))

    for (j in seq_len(nrow(canonical_grid))) {
      c_start <- canonical_grid$bin_start[j]
      c_end   <- canonical_grid$bin_end[j]

      overlaps <- pmax(
        0,
        pmin(stbl$bin_end, c_end) - pmax(stbl$bin_start, c_start)
      )

      keep <- overlaps > 0

      if (!any(keep)) {
        row_j <- tibble::tibble(
          site_id = sid,
          horizon_bin = canonical_grid$horizon_bin[j],
          bin_start = c_start,
          bin_end = c_end
        )

        for (nm in value_cols) row_j[[nm]] <- 0
        for (nm in context_cols) row_j[[nm]] <- first_non_missing(stbl[[nm]])

        out_parts[[j]] <- row_j
        next
      }

      widths <- stbl$bin_end[keep] - stbl$bin_start[keep]
      frac <- overlaps[keep] / widths

      row_j <- tibble::tibble(
        site_id = sid,
        horizon_bin = canonical_grid$horizon_bin[j],
        bin_start = c_start,
        bin_end = c_end
      )

      for (nm in value_cols) {
        row_j[[nm]] <- sum(stbl[[nm]][keep] * frac, na.rm = TRUE)
      }

      for (nm in context_cols) {
        row_j[[nm]] <- first_non_missing(stbl[[nm]])
      }

      out_parts[[j]] <- row_j
    }

    out_all[[s]] <- dplyr::bind_rows(out_parts)
  }

  dplyr::bind_rows(out_all)
}

summarize_site_bin_across_replicates <- function(site_bin_replicates,
                                                 canonical_grid = NULL) {
  if (is.null(site_bin_replicates) || length(site_bin_replicates) == 0) {
    rlang::abort("No site-bin replicate tables available.")
  }

  reps <- lapply(seq_along(site_bin_replicates), function(i) {
    sb <- tibble::as_tibble(site_bin_replicates[[i]])

    if (!is.null(canonical_grid)) {
      sb <- site_bin_to_canonical_grid(sb, canonical_grid)
    }

    sb$replicate_id <- i
    sb
  })

  dat <- dplyr::bind_rows(reps)

  .require_cols_plot(
    dat,
    c("site_id", "horizon_bin", "bin_start", "bin_end"),
    "site_bin_replicates"
  )

  numeric_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
  exclude_numeric <- c("replicate_id")
  avg_cols <- setdiff(numeric_cols, exclude_numeric)

  context_cols <- intersect(
    c(
      "site_name", "coord_y", "coord_x", "site_country", "site_region",
      "site_admin", "site_start", "site_end", "site_size", "site_dig_date",
      "normalization_mode"
    ),
    names(dat)
  )

  first_non_missing <- function(x) {
    idx <- which(!is.na(x) & x != "")
    if (length(idx) == 0) {
      if (length(x) == 0) return(NA)
      return(x[1])
    }
    x[idx[1]]
  }

  avg_exprs <- stats::setNames(
    lapply(avg_cols, function(nm) rlang::expr(mean(.data[[!!nm]], na.rm = TRUE))),
    avg_cols
  )

  context_exprs <- stats::setNames(
    lapply(context_cols, function(nm) rlang::expr(first_non_missing(.data[[!!nm]]))),
    context_cols
  )

  dat |>
    dplyr::group_by(site_id, horizon_bin, bin_start, bin_end) |>
    dplyr::summarise(
      !!!avg_exprs,
      !!!context_exprs,
      n_replicates = dplyr::n_distinct(replicate_id),
      .groups = "drop"
    ) |>
    tibble::as_tibble()
}

.extract_site_bin_for_spatial_plot <- function(x) {
  reps <- .extract_site_bin_replicates_for_plot(x)

  if (is.null(reps)) {
    rlang::abort(
      "Could not extract site-bin data for spatial plotting. Supply an analysis object, an ensemble object, a finalized ensemble object, or a site_bin table."
    )
  }

  if (length(reps) == 1L) {
    return(tibble::as_tibble(reps[[1]]))
  }

  canonical_grid <- NULL
  if (is.list(x) && !is.null(x$pooled) && !is.null(x$pooled$canonical_grid)) {
    canonical_grid <- tibble::as_tibble(x$pooled$canonical_grid)
  }

  summarize_site_bin_across_replicates(reps, canonical_grid = canonical_grid)
}
.require_cols_plot <- function(tbl, cols, name = "table") {
  validate_required_fields(tbl, cols, name)
  invisible(TRUE)
}

#' Function for plotting the already constructed harmonization table for
#' publication purposes.

plot_harmonization_table_merged <- function(chronology,
                                            save_path = NULL,
                                            base_size = 11,
                                            drop_systems = NULL,
                                            alternate_bands = TRUE,
                                            band_alpha = 0.15,
                                            show_grid = TRUE,
                                            grid_width = 10,
                                            grid_offset = 0,
                                            wrap_labels = TRUE,
                                            wrap_phu_per_char = 1.8,
                                            bar_linewidth = 4.5,
                                            fade_linewidth = 1.0,
                                            font_family = "",
                                            out_width = NULL,
                                            out_height = NULL,
                                            stagger_labels = TRUE,
                                            y_stagger = 0.18,
                                            label_size = 3.1,
                                            label_halo_size = 3.6,
                                            hide_short_labels = FALSE,
                                            min_label_span = 6) {
  chr <- .extract_tbl(chronology, path = c("data", "chronology"))
  if (is.null(chr)) chr <- .extract_tbl(chronology)

  if (is.null(chr)) {
    rlang::abort("Could not extract a chronology table.")
  }

  .require_cols_plot(
    chr,
    c("phase_id", "system_name", "phase_name",
      "horizon_start", "horizon_end", "fade_in_start", "fade_out_end"),
    "chronology"
  )

  if (!is.null(drop_systems) && length(drop_systems) > 0) {
    chr <- chr |>
      dplyr::filter(!.data$system_name %in% drop_systems)
  }

  if (nrow(chr) == 0) {
    rlang::abort("No chronology rows remain after filtering.")
  }

  chr <- chr |>
    dplyr::mutate(
      horizon_start = as.numeric(.data$horizon_start),
      horizon_end = as.numeric(.data$horizon_end),
      fade_in_start = as.numeric(.data$fade_in_start),
      fade_out_end = as.numeric(.data$fade_out_end),
      phase_name = as.character(.data$phase_name)
    ) |>
    dplyr::filter(is.finite(.data$horizon_start), is.finite(.data$horizon_end))

  if (nrow(chr) == 0) {
    rlang::abort("No valid chronology rows with finite horizon_start/horizon_end.")
  }

  sys_order <- chr |>
    dplyr::group_by(.data$system_name) |>
    dplyr::summarise(
      earliest_start = min(.data$horizon_start, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$earliest_start, .data$system_name) |>
    dplyr::mutate(system_numeric = dplyr::row_number())

  dat <- chr |>
    dplyr::left_join(sys_order, by = "system_name") |>
    dplyr::group_by(.data$system_name) |>
    dplyr::arrange(.data$horizon_start, .by_group = TRUE) |>
    dplyr::mutate(
      phase_idx = dplyr::row_number(),
      span_phu = pmax(0, .data$horizon_end - .data$horizon_start),
      x_mid = (.data$horizon_start + .data$horizon_end) / 2
    ) |>
    dplyr::ungroup()

  if (wrap_labels) {
    char_budget <- pmax(4, floor(dat$span_phu / wrap_phu_per_char))
    dat <- dat |>
      dplyr::mutate(
        label = purrr::map2_chr(
          .data$phase_name,
          char_budget,
          ~ stringr::str_wrap(.x %||% "", width = as.integer(.y))
        )
      )
  } else {
    dat <- dat |>
      dplyr::mutate(label = .data$phase_name)
  }

  if (hide_short_labels) {
    dat <- dat |>
      dplyr::mutate(
        label = dplyr::if_else(.data$span_phu < min_label_span, "", .data$label)
      )
  }

  dat <- dat |>
    dplyr::mutate(
      y_label = .data$system_numeric
    )

  if (stagger_labels) {
    dat <- dat |>
      dplyr::mutate(
        y_label = .data$system_numeric + ifelse(.data$phase_idx %% 2 == 0, y_stagger, -y_stagger)
      )
  }

  band_df <- sys_order |>
    dplyr::mutate(
      is_alt = (.data$system_numeric %% 2 == 0),
      y0 = .data$system_numeric - 0.5,
      y1 = .data$system_numeric + 0.5
    )

  gs <- min(dat$horizon_start, na.rm = TRUE)
  ge <- max(dat$horizon_end, na.rm = TRUE)

  b0 <- floor((gs - grid_offset) / grid_width) * grid_width + grid_offset
  b1 <- floor(((ge - 1e-9) - grid_offset) / grid_width) * grid_width + grid_offset
  grid_x <- seq(b0, b1, by = grid_width)

  sep <- dat |>
    dplyr::group_by(.data$system_name, .data$system_numeric) |>
    dplyr::arrange(.data$horizon_start, .by_group = TRUE) |>
    dplyr::mutate(next_start = dplyr::lead(.data$horizon_start)) |>
    dplyr::filter(!is.na(.data$next_start)) |>
    dplyr::transmute(
      system_numeric = .data$system_numeric,
      sep_x = .data$horizon_end
    ) |>
    dplyr::ungroup()

  p <- ggplot2::ggplot() +
    .ocp_theme(base_size = base_size, base_family = font_family) +
    ggplot2::labs(
      title = "Chronology harmonization",
      #subtitle = sprintf("Grid width=%s · offset=%s", grid_width, grid_offset),
      x = "Pseudo-horizon units (PHU)",
      y = "Systems"
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = base_size + 0.5, colour = .ocp_palette$dark),
      axis.text.x = ggplot2::element_text(size = base_size - 0.5, colour = .ocp_palette$dark),
      axis.title = ggplot2::element_text(size = base_size + 1, colour = .ocp_palette$dark),
      plot.title = ggplot2::element_text(face = "bold", colour = .ocp_palette$dark),
      plot.subtitle = ggplot2::element_text(colour = .ocp_palette$grey),
      panel.border = ggplot2::element_rect(colour = "grey80", fill = NA),
      plot.margin = grid::unit(c(6, 8, 6, 8), "pt")
    ) +
    ggplot2::scale_y_continuous(
      breaks = sys_order$system_numeric,
      labels = sys_order$system_name,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::coord_cartesian(xlim = c(gs, ge), clip = "off")

  if (alternate_bands && nrow(band_df) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = dplyr::filter(band_df, .data$is_alt),
        ggplot2::aes(xmin = gs, xmax = ge, ymin = y0, ymax = y1),
        inherit.aes = FALSE,
        fill = .ocp_palette$pale,
        alpha = band_alpha,
        colour = NA
      )
  }

  if (show_grid && length(grid_x) > 0) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = grid_x,
        colour = "grey85",
        linewidth = 0.3
      )
  }

  fade_in_dat <- dat |>
    dplyr::filter(is.finite(.data$fade_in_start))

  if (nrow(fade_in_dat) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = fade_in_dat,
        ggplot2::aes(
          x = .data$fade_in_start,
          xend = pmin(.data$horizon_start, .data$horizon_end),
          y = .data$system_numeric,
          yend = .data$system_numeric
        ),
        linewidth = fade_linewidth,
        colour = .ocp_palette$pale
      )
  }

  fade_out_dat <- dat |>
    dplyr::filter(is.finite(.data$fade_out_end))

  if (nrow(fade_out_dat) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = fade_out_dat,
        ggplot2::aes(
          x = pmax(.data$horizon_start, .data$horizon_end),
          xend = .data$fade_out_end,
          y = .data$system_numeric,
          yend = .data$system_numeric
        ),
        linewidth = fade_linewidth,
        colour = .ocp_palette$pale
      )
  }

  p <- p +
    ggplot2::geom_segment(
      data = dat,
      ggplot2::aes(
        x = .data$horizon_start,
        xend = .data$horizon_end,
        y = .data$system_numeric,
        yend = .data$system_numeric
      ),
      linewidth = bar_linewidth,
      colour = .ocp_palette$pale,
      lineend = "butt"
    )

  if (nrow(sep) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = sep,
        ggplot2::aes(
          x = .data$sep_x,
          xend = .data$sep_x,
          y = .data$system_numeric - 0.22,
          yend = .data$system_numeric + 0.22
        ),
        linewidth = 0.5,
        colour = .ocp_palette$teal
      )
  }

  label_dat <- dat |>
    dplyr::filter(.data$label != "")

  if (nrow(label_dat) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = label_dat,
        ggplot2::aes(
          x = .data$x_mid,
          y = .data$y_label,
          label = .data$label
        ),
        size = label_halo_size,
        colour = "white",
        lineheight = 0.95,
        family = font_family
      ) +
      ggplot2::geom_text(
        data = label_dat,
        ggplot2::aes(
          x = .data$x_mid,
          y = .data$y_label,
          label = .data$label
        ),
        size = label_size,
        fontface = "bold",
        colour = .ocp_palette$dark,
        lineheight = 0.95,
        family = font_family
      )
  }

  if (is.null(out_width) || is.null(out_height)) {
    n_sys <- nrow(sys_order)
    ph_per <- dat |>
      dplyr::count(.data$system_name, name = "k") |>
      dplyr::summarise(m = max(.data$k)) |>
      dplyr::pull(.data$m)

    ph_per <- ifelse(is.finite(ph_per), ph_per, 10)

    out_width <- out_width %||% (8 + 0.30 * ph_per)
    out_height <- out_height %||% (2.8 + 0.55 * n_sys)
  }

  .save_plot_optional(p, save_path = save_path, width = out_width, height = out_height)
  p
}
# -------------------------------------------------------------------------
# Global occupancy curve
# -------------------------------------------------------------------------

#' Plot pooled occupancy curve
#'
#' @param pooled_estimates Pooled estimates table from finalize_ensemble(), or tibble.
#' @param save_path Optional file path for saving.
#' @param show_points Show points.
#' @param show_ci Show CI ribbon/error bars if available.
#' @param add_loess Add descriptive LOESS smoother (not an inferential band).
#' @param loess_span LOESS span for descriptive smoothing.
#'
#' @return ggplot object.
plot_occupancy_curve <- function(pooled_estimates,
                                 save_path = NULL,
                                 show_points = TRUE,
                                 show_ci = TRUE,
                                 show_ci_ribbon = FALSE,
                                 show_line = FALSE,
                                 add_loess = TRUE,
                                 loess_span = 0.75,
                                 subtitle = NULL) {
  dat <- .extract_tbl(pooled_estimates, path = c("pooled", "estimates"))
  if (is.null(dat)) dat <- .extract_tbl(pooled_estimates)

  if (is.null(dat)) {
    rlang::abort("Could not extract pooled estimates.")
  }

  if ("estimate" %in% names(dat)) {
    ycol <- "estimate"
  } else if ("mean_value" %in% names(dat)) {
    ycol <- "mean_value"
  } else if ("total_value" %in% names(dat)) {
    ycol <- "total_value"
  } else {
    rlang::abort("No plottable estimate column found.")
  }

  .require_cols_plot(dat, c("horizon_bin", ycol), "pooled_estimates")

  dat <- tibble::as_tibble(dat) |>
    dplyr::mutate(horizon_bin = as.numeric(.data$horizon_bin)) |>
    dplyr::filter(is.finite(.data$horizon_bin), is.finite(.data[[ycol]]))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$horizon_bin, y = .data[[ycol]]))

  if (show_ci_ribbon && all(c("lower", "upper") %in% names(dat))) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        fill = .ocp_palette$pale,
        alpha = 0.35
      )
  }

  if (show_ci && all(c("lower", "upper") %in% names(dat))) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        width = 0.2,
        alpha = 0.7,
        colour = .ocp_palette$dark
      )
  }

  if (show_points) {
    p <- p + ggplot2::geom_point(colour = .ocp_palette$dark)
  }

  if (show_line) {
    p <- p + ggplot2::geom_line(linewidth = 0.8, colour = .ocp_palette$coral)
  }

  # Default ggplot LOESS bands describe the smoother, not the full
  # replicate uncertainty. Use replicate-wise LOESS bands when available.
  if (add_loess) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        se = TRUE,
        span = loess_span,
        colour = .ocp_palette$coral,
        fill = .ocp_palette$pale,
        alpha = 0.5
      )
  }

  p <- .ocp_occupancy_curve_style(
    p,
    title = "Normalized cemetery occupancy curve across regions",
    subtitle = subtitle,
    y_lab = "Normalized occupancy"
  )

  .save_plot_optional(p, save_path = save_path, width = 10, height = 6)
  p
}

# -------------------------------------------------------------------------
# Regional curves
# -------------------------------------------------------------------------

#' Plot regional curves
#'
#' @param region_bin Region-by-bin table, or object containing it.
#' @param region_col Region column name.
#' @param facet If TRUE, facet regions in one figure.
#' @param save_path Optional file path for saving the faceted plot.
#' @param save_dir Optional directory for one-file-per-region export.
#' @param add_loess Add descriptive LOESS smoother (not an inferential band).
#' @param loess_span LOESS span for descriptive smoothing.
#'
#' @return ggplot object or list of ggplots if facet = FALSE and save_dir is NULL.
plot_regional_curves <- function(region_bin,
                                 region_col = "site_region",
                                 facet = TRUE,
                                 save_path = NULL,
                                 save_dir = NULL,
                                 show_ci = TRUE,
                                 show_ci_ribbon = FALSE,
                                 show_points = TRUE,
                                 show_line = FALSE,
                                 add_loess = TRUE,
                                 loess_span = 0.75,
                                 subtitle = NULL) {
  dat <- .extract_tbl(region_bin, path = c("pooled", "region_estimates"))
  if (is.null(dat)) dat <- .extract_tbl(region_bin, path = c("data", "region_bin"))
  if (is.null(dat)) dat <- .extract_tbl(region_bin)

  if (is.null(dat)) {
    rlang::abort("Could not extract region_bin or pooled region_estimates.")
  }

  if (!region_col %in% names(dat)) {
    rlang::abort(paste0("Region column `", region_col, "` not found."))
  }

  if ("estimate" %in% names(dat)) {
    ycol <- "estimate"
  } else if ("mean_value" %in% names(dat)) {
    ycol <- "mean_value"
  } else if ("total_value" %in% names(dat)) {
    ycol <- "total_value"
  } else {
    rlang::abort("No plottable value column found for regional curves.")
  }

  .require_cols_plot(dat, c("horizon_bin", region_col, ycol), "region_bin")

  dat <- tibble::as_tibble(dat) |>
    dplyr::filter(!is.na(.data[[region_col]]), .data[[region_col]] != "") |>
    dplyr::mutate(horizon_bin = as.numeric(.data$horizon_bin)) |>
    dplyr::filter(is.finite(.data$horizon_bin), is.finite(.data[[ycol]]))

  if (nrow(dat) == 0) {
    rlang::abort("No non-missing regional rows available to plot.")
  }

  if (facet) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$horizon_bin, y = .data[[ycol]]))

    if (show_ci_ribbon && all(c("lower", "upper") %in% names(dat))) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          fill = .ocp_palette$pale,
          alpha = 0.35
        )
    }

    if (show_ci && all(c("lower", "upper") %in% names(dat))) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2,
          alpha = 0.7,
          colour = .ocp_palette$dark
        )
    }

    if (show_points) {
      p <- p + ggplot2::geom_point(size = 1.6, colour = .ocp_palette$dark)
    }

    if (show_line) {
      p <- p + ggplot2::geom_line(linewidth = 0.7, colour = .ocp_palette$teal)
    }

    if (add_loess) {
      p <- p +
        ggplot2::geom_smooth(
          method = "loess",
          se = TRUE,
          span = loess_span,
          colour = .ocp_palette$teal,
          fill = .ocp_palette$pale,
          alpha = 0.5
        )
    }

    p <- p +
      ggplot2::facet_wrap(stats::as.formula(paste("~", region_col)), scales = "free_y")

    p <- .ocp_occupancy_curve_style(
      p,
      title = "Regional normalized cemetery occupancy curve",
      subtitle = subtitle,
      y_lab = "Normalized occupancy"
    )

    .save_plot_optional(p, save_path = save_path, width = 12, height = 8)
    return(p)
  }

  regs <- sort(unique(dat[[region_col]]))

  plots <- lapply(regs, function(reg) {
    d <- dat[dat[[region_col]] == reg, , drop = FALSE]

    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$horizon_bin, y = .data[[ycol]]))

    if (show_ci_ribbon && all(c("lower", "upper") %in% names(d))) {
      p <- p +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          fill = .ocp_palette$pale,
          alpha = 0.35
        )
    }

    if (show_ci && all(c("lower", "upper") %in% names(d))) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
          width = 0.2,
          alpha = 0.7,
          colour = .ocp_palette$dark
        )
    }

    if (show_points) {
      p <- p + ggplot2::geom_point(size = 1.6, colour = .ocp_palette$dark)
    }

    if (show_line) {
      p <- p + ggplot2::geom_line(linewidth = 0.7, colour = .ocp_palette$teal)
    }

    if (add_loess) {
      p <- p +
        ggplot2::geom_smooth(
          method = "loess",
          se = TRUE,
          span = loess_span,
          colour = .ocp_palette$teal,
          fill = .ocp_palette$pale,
          alpha = 0.5
        )
    }

    p <- .ocp_occupancy_curve_style(
      p,
      title = paste("Regional normalized cemetery occupancy curve:", reg),
      subtitle = subtitle,
      y_lab = "Normalized occupancy"
    )

    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      fn <- file.path(
        save_dir,
        paste0("region_curve_", gsub("[^A-Za-z0-9_]+", "_", reg), ".png")
      )
      .save_plot_optional(p, save_path = fn, width = 10, height = 6)
    }

    p
  })

  names(plots) <- regs
  plots
}

#' Plot fallback share against between/within uncertainty ratio
#'
#' @param diag_tbl Output from `flag_uncertain_fallback_bins()`.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_fallback_uncertainty_diagnostic <- function(diag_tbl,
                                                 save_path = NULL,
                                                 bar_width = 9) {
  diag_tbl <- tibble::as_tibble(diag_tbl)

  .require_cols_plot(
    diag_tbl,
    c("horizon_bin", "fallback_share", "BW_ratio"),
    "diag_tbl"
  )

  p <- ggplot2::ggplot(
    diag_tbl,
    ggplot2::aes(
      x = horizon_bin,
      y = BW_ratio,
      fill = fallback_share
    )
  ) +
    ggplot2::geom_col(width = bar_width) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = .ocp_palette$dark
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      name = "Fallback share"
    ) +
    ggplot2::labs(
      title = "Fallback-driven uncertainty diagnostic",
      x = "Horizon bin",
      y = "Between / Within variance ratio (B/W)"
    ) +
    .ocp_theme()

  .save_plot_optional(p, save_path = save_path, width = 11, height = 6)

  p
}

# -------------------------------------------------------------------------
# Uncertainty components
# -------------------------------------------------------------------------

#' Plot uncertainty components
#'
#' @param components Output from summarize_uncertainty_components(), or tibble.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_uncertainty_components <- function(components,
                                        save_path = NULL) {
  dat <- .extract_tbl(components, path = c("pooled", "components"))
  if (is.null(dat)) dat <- .extract_tbl(components)

  if (is.null(dat)) {
    rlang::abort("Could not extract uncertainty components.")
  }

  .require_cols_plot(dat, c("horizon_bin", "component", "variance"), "components")

  comp_cols <- c(
    within = .ocp_palette$teal,
    between = .ocp_palette$orange,
    total = .ocp_palette$coral
  )

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = horizon_bin, y = variance, colour = component, fill = component)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_area(alpha = 0.18, position = "identity") +
    ggplot2::scale_colour_manual(values = comp_cols) +
    ggplot2::scale_fill_manual(values = comp_cols) +
    ggplot2::labs(
      title = "Uncertainty components",
      x = "Horizon bin",
      y = "Variance"
    ) +
    .ocp_theme()

  .save_plot_optional(p, save_path = save_path, width = 10, height = 6)
  p
}

# -------------------------------------------------------------------------
# Provenance contributions
# -------------------------------------------------------------------------

#' Plot provenance contributions across bins
#'
#' @param site_bin Site-bin table, or object containing it.
#' @param normalized If TRUE, use normalized provenance columns where present.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_provenance_contributions <- function(site_bin,
                                          normalized = TRUE,
                                          save_path = NULL) {

  # Use the same replicate-aware extractor used for spatial plotting
  sb <- .extract_site_bin_for_spatial_plot(site_bin)

  .require_cols_plot(sb, c("horizon_bin"), "site_bin")

  choose_cols <- function(tbl, normalized = TRUE) {

    if (normalized &&
        all(c("value_phase_system_norm",
              "value_site_fallback_norm",
              "value_synthetic_norm") %in% names(tbl))) {

      return(list(
        phase = "value_phase_system_norm",
        fallback = "value_site_fallback_norm",
        synthetic = "value_synthetic_norm"
      ))
    }

    if (all(c("value_phase_system",
              "value_site_fallback",
              "value_synthetic") %in% names(tbl))) {

      return(list(
        phase = "value_phase_system",
        fallback = "value_site_fallback",
        synthetic = "value_synthetic"
      ))
    }

    rlang::abort("Required provenance contribution columns were not found.")
  }

  cols <- choose_cols(sb, normalized = normalized)

  agg <- sb |>
    dplyr::group_by(horizon_bin) |>
    dplyr::summarise(
      phase_system = sum(.data[[cols$phase]], na.rm = TRUE),
      site_fallback = sum(.data[[cols$fallback]], na.rm = TRUE),
      synthetic = sum(.data[[cols$synthetic]], na.rm = TRUE),
      .groups = "drop"
    )

  long <- dplyr::bind_rows(
    tibble::tibble(
      horizon_bin = agg$horizon_bin,
      source = "phase_system",
      value = agg$phase_system
    ),
    tibble::tibble(
      horizon_bin = agg$horizon_bin,
      source = "site_fallback",
      value = agg$site_fallback
    ),
    tibble::tibble(
      horizon_bin = agg$horizon_bin,
      source = "synthetic",
      value = agg$synthetic
    )
  )

  src_cols <- c(
    phase_system = .ocp_palette$teal,
    site_fallback = .ocp_palette$orange,
    synthetic = .ocp_palette$red
  )

  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = horizon_bin,
                                    y = value,
                                    fill = source)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = src_cols) +
    ggplot2::labs(
      title = "Provenance contributions by bin",
      x = "Horizon bin",
      y = "Contribution"
    ) +
    .ocp_theme()

  .save_plot_optional(p, save_path = save_path, width = 10, height = 6)

  p
}

# -------------------------------------------------------------------------
# Spatial slices
# -------------------------------------------------------------------------

#' Plot spatial slices for selected bins
#'
#' @param site_bin Site-bin table, or object containing it.
#' @param selected_bins Numeric vector of bin IDs to show.
#' @param use_normalized If TRUE, prefer normalized value columns.
#' @param point_colour Point outline colour.
#' @param point_alpha Point alpha.
#' @param bbox Optional numeric vector c(xmin, xmax, ymin, ymax).
#' @param facet If TRUE, facet selected bins in one figure.
#' @param save_path Optional file path for saving the faceted figure.
#' @param save_dir Optional directory for one-file-per-bin output.
#'
#' @return ggplot object or list of ggplots.
plot_spatial_slices <- function(site_bin,
                                selected_bins,
                                use_normalized = FALSE,
                                point_colour = .ocp_palette$point_red,
                                point_alpha = 0.6,
                                bbox = NULL,
                                facet = TRUE,
                                save_path = NULL,
                                save_dir = NULL,
                                bins = 7) {
  sb <- .extract_site_bin_for_spatial_plot(site_bin)

  .require_cols_plot(sb, c("site_id", "horizon_bin", "coord_x", "coord_y"), "site_bin")

  value_col <- NULL
  if (use_normalized && "value_norm" %in% names(sb)) {
    value_col <- "value_norm"
  } else if ("value" %in% names(sb)) {
    value_col <- "value"
  } else if ("mean_value" %in% names(sb)) {
    value_col <- "mean_value"
  } else if ("estimate" %in% names(sb)) {
    value_col <- "estimate"
  } else {
    rlang::abort("No plottable value column found for spatial plotting.")
  }

  dat <- sb |>
    dplyr::filter(.data$horizon_bin %in% selected_bins) |>
    dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) |>
    dplyr::filter(is.finite(.data[[value_col]])) |>
    dplyr::filter(.data[[value_col]] > 0)

  if (nrow(dat) == 0) {
    rlang::abort("No rows available for the requested spatial bins.")
  }

  bg <- .get_world_basemap(bbox = bbox)

  density_aes <- ggplot2::aes(
    x = coord_x,
    y = coord_y,
    weight = .data[[value_col]],
    fill = ggplot2::after_stat(level),
    alpha = ggplot2::after_stat(level)
  )

  base_plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = bg,
      fill = "grey95",
      colour = "grey80",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::stat_density_2d(
      data = dat,
      mapping = density_aes,
      geom = "polygon",
      contour = TRUE,
      contour_var = "ndensity",
      bins = bins,
      show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
    ) +
    ggplot2::geom_point(
      data = dat,
      ggplot2::aes(x = coord_x, y = coord_y, size = .data[[value_col]]),
      colour = point_colour,
      alpha = point_alpha,
      show.legend = c(size = TRUE)
    ) +
    .ocp_spatial_map_style(
      use_normalized = use_normalized,
      point_colour = point_colour,
      alpha_range = c(0.2, 0.7)
    )

  if (!is.null(bbox) && length(bbox) == 4L) {
    base_plot <- base_plot +
      ggplot2::coord_sf(
        xlim = c(bbox[1], bbox[2]),
        ylim = c(bbox[3], bbox[4]),
        expand = FALSE
      )
  } else {
    base_plot <- base_plot + ggplot2::coord_sf(expand = FALSE)
  }

  if (facet) {
    p <- base_plot +
      ggplot2::facet_wrap(~ horizon_bin) +
      ggplot2::labs(title = "Spatial density by horizon bin")

    .save_plot_optional(p, save_path = save_path, width = 12, height = 8)
    return(p)
  }

  plots <- lapply(sort(unique(dat$horizon_bin)), function(hb) {
    d <- dat[dat$horizon_bin == hb, , drop = FALSE]

    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(
        data = bg,
        fill = "grey95",
        colour = "grey80",
        linewidth = 0.2,
        inherit.aes = FALSE
      ) +
      ggplot2::stat_density_2d(
        data = d,
        mapping = density_aes,
        geom = "polygon",
        contour = TRUE,
        contour_var = "ndensity",
        bins = bins,
        show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
      ) +
      ggplot2::geom_point(
        data = d,
        ggplot2::aes(x = coord_x, y = coord_y, size = .data[[value_col]]),
        colour = point_colour,
        alpha = point_alpha,
        show.legend = c(size = TRUE)
      ) +
      ggplot2::labs(title = sprintf("Horizon Bin: %s", hb)) +
      .ocp_spatial_map_style(
        use_normalized = use_normalized,
        point_colour = point_colour,
        alpha_range = c(0.2, 0.7)
      )

    if (!is.null(bbox) && length(bbox) == 4L) {
      p <- p +
        ggplot2::coord_sf(
          xlim = c(bbox[1], bbox[2]),
          ylim = c(bbox[3], bbox[4]),
          expand = FALSE
        )
    } else {
      p <- p + ggplot2::coord_sf(expand = FALSE)
    }

    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      fn <- file.path(save_dir, paste0("spatial_slice_", hb, ".png"))
      .save_plot_optional(p, save_path = fn, width = 8, height = 6)
    }

    p
  })

  names(plots) <- as.character(sort(unique(dat$horizon_bin)))
  plots
}

# -------------------------------------------------------------------------
# Spatial animation
# -------------------------------------------------------------------------

#' Animate spatial density over horizon bins
#'
#' @param site_bin Site-bin table, or object containing it.
#' @param use_normalized If TRUE, prefer normalized value columns.
#' @param bin_group_size Optional integer bin size to aggregate neighboring bins.
#' @param point_colour Point colour.
#' @param point_alpha Point alpha.
#' @param bbox Optional numeric vector c(xmin, xmax, ymin, ymax).
#' @param save_path Optional GIF/video path.
#' @param width Plot width in px for animation.
#' @param height Plot height in px for animation.
#' @param fps Frames per second.
#' @param duration Total animation duration in seconds.
#' @param end_pause Number of frames to pause on the final frame.
#' @param transition_length Relative transition length between states.
#' @param state_length Relative hold length within each state.
#' @param bins Number of density bins.
#'
#' @return Animated plot object.
animate_spatial_density <- function(site_bin,
                                    use_normalized = FALSE,
                                    bin_group_size = NULL,
                                    point_colour = .ocp_palette$point_red,
                                    point_alpha = 0.6,
                                    bbox = NULL,
                                    save_path = NULL,
                                    width = 900,
                                    height = 650,
                                    fps = 12,
                                    duration = 8,
                                    end_pause = 12,
                                    transition_length = 15,
                                    state_length = 30,
                                    bins = 7) {
  if (!requireNamespace("gganimate", quietly = TRUE)) {
    rlang::abort("`gganimate` is required for `animate_spatial_density()`.")
  }

  sb <- .extract_site_bin_for_spatial_plot(site_bin)

  .require_cols_plot(sb, c("site_id", "horizon_bin", "coord_x", "coord_y"), "site_bin")

  value_col <- NULL
  if (use_normalized && "value_norm" %in% names(sb)) {
    value_col <- "value_norm"
  } else if ("value" %in% names(sb)) {
    value_col <- "value"
  } else if ("mean_value" %in% names(sb)) {
    value_col <- "mean_value"
  } else if ("estimate" %in% names(sb)) {
    value_col <- "estimate"
  } else {
    rlang::abort("No plottable value column found for spatial animation.")
  }

  dat <- sb |>
    dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) |>
    dplyr::filter(is.finite(.data[[value_col]])) |>
    dplyr::filter(.data[[value_col]] > 0)

  if (nrow(dat) == 0) {
    rlang::abort("No rows available for spatial animation.")
  }

  if (!is.null(bin_group_size) && is.numeric(bin_group_size) && bin_group_size > 1) {
    dat <- dat |>
      dplyr::mutate(combined_bin = floor(.data$horizon_bin / bin_group_size) * bin_group_size)
  } else {
    dat <- dat |>
      dplyr::mutate(combined_bin = .data$horizon_bin)
  }

  dat <- dat |>
    dplyr::group_by(.data$combined_bin) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup()

  if (nrow(dat) == 0) {
    rlang::abort("Not enough points to animate density.")
  }

  bg <- .get_world_basemap(bbox = bbox)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = coord_x, y = coord_y)) +
    ggplot2::geom_sf(
      data = bg,
      fill = "grey95",
      colour = "grey80",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::stat_density_2d(
      ggplot2::aes(
        weight = .data[[value_col]],
        fill = ggplot2::after_stat(level),
        alpha = ggplot2::after_stat(level)
      ),
      geom = "polygon",
      contour = TRUE,
      contour_var = "ndensity",
      bins = bins,
      show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[[value_col]]),
      colour = point_colour,
      alpha = point_alpha,
      show.legend = c(size = TRUE)
    ) +
    .ocp_spatial_map_style(
      use_normalized = use_normalized,
      point_colour = point_colour,
      alpha_range = c(0.1, 0.9)
    ) +
    ggplot2::labs(title = "Horizon Bin: {closest_state}") +
    gganimate::transition_states(
      states = combined_bin,
      transition_length = transition_length,
      state_length = state_length
    ) +
    gganimate::ease_aes("circular-in-out")

  if (!is.null(bbox) && length(bbox) == 4L) {
    p <- p +
      ggplot2::coord_sf(
        xlim = c(bbox[1], bbox[2]),
        ylim = c(bbox[3], bbox[4]),
        expand = FALSE
      )
  } else {
    p <- p + ggplot2::coord_sf(expand = FALSE)
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)

    ext <- tolower(tools::file_ext(save_path))

    renderer <- switch(
      ext,
      gif = gganimate::gifski_renderer(save_path),
      mp4 = gganimate::ffmpeg_renderer(save_path),
      mov = gganimate::ffmpeg_renderer(save_path),
      avi = gganimate::ffmpeg_renderer(save_path),
      rlang::abort("Unsupported animation file extension. Use .gif or .mp4.")
    )

    anim <- gganimate::animate(
      p,
      width = width,
      height = height,
      fps = fps,
      duration = duration,
      end_pause = end_pause,
      renderer = renderer
    )
  } else {
    anim <- gganimate::animate(
      p,
      width = width,
      height = height,
      fps = fps,
      duration = duration,
      end_pause = end_pause
    )
  }

  anim
}


# -------------------------------------------------------------------------
# Legacy aliases
# -------------------------------------------------------------------------
# The package now uses "bin" terminology. These wrappers keep older scripts
# from failing immediately while making new output use the bin-based API.
flag_uncertain_fallback_buckets <- function(...) flag_uncertain_fallback_bins(...)
site_bucket_to_canonical_grid <- function(site_bucket_tbl, canonical_grid) {
  site_bin_to_canonical_grid(site_bucket_tbl, canonical_grid)
}
summarize_site_bucket_across_replicates <- function(site_bucket_replicates, canonical_grid = NULL) {
  summarize_site_bin_across_replicates(site_bucket_replicates, canonical_grid = canonical_grid)
}
