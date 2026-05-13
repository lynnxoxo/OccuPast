# mortAAR_bridge.R
#
# Bridge helpers between OccuPast and mortAAR.
#
# Purpose:
# - prepare burial-level mortAAR input from ensemble output
# - compute site-level mortAAR sex-relation metrics
# - compute the same metrics across sliding temporal windows
# - optionally build mortAAR life-table objects for those windows
#
# Design:
# - Time assignment is based on `median_bin` produced by prepare_for_mortAAR().
# - Indeterminate-sex burials are RETAINED in the burial table and general
#   life-table computations, but sex-relation metrics are computed only from
#   sexed individuals ("f" / "m").
# - Burial-level representative bins are snapped to the canonical grid stored
#   in `final$pooled$canonical_grid`.
# - mortAAR should not receive `group = NULL`; the `group` argument must be
#   omitted when no grouping is desired.
#
# Expected upstream functions:
# - validate_required_fields()
# - flag_uncertain_fallback_bins()

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Compute sex-stratified mortAAR life tables across sliding windows
#'
#' @description
#' Builds life tables by window with sex as the grouping variable. Only sexed
#' individuals ("f", "m") are used.
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#' @param window_width Width of each window in bin units.
#' @param step_size Step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#' @param method Method passed to `mortAAR::prep.life.table()`.
#'
#' @return Named list of grouped mortAAR life-table objects.
compute_life_tables_windows_by_sex <- function(final,
                                               prepared,
                                               window_width,
                                               step_size,
                                               start_bin = NULL,
                                               end_bin = NULL,
                                               method = "Standard") {
  .assert_package_mortAAR()

  mort_tbl <- build_mortaar_burial_table(final, prepared) %>%
    dplyr::filter(!is.na(sex_gender))

  windows <- make_mortaar_windows(
    mort_tbl = mort_tbl,
    window_width = window_width,
    step_size = step_size,
    start_bin = start_bin,
    end_bin = end_bin
  )

  out <- purrr::map(seq_len(nrow(windows)), function(i) {
    w <- windows[i, ]

    mort_tbl_w <- filter_mortaar_by_window(
      mort_tbl,
      window_start = w$window_start,
      window_end = w$window_end
    )

    if (nrow(mort_tbl_w) == 0) {
      return(NULL)
    }

    prep_obj <- .run_prep_life_table(
      x = mort_tbl_w,
      dec = NA,
      agebeg = "from",
      ageend = "to",
      group = "sex_gender",
      method = method,
      agerange = "included"
    )

    mortAAR::life.table(prep_obj)
  })

  names(out) <- paste0(windows$window_start, "_", windows$window_end)
  out
}

.split_grouped_life_table <- function(lt_obj) {
  if (is.null(lt_obj)) return(NULL)

  # Common mortAAR grouped output case: named list of tables
  if (is.list(lt_obj) && !inherits(lt_obj, "data.frame")) {
    keep <- vapply(lt_obj, function(x) {
      is.data.frame(x) || inherits(x, "data.frame")
    }, logical(1))

    if (any(keep)) {
      return(lt_obj[keep])
    }
  }

  # Fallback: data frame with group column
  lt_tbl <- try(as.data.frame(lt_obj), silent = TRUE)
  if (!inherits(lt_tbl, "try-error") && "group" %in% names(lt_tbl)) {
    return(split(lt_tbl, lt_tbl$group))
  }

  NULL
}

#' Extract one grouped life-table measure across sliding windows
#'
#' @param life_tables_windows_by_group Named list returned by
#'   `compute_life_tables_windows_by_sex()`.
#' @param measure One of "qx", "dx", "lx", "ex", or "rel_popx".
#' @param groups Character vector of groups to keep.
#'
#' @return Tibble with columns `window_label`, `window_start`, `window_end`,
#'   `window_mid`, `sex_gender`, `age`, and `value`.
extract_life_table_measure_windows_by_group <- function(life_tables_windows_by_group,
                                                        measure = c("qx", "dx", "lx", "ex", "rel_popx"),
                                                        groups = c("f", "m", "indet")) {
  measure <- match.arg(measure)

  if (is.null(life_tables_windows_by_group) || length(life_tables_windows_by_group) == 0) {
    stop("`life_tables_windows_by_group` is empty.", call. = FALSE)
  }

  out <- purrr::map_dfr(names(life_tables_windows_by_group), function(nm) {
    lt <- life_tables_windows_by_group[[nm]]
    if (is.null(lt)) return(NULL)

    lt_groups <- .split_grouped_life_table(lt)
    if (is.null(lt_groups)) return(NULL)

    lt_groups <- lt_groups[names(lt_groups) %in% groups]
    if (length(lt_groups) == 0) return(NULL)

    parts <- strsplit(nm, "_")[[1]]
    window_start <- suppressWarnings(as.numeric(parts[1]))
    window_end   <- suppressWarnings(as.numeric(parts[2]))
    window_mid   <- if (is.finite(window_start) && is.finite(window_end)) {
      (window_start + window_end) / 2
    } else {
      NA_real_
    }

    purrr::map_dfr(names(lt_groups), function(g) {
      lt_tbl <- try(tibble::as_tibble(as.data.frame(lt_groups[[g]])), silent = TRUE)
      if (inherits(lt_tbl, "try-error")) return(NULL)
      if (!measure %in% names(lt_tbl)) return(NULL)

      age_col <- names(lt_tbl)[names(lt_tbl) %in% c("x", "age", "agebeg", "Age", "X")]
      if (length(age_col) == 0) {
        lt_tbl$age_proxy <- seq_len(nrow(lt_tbl))
        age_col <- "age_proxy"
      } else {
        age_col <- age_col[1]
      }

      tibble::tibble(
        window_label = nm,
        window_start = window_start,
        window_end = window_end,
        window_mid = window_mid,
        sex_gender = g,
        age = lt_tbl[[age_col]],
        value = lt_tbl[[measure]]
      )
    })
  })

  if (nrow(out) == 0) {
    stop("No grouped life-table data could be extracted for the requested measure.", call. = FALSE)
  }

  out
}

#' Compute pooled mortAAR metrics across sliding windows
#'
#' @description
#' Computes MI, Ratio_F_M, and MMR2 for the full pooled dataset within each
#' sliding window, ignoring site boundaries.
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#' @param window_width Width of each window in bin units.
#' @param step_size Step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#'
#' @return Tibble with one row per window.
compute_mortaar_pooled_metrics_windows <- function(final,
                                                   prepared,
                                                   window_width,
                                                   step_size,
                                                   start_bin = NULL,
                                                   end_bin = NULL) {
  mort_tbl <- build_mortaar_burial_table(final, prepared)

  windows <- make_mortaar_windows(
    mort_tbl = mort_tbl,
    window_width = window_width,
    step_size = step_size,
    start_bin = start_bin,
    end_bin = end_bin
  )

  purrr::map_dfr(seq_len(nrow(windows)), function(i) {
    w <- windows[i, ]

    mort_tbl_w <- filter_mortaar_by_window(
      mort_tbl,
      window_start = w$window_start,
      window_end = w$window_end
    )

    if (nrow(mort_tbl_w) == 0) {
      return(tibble::tibble(
        window_id = w$window_id,
        window_start = w$window_start,
        window_end = w$window_end,
        window_mid = w$window_mid,
        n_burials = 0,
        n_female = 0,
        n_male = 0,
        n_indet = 0,
        n_female_adult = 0,
        n_male_adult = 0,
        MI = NA_real_,
        Ratio_F_M = NA_real_,
        MMR2 = NA_real_,
        MI_possible = FALSE,
        Ratio_F_M_possible = FALSE,
        MMR2_possible = FALSE,
        any_metric_possible = FALSE,
        note_MI = "No burials in window.",
        note_Ratio_F_M = "No burials in window.",
        note_MMR2 = "No burials in window."
      ))
    }

    compute_site_mortaar_metrics(mort_tbl_w) %>%
      dplyr::mutate(
        window_id = w$window_id,
        window_start = w$window_start,
        window_end = w$window_end,
        window_mid = w$window_mid
      ) %>%
      dplyr::relocate(window_id, window_start, window_end, window_mid, .before = 1)
  })
}

#' Plot dx as a ridge-style temporal series
#'
#' @description
#' Creates a ridge-style plot of the death distribution (`dx`) across sliding
#' windows. Each window is vertically offset on an integer stack and labelled
#' by its window range.
#'
#' @param life_tables_windows Named list returned by `compute_life_tables_windows()`.
#' @param normalize If TRUE, rescale each window's dx to a maximum of 1 before plotting.
#' @param ridge_height Relative ridge height. Values around 0.7-0.9 work well.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_dx_ridge_windows <- function(life_tables_windows,
                                  normalize = TRUE,
                                  ridge_height = 0.8,
                                  save_path = NULL) {
  dat <- extract_life_table_measure_windows(
    life_tables_windows = life_tables_windows,
    measure = "dx"
  )

  dat <- dat %>%
    dplyr::filter(is.finite(value))

  if (nrow(dat) == 0) {
    stop("No finite dx values available for ridge plotting.", call. = FALSE)
  }

  # Robust age parsing
  age_chr <- as.character(dat$age)
  age_num <- suppressWarnings(as.numeric(age_chr))

  need_parse <- !is.finite(age_num)
  age_num[need_parse] <- suppressWarnings(
    as.numeric(sub("-.*", "", age_chr[need_parse]))
  )

  # Final fallback: sequence within each window if age labels cannot be parsed
  if (any(!is.finite(age_num))) {
    dat <- dat %>%
      dplyr::group_by(window_label) %>%
      dplyr::mutate(.age_fallback = dplyr::row_number()) %>%
      dplyr::ungroup()

    age_num[!is.finite(age_num)] <- dat$.age_fallback[!is.finite(age_num)]
  }

  dat <- dat %>%
    dplyr::mutate(age_num = age_num) %>%
    dplyr::arrange(window_start, age_num)

  # Normalize within window if requested
  if (isTRUE(normalize)) {
    dat <- dat %>%
      dplyr::group_by(window_label) %>%
      dplyr::mutate(
        .window_max = max(value, na.rm = TRUE),
        value_plot = dplyr::if_else(
          .window_max > 0,
          value / .window_max,
          0 * value
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.window_max)
  } else {
    dat <- dat %>%
      dplyr::mutate(value_plot = value)
  }

  # Integer stacking index for robust ridge drawing
  window_levels <- dat %>%
    dplyr::distinct(window_label, window_start, window_end, window_mid) %>%
    dplyr::arrange(window_start, window_end) %>%
    dplyr::mutate(window_index = dplyr::row_number())

  dat <- dat %>%
    dplyr::left_join(window_levels, by = c("window_label", "window_start", "window_end", "window_mid")) %>%
    dplyr::mutate(
      y_base = window_index,
      y_top = window_index + value_plot * ridge_height
    )

  axis_breaks <- window_levels$window_index
  axis_labels <- paste0(window_levels$window_start, "-", window_levels$window_end)

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      data = window_levels,
      ggplot2::aes(yintercept = window_index),
      colour = "grey85",
      linewidth = 0.3
    ) +
    ggplot2::geom_ribbon(
      data = dat,
      ggplot2::aes(
        x = age_num,
        ymin = y_base,
        ymax = y_top,
        group = window_label
      ),
      fill = "#2C3E50",
      alpha = 0.55
    ) +
    ggplot2::geom_line(
      data = dat,
      ggplot2::aes(
        x = age_num,
        y = y_top,
        group = window_label
      ),
      colour = "#1B2631",
      linewidth = 0.5
    ) +
    ggplot2::scale_y_continuous(
      breaks = axis_breaks,
      labels = axis_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.05))
    ) +
    ggplot2::labs(
      title = "dx ridge plot across sliding windows",
      x = "Age",
      y = "Window"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  }

  p
}

#' Plot changes in a life-table measure by sex_gender
#'
#' @description
#' Visualizes how one mortAAR life-table measure changes across sliding windows
#' separately for female, male, and indeterminate life tables.
#'
#' @param life_tables_windows_by_sex Named list returned by
#'   `compute_life_tables_windows_by_sex()`.
#' @param measure One of "qx", "dx", "lx", "ex", or "rel_popx".
#' @param mode Either "heatmap" or "line".
#' @param groups Character vector of groups to include.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_life_table_changes_by_sex <- function(life_tables_windows_by_sex,
                                           measure = c("qx", "dx", "lx", "ex", "rel_popx"),
                                           mode = c("heatmap", "line"),
                                           groups = c("f", "m", "indet"),
                                           save_path = NULL) {
  measure <- match.arg(measure)
  mode <- match.arg(mode)

  dat <- extract_life_table_measure_windows_by_group(
    life_tables_windows_by_group = life_tables_windows_by_sex,
    measure = measure,
    groups = groups
  )

  dat <- dat %>%
    dplyr::filter(is.finite(value))

  if (nrow(dat) == 0) {
    stop("No finite grouped values available for plotting.", call. = FALSE)
  }

  dat <- dat %>%
    dplyr::mutate(
      age_lower = suppressWarnings(as.numeric(sub("-.*", "", as.character(age)))),
      sex_gender = factor(sex_gender, levels = c("f", "m", "indet"))
    ) %>%
    dplyr::arrange(sex_gender, window_mid, age_lower)

  age_levels <- unique(dat$age[order(dat$age_lower)])
  dat <- dat %>%
    dplyr::mutate(age = factor(age, levels = age_levels))

  # Show every second age label for readability
  age_breaks <- age_levels[seq(1, length(age_levels), by = 2)]

  y_lab <- dplyr::case_when(
    measure == "qx" ~ "Probability of death (qx)",
    measure == "dx" ~ "Proportion of deaths (dx)",
    measure == "lx" ~ "Survivorship (lx)",
    measure == "ex" ~ "Life expectancy (ex)",
    measure == "rel_popx" ~ "Relative population structure (rel_popx)",
    TRUE ~ measure
  )

  if (mode == "heatmap") {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = age, y = window_mid, fill = value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
      ggplot2::scale_x_discrete(breaks = age_breaks) +
      ggplot2::facet_wrap(~ sex_gender) +
      ggplot2::labs(
        title = paste0("Changes in ", measure, " across sliding windows by sex_gender"),
        x = "Age",
        y = "Window midpoint",
        fill = measure
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
      )
  } else {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = age, y = value, group = window_label, colour = window_mid)
    ) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::scale_colour_viridis_c(option = "C") +
      ggplot2::scale_x_discrete(breaks = age_breaks) +
      ggplot2::facet_wrap(~ sex_gender, scales = "free_y") +
      ggplot2::labs(
        title = paste0("Changes in ", measure, " across sliding windows by sex_gender"),
        x = "Age",
        y = y_lab,
        colour = "Window midpoint"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
      )
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, plot = p, width = 13, height = 7.5, dpi = 300, bg = "white")
  }

  p
}

# -------------------------------------------------------------------------
# Life-table trajectory extraction and plotting
# -------------------------------------------------------------------------

#' Plot changes in a life-table measure across sliding windows
#'
#' @description
#' Visualizes how one mortAAR life-table measure changes across sliding windows.
#' Can be shown either as multiple lines (one per window) or as a heatmap.
#'
#' @param life_tables_windows Named list returned by `compute_life_tables_windows()`.
#' @param measure One of "qx", "dx", "lx", "ex", or "rel_popx".
#' @param mode Either "heatmap" or "line".
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_life_table_changes <- function(life_tables_windows,
                                    measure = c("qx", "dx", "lx", "ex", "rel_popx"),
                                    mode = c("heatmap", "line"),
                                    save_path = NULL) {
  measure <- match.arg(measure)
  mode <- match.arg(mode)

  dat <- extract_life_table_measure_windows(
    life_tables_windows = life_tables_windows,
    measure = measure
  )

  y_lab <- dplyr::case_when(
    measure == "qx" ~ "Probability of death (qx)",
    measure == "dx" ~ "Proportion of deaths (dx)",
    measure == "lx" ~ "Survivorship (lx)",
    measure == "ex" ~ "Life expectancy (ex)",
    measure == "rel_popx" ~ "Relative population structure (rel_popx)",
    TRUE ~ measure
  )

  if (mode == "heatmap") {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = age, y = window_mid, fill = value)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
      ggplot2::labs(
        title = paste0("Changes in ", measure, " across sliding windows"),
        x = "Age",
        y = "Window midpoint",
        fill = measure
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )
  } else {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = age, y = value, group = window_label, colour = window_mid)
    ) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::scale_colour_viridis_c(option = "C") +
      ggplot2::labs(
        title = paste0("Changes in ", measure, " across sliding windows"),
        x = "Age",
        y = y_lab,
        colour = "Window midpoint"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
  }

  p
}

#' Plot all standard life-table change panels
#'
#' @param life_tables_windows Named list returned by `compute_life_tables_windows()`.
#' @param mode Either "heatmap" or "line".
#' @param save_dir Optional directory for saving all plots.
#' @param include_dx_ridge If TRUE, also include a ridge-style dx plot.
#' @param life_tables_windows_by_sex Optional grouped life-table list returned by
#'   `compute_life_tables_windows_by_sex()`.
#' @param include_by_sex If TRUE, also include sex-stratified plots for all measures.
#' @param sex_mode Plot mode for sex-stratified plots, either "heatmap" or "line".
#'
#' @return Named list of ggplot objects.
plot_all_life_table_changes <- function(life_tables_windows,
                                        mode = c("heatmap", "line"),
                                        save_dir = NULL,
                                        include_dx_ridge = FALSE,
                                        life_tables_windows_by_sex = NULL,
                                        include_by_sex = FALSE,
                                        sex_mode = c("heatmap", "line")) {
  mode <- match.arg(mode)
  sex_mode <- match.arg(sex_mode)

  measures <- c("qx", "dx", "lx", "ex", "rel_popx")

  plots <- lapply(measures, function(m) {
    save_path <- NULL
    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      save_path <- file.path(save_dir, paste0("life_table_change_", m, ".png"))
    }

    plot_life_table_changes(
      life_tables_windows = life_tables_windows,
      measure = m,
      mode = mode,
      save_path = save_path
    )
  })

  names(plots) <- measures

  if (isTRUE(include_dx_ridge)) {
    save_path <- NULL
    if (!is.null(save_dir)) {
      save_path <- file.path(save_dir, "life_table_change_dx_ridge.png")
    }

    plots$dx_ridge <- plot_dx_ridge_windows(
      life_tables_windows = life_tables_windows,
      save_path = save_path
    )
  }

  if (isTRUE(include_by_sex)) {
    if (is.null(life_tables_windows_by_sex)) {
      stop("`life_tables_windows_by_sex` must be supplied when `include_by_sex = TRUE`.", call. = FALSE)
    }

    sex_plots <- lapply(measures, function(m) {
      save_path <- NULL
      if (!is.null(save_dir)) {
        save_path <- file.path(save_dir, paste0("life_table_change_", m, "_by_sex.png"))
      }

      plot_life_table_changes_by_sex(
        life_tables_windows_by_sex = life_tables_windows_by_sex,
        measure = m,
        mode = sex_mode,
        save_path = save_path
      )
    })

    names(sex_plots) <- paste0(measures, "_by_sex")
    plots <- c(plots, sex_plots)
  }

  plots
}

.assert_package_mortAAR <- function() {
  if (!requireNamespace("mortAAR", quietly = TRUE)) {
    stop("Package `mortAAR` is required but not installed.", call. = FALSE)
  }
  invisible(TRUE)
}

.run_prep_life_table <- function(x,
                                 agebeg = "from",
                                 ageend = "to",
                                 group = NULL,
                                 method = "Standard",
                                 dec = NA,
                                 agerange = "included") {
  .assert_package_mortAAR()

  if (is.null(group) || length(group) == 0 || identical(group, "none")) {
    mortAAR::prep.life.table(
      x = x,
      dec = dec,
      agebeg = agebeg,
      ageend = ageend,
      method = method,
      agerange = agerange
    )
  } else {
    mortAAR::prep.life.table(
      x = x,
      dec = dec,
      agebeg = agebeg,
      ageend = ageend,
      group = group,
      method = method,
      agerange = agerange
    )
  }
}

.extract_tbl_mortaar <- function(x, path = NULL) {
  if (is.null(x)) return(NULL)

  obj <- x
  if (!is.null(path)) {
    for (nm in path) {
      if (!is.list(obj) || !nm %in% names(obj)) return(NULL)
      obj <- obj[[nm]]
    }
  }

  if (is.data.frame(obj)) {
    return(tibble::as_tibble(obj))
  }

  NULL
}

.snap_to_canonical_bin <- function(x, canonical_bins) {
  canonical_bins <- sort(unique(stats::na.omit(as.numeric(canonical_bins))))

  if (length(canonical_bins) == 0) {
    rlang::abort("`canonical_bins` must contain at least one finite value.")
  }

  vapply(as.numeric(x), function(xx) {
    if (is.na(xx) || !is.finite(xx)) return(NA_real_)
    canonical_bins[which.min(abs(canonical_bins - xx))]
  }, numeric(1))
}


# Backwards-compatible alias. New code should use .snap_to_canonical_bin().
.snap_to_canonical_bucket <- function(x, canonical_buckets) {
  .snap_to_canonical_bin(x, canonical_buckets)
}

# Standardise legacy bin names only at the boundaries of this bridge. The main
# codebase now uses horizon_bin / median_bin, but older saved ensemble objects
# may still contain horizon_bucket / median_bucket.
.mortaar_standardize_bin_columns <- function(df) {
  if (exists(".standardize_bin_columns", mode = "function")) {
    df <- .standardize_bin_columns(df)
  }

  legacy <- c(
    horizon_bucket = "horizon_bin",
    median_bucket = "median_bin",
    bucket_mean = "bin_mean",
    bucket_sd = "bin_sd",
    bucket_q25 = "bin_q25",
    bucket_q75 = "bin_q75",
    bucket_min = "bin_min",
    bucket_max = "bin_max"
  )

  for (old in names(legacy)) {
    new <- legacy[[old]]
    if (old %in% names(df) && !new %in% names(df)) {
      names(df)[names(df) == old] <- new
    }
  }

  df
}

.mortaar_fallback_diagnostics <- function(final, ...) {
  if (exists("flag_uncertain_fallback_bins", mode = "function")) {
    return(flag_uncertain_fallback_bins(final = final, ...))
  }
  if (exists("flag_uncertain_fallback_buckets", mode = "function")) {
    return(flag_uncertain_fallback_buckets(final = final, ...))
  }
  NULL
}

# -------------------------------------------------------------------------
# Prepare burial-level output for mortAAR
# -------------------------------------------------------------------------

#' Prepare burial-level output for mortAAR
#'
#' @description
#' Creates a burial-level table suitable for downstream mortAAR workflows.
#' One row is returned per burial (`UID`). The table preserves as many prepared
#' mortuary and metadata fields as possible, adds a representative temporal
#' assignment defined as the median bin across ensemble replicates, and
#' attaches uncertainty and diagnostic information for that representative bin.
#'
#' Burial-level temporal summaries are computed from replicate allocations and
#' then snapped to the canonical bin grid stored in
#' `final$pooled$canonical_grid`, ensuring consistency with pooled ensemble
#' outputs.
#'
#' This helper requires burial-level replicate allocation data to have been
#' retained during `run_temporal_ensemble()` via
#' `retain_burial_allocations = TRUE`.
#'
#' @param final Finalized ensemble object from `finalize_ensemble()`.
#' @param prepared Optional prepared object from `assemble_prepared_inputs()`.
#' @param include_bin_uncertainty If TRUE, attach pooled uncertainty for the
#'   burial's median bin.
#' @param include_fallback_diagnostic If TRUE, attach fallback/uncertainty
#'   diagnostic fields when available.
#'
#' @return Tibble with one row per burial.
prepare_for_mortAAR <- function(final,
                                prepared = NULL,
                                include_bin_uncertainty = TRUE,
                                include_fallback_diagnostic = TRUE) {
  if (!is.list(final) ||
      is.null(final$replicate_data) ||
      is.null(final$replicate_data$burial_allocations)) {
    rlang::abort(
      paste0(
        "`final` does not contain burial-level replicate allocations. ",
        "Re-run `run_temporal_ensemble(..., retain_burial_allocations = TRUE)`."
      )
    )
  }

  alloc_reps <- final$replicate_data$burial_allocations
  alloc_reps <- Filter(Negate(is.null), alloc_reps)

  if (length(alloc_reps) == 0) {
    rlang::abort("No burial-level replicate allocations available.")
  }

  canonical_grid <- NULL
  if (!is.null(final$pooled) && !is.null(final$pooled$canonical_grid)) {
    canonical_grid <- tibble::as_tibble(final$pooled$canonical_grid) %>%
      .mortaar_standardize_bin_columns()
  }

  if (is.null(canonical_grid) || !"horizon_bin" %in% names(canonical_grid)) {
    rlang::abort("`final$pooled$canonical_grid` is required for mortAAR preparation.")
  }

  canonical_bins <- sort(unique(canonical_grid$horizon_bin))

  alloc_long <- dplyr::bind_rows(
    lapply(seq_along(alloc_reps), function(i) {
      x <- tibble::as_tibble(alloc_reps[[i]]) %>%
        .mortaar_standardize_bin_columns()
      x$replicate_id <- i
      x
    })
  )

  validate_required_fields(
    alloc_long,
    c("UID", "horizon_bin"),
    "burial_allocations"
  )

  raw_q25 <- function(x) stats::quantile(x, probs = 0.25, na.rm = TRUE, names = FALSE)
  raw_q75 <- function(x) stats::quantile(x, probs = 0.75, na.rm = TRUE, names = FALSE)

  burial_summary <- alloc_long %>%
    dplyr::group_by(UID) %>%
    dplyr::summarise(
      median_bin_raw = stats::median(as.numeric(horizon_bin), na.rm = TRUE),
      bin_mean_raw = mean(as.numeric(horizon_bin), na.rm = TRUE),
      bin_sd = stats::sd(as.numeric(horizon_bin), na.rm = TRUE),
      bin_q25_raw = raw_q25(as.numeric(horizon_bin)),
      bin_q75_raw = raw_q75(as.numeric(horizon_bin)),
      bin_min_raw = min(as.numeric(horizon_bin), na.rm = TRUE),
      bin_max_raw = max(as.numeric(horizon_bin), na.rm = TRUE),
      bin_n_replicates = dplyr::n_distinct(replicate_id),
      chronology_source = dplyr::first(chronology_source),
      is_synthetic = dplyr::first(is_synthetic),
      input_type = dplyr::first(input_type),
      site_id = dplyr::first(site_id),
      burial_id = dplyr::first(burial_id),
      record_id = dplyr::first(record_id),
      phase_id = dplyr::first(phase_id),
      system_name = dplyr::first(system_name),
      phase_name = dplyr::first(phase_name),
      sex_gender = dplyr::first(sex_gender),
      age = dplyr::first(age),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      median_bin = .snap_to_canonical_bin(median_bin_raw, canonical_bins),
      bin_q25 = .snap_to_canonical_bin(bin_q25_raw, canonical_bins),
      bin_q75 = .snap_to_canonical_bin(bin_q75_raw, canonical_bins),
      bin_min = .snap_to_canonical_bin(bin_min_raw, canonical_bins),
      bin_max = .snap_to_canonical_bin(bin_max_raw, canonical_bins),
      bin_mean = bin_mean_raw
    )

  out <- burial_summary

  mortuary_tbl <- .extract_tbl_mortaar(prepared, c("data", "mortuary"))
  site_meta_tbl <- .extract_tbl_mortaar(prepared, c("data", "site_metadata"))

  if (!is.null(mortuary_tbl) && "UID" %in% names(mortuary_tbl)) {
    keep_mort_cols <- setdiff(names(mortuary_tbl), names(out))
    out <- dplyr::left_join(
      out,
      mortuary_tbl[, c("UID", keep_mort_cols), drop = FALSE],
      by = "UID"
    )
  }

  if (!is.null(site_meta_tbl) && "site_id" %in% names(site_meta_tbl)) {
    keep_site_cols <- setdiff(names(site_meta_tbl), names(out))
    out <- dplyr::left_join(
      out,
      site_meta_tbl[, c("site_id", keep_site_cols), drop = FALSE],
      by = "site_id"
    )
  }

  if (isTRUE(include_bin_uncertainty)) {
    pooled_tbl <- .extract_tbl_mortaar(final, c("pooled", "estimates"))

    if (!is.null(pooled_tbl)) {
      pooled_tbl <- .mortaar_standardize_bin_columns(pooled_tbl)
      bin_unc <- pooled_tbl %>%
        dplyr::select(horizon_bin, W, B, T, se_total, lower, upper) %>%
        dplyr::rename(
          bin_W = W,
          bin_B = B,
          bin_T = T,
          bin_se_total = se_total,
          bin_lower = lower,
          bin_upper = upper
        )

      names(bin_unc)[1] <- "median_bin"

      out <- dplyr::left_join(out, bin_unc, by = "median_bin")
    }
  }

  if (isTRUE(include_fallback_diagnostic)) {
    diag_tbl <- try(
      .mortaar_fallback_diagnostics(
        final = final,
        provenance_source = final,
        fallback_threshold = 0.25,
        ratio_threshold = 1,
        normalized = TRUE
      ),
      silent = TRUE
    )

    if (!inherits(diag_tbl, "try-error") && !is.null(diag_tbl)) {
      diag_tbl <- tibble::as_tibble(diag_tbl) %>%
        .mortaar_standardize_bin_columns() %>%
        dplyr::select(horizon_bin, fallback_share, BW_ratio, BT_share, flagged) %>%
        dplyr::rename(
          median_bin = horizon_bin,
          fallback_share_at_bin = fallback_share,
          BW_ratio_at_bin = BW_ratio,
          BT_share_at_bin = BT_share,
          flagged_fallback_uncertainty = flagged
        )

      out <- dplyr::left_join(out, diag_tbl, by = "median_bin")
    }
  }

  tibble::as_tibble(out)
}

# -------------------------------------------------------------------------
# Age parsing
# -------------------------------------------------------------------------

#' Parse archaeological age labels into numeric lower/upper bounds
#'
#' @description
#' Converts a mixed age field into numeric `from` / `to` bounds suitable for
#' mortAAR. Explicit numeric ranges such as `"20-30"` are parsed directly.
#' Ontology labels such as `"inf1"` or `"juv"` are mapped onto standard ranges.
#' Single numeric ages are placed into the same standard brackets.
#'
#' @param age_vec Character or numeric vector of age labels.
#'
#' @return Tibble with columns `from` and `to`.
parse_mortaar_age_bounds <- function(age_vec) {
  age_chr <- stringr::str_trim(as.character(age_vec))

  is_range <- grepl("^\\d+\\s*-\\s*\\d+$", age_chr)

  from <- rep(NA_real_, length(age_chr))
  to   <- rep(NA_real_, length(age_chr))

  if (any(is_range)) {
    parts <- strsplit(gsub("\\s+", "", age_chr[is_range]), "-")
    from[is_range] <- as.numeric(vapply(parts, `[`, 1, FUN.VALUE = character(1)))
    to[is_range]   <- as.numeric(vapply(parts, `[`, 2, FUN.VALUE = character(1)))
  }

  idx <- !is_range

  mapped <- dplyr::case_when(
    age_chr[idx] == "inf1" ~ "0-6",
    age_chr[idx] == "inf2" ~ "7-14",
    age_chr[idx] == "juv"  ~ "15-20",
    age_chr[idx] == "fad"  ~ "21-29",
    age_chr[idx] == "ad"   ~ "30-40",
    age_chr[idx] == "mat"  ~ "41-60",
    age_chr[idx] == "sen"  ~ "61-85",
    TRUE ~ age_chr[idx]
  )

  mapped_num <- suppressWarnings(as.numeric(mapped))
  mapped <- dplyr::case_when(
    !is.na(mapped_num) & mapped_num >= 0  & mapped_num <= 6  ~ "0-6",
    !is.na(mapped_num) & mapped_num >= 7  & mapped_num <= 14 ~ "7-14",
    !is.na(mapped_num) & mapped_num >= 15 & mapped_num <= 20 ~ "15-20",
    !is.na(mapped_num) & mapped_num >= 21 & mapped_num <= 29 ~ "21-29",
    !is.na(mapped_num) & mapped_num >= 30 & mapped_num <= 40 ~ "30-40",
    !is.na(mapped_num) & mapped_num >= 41 & mapped_num <= 60 ~ "41-60",
    !is.na(mapped_num) & mapped_num >= 61 & mapped_num <= 85 ~ "61-85",
    TRUE ~ mapped
  )

  split2 <- strsplit(gsub("\\s+", "", mapped), "-")
  good <- lengths(split2) == 2

  from[idx][good] <- as.numeric(vapply(split2[good], `[`, 1, FUN.VALUE = character(1)))
  to[idx][good]   <- as.numeric(vapply(split2[good], `[`, 2, FUN.VALUE = character(1)))

  tibble::tibble(from = from, to = to)
}

# -------------------------------------------------------------------------
# Sex normalization
# -------------------------------------------------------------------------

#' Normalize sex labels for mortAAR
#'
#' @param x Character vector of sex labels.
#'
#' @return Character vector with values "f", "m", "indet", or NA.
normalize_mortaar_sex <- function(x) {
  x_chr <- stringr::str_trim(stringr::str_to_lower(as.character(x)))

  dplyr::case_when(
    x_chr %in% c("f", "female", "frau", "w", "weiblich") ~ "f",
    x_chr %in% c("m", "male", "mann", "maennlich", "männlich") ~ "m",
    x_chr %in% c("indet", "indeterminate", "indeterminable", "undet", "unknown", "u") ~ "indet",
    TRUE ~ NA_character_
  )
}

# -------------------------------------------------------------------------
# Build mortAAR-ready burial table
# -------------------------------------------------------------------------

#' Build a mortAAR-ready burial table
#'
#' @description
#' Takes the burial-level output from `prepare_for_mortAAR()` and converts it
#' into a mortAAR-ready table by adding numeric age bounds and normalized sex.
#' Indeterminate-sex burials are retained.
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#'
#' @return Tibble with burial-level rows and numeric `from` / `to` age bounds.
build_mortaar_burial_table <- function(final, prepared) {
  mortAAR_input <- prepare_for_mortAAR(
    final = final,
    prepared = prepared,
    include_bin_uncertainty = TRUE,
    include_fallback_diagnostic = TRUE
  )

  age_bounds <- parse_mortaar_age_bounds(mortAAR_input$age)

  mortAAR_input %>%
    dplyr::bind_cols(age_bounds) %>%
    dplyr::mutate(
      sex_gender = normalize_mortaar_sex(sex_gender)
    ) %>%
    dplyr::filter(
      !is.na(site_name),
      !is.na(from),
      !is.na(to),
      !is.na(sex_gender)
    ) %>%
    tibble::as_tibble()
}

# -------------------------------------------------------------------------
# Extract mortAAR lt.sexrelation() output
# -------------------------------------------------------------------------

#' Extract MI, female:male ratio, MMR1, and MMR2 from lt.sexrelation()
#'
#' @param rel_obj Result from `mortAAR::lt.sexrelation()`.
#'
#' @return Tibble with columns `MI`, `Ratio_F_M`, `MMR1`, and `MMR2`.
extract_sexrelation_metrics <- function(rel_obj) {
  rel_df <- tibble::as_tibble(as.data.frame(rel_obj))

  # Standard mortAAR output uses columns: method, value, description
  if (all(c("method", "value") %in% names(rel_df))) {
    rel_df <- rel_df %>%
      dplyr::mutate(
        method = as.character(method),
        value = suppressWarnings(as.numeric(value))
      )

    mi_val <- rel_df %>%
      dplyr::filter(method == "MI") %>%
      dplyr::pull(value)
    mi_val <- if (length(mi_val) > 0) mi_val[[1]] else NA_real_

    ratio_fm_val <- rel_df %>%
      dplyr::filter(method == "Ratio_F_M") %>%
      dplyr::pull(value)
    ratio_fm_val <- if (length(ratio_fm_val) > 0) ratio_fm_val[[1]] else NA_real_

    mmr1_val <- rel_df %>%
      dplyr::filter(method == "MMR1") %>%
      dplyr::pull(value)
    mmr1_val <- if (length(mmr1_val) > 0) mmr1_val[[1]] else NA_real_

    mmr2_val <- rel_df %>%
      dplyr::filter(method == "MMR2") %>%
      dplyr::pull(value)
    mmr2_val <- if (length(mmr2_val) > 0) mmr2_val[[1]] else NA_real_

    return(tibble::tibble(
      MI = mi_val,
      Ratio_F_M = ratio_fm_val,
      MMR1 = mmr1_val,
      MMR2 = mmr2_val
    ))
  }

  # Fallback for differently structured outputs
  rel_df <- rel_df %>%
    tibble::rownames_to_column("metric_raw") %>%
    dplyr::mutate(metric = stringr::str_to_lower(metric_raw))

  num_cols <- names(rel_df)[vapply(rel_df, is.numeric, logical(1))]
  if (length(num_cols) == 0) {
    return(tibble::tibble(MI = NA_real_, Ratio_F_M = NA_real_, MMR1 = NA_real_, MMR2 = NA_real_))
  }

  value_col <- num_cols[1]
  rel_df$value <- rel_df[[value_col]]

  mi_val <- rel_df %>%
    dplyr::filter(stringr::str_detect(metric, "mi|mascul")) %>%
    dplyr::pull(value)
  mi_val <- if (length(mi_val) > 0) mi_val[[1]] else NA_real_

  ratio_fm_val <- rel_df %>%
    dplyr::filter(stringr::str_detect(metric, "ratio") &
                    stringr::str_detect(metric, "f") &
                    stringr::str_detect(metric, "m")) %>%
    dplyr::pull(value)
  ratio_fm_val <- if (length(ratio_fm_val) > 0) ratio_fm_val[[1]] else NA_real_

  mmr1_val <- rel_df %>%
    dplyr::filter(stringr::str_detect(metric, "mmr1")) %>%
    dplyr::pull(value)
  mmr1_val <- if (length(mmr1_val) > 0) mmr1_val[[1]] else NA_real_

  mmr2_val <- rel_df %>%
    dplyr::filter(stringr::str_detect(metric, "mmr2|maternal")) %>%
    dplyr::pull(value)
  mmr2_val <- if (length(mmr2_val) > 0) mmr2_val[[1]] else NA_real_

  tibble::tibble(
    MI = mi_val,
    Ratio_F_M = ratio_fm_val,
    MMR1 = mmr1_val,
    MMR2 = mmr2_val
  )
}

# -------------------------------------------------------------------------
# Site-level mortAAR metrics
# -------------------------------------------------------------------------

#' Compute site-level mortAAR sex-relation metrics
#'
#' @description
#' Builds female and male life tables for one site and extracts:
#' - Masculinity Index (MI)
#' - Ratio of Females to Males (Ratio_F_M)
#' - Motherhood Mortality Rate estimator 2 (MMR2)
#'
#' Indeterminate-sex burials are retained in the input table but excluded from
#' the sex-relation calculations themselves.
#'
#' @param df_site MortAAR-ready burial table for a single site.
#'
#' @return Tibble with site-level metrics, per-metric feasibility flags, and notes.
compute_site_mortaar_metrics <- function(df_site) {
  .assert_package_mortAAR()

  n_burials <- nrow(df_site)
  n_female <- sum(df_site$sex_gender == "f", na.rm = TRUE)
  n_male   <- sum(df_site$sex_gender == "m", na.rm = TRUE)
  n_indet  <- sum(df_site$sex_gender == "indet", na.rm = TRUE)

  n_female_adult <- sum(df_site$sex_gender == "f" & !is.na(df_site$from) & df_site$from >= 15, na.rm = TRUE)
  n_male_adult   <- sum(df_site$sex_gender == "m" & !is.na(df_site$from) & df_site$from >= 15, na.rm = TRUE)

  MI <- NA_real_
  Ratio_F_M <- NA_real_
  MMR1 <- NA_real_
  MMR2 <- NA_real_

  MI_possible <- FALSE
  Ratio_F_M_possible <- FALSE
  MMR1_possible <- FALSE
  MMR2_possible <- FALSE

  note_MI <- NA_character_
  note_Ratio_F_M <- NA_character_
  note_MMR1 <- NA_character_
  note_MMR2 <- NA_character_

  df_site_sexed <- df_site %>%
    dplyr::filter(sex_gender %in% c("f", "m"))

  rel_obj <- NULL

  if (all(c("f", "m") %in% unique(df_site_sexed$sex_gender))) {
    prep_obj <- try(
      .run_prep_life_table(
        x = df_site_sexed,
        group = "sex_gender",
        agebeg = "from",
        ageend = "to",
        agerange = "included"
      ),
      silent = TRUE
    )

    if (!inherits(prep_obj, "try-error")) {
      life_obj <- try(mortAAR::life.table(prep_obj), silent = TRUE)

      if (!inherits(life_obj, "try-error") && !is.null(life_obj$f) && !is.null(life_obj$m)) {
        rel_obj <- try(mortAAR::lt.sexrelation(life_obj$f, life_obj$m), silent = TRUE)

        if (!inherits(rel_obj, "try-error")) {
          metrics <- extract_sexrelation_metrics(rel_obj)

          if (!is.na(metrics$MI) && is.finite(metrics$MI)) {
            MI <- metrics$MI
            MI_possible <- TRUE
          }

          if (!is.na(metrics$Ratio_F_M) && is.finite(metrics$Ratio_F_M)) {
            Ratio_F_M <- metrics$Ratio_F_M
            Ratio_F_M_possible <- TRUE
          }

          if (!is.na(metrics$MMR1) && is.finite(metrics$MMR1) && metrics$MMR1 >= 0) {
            MMR1 <- metrics$MMR1
            MMR1_possible <- TRUE
          }

          if (!is.na(metrics$MMR2) && is.finite(metrics$MMR2) && metrics$MMR2 >= 0) {
            MMR2 <- metrics$MMR2
            MMR2_possible <- TRUE
          }
        }
      }
    }
  }

  if (!MI_possible) {
    if (n_female_adult > 0 && n_male_adult > 0) {
      MI <- n_male_adult / n_female_adult
      MI_possible <- TRUE
      note_MI <- "Computed from adult (15+) sexed counts because mortAAR did not return MI."
    } else {
      note_MI <- "Insufficient adult (15+) female and male counts for MI."
    }
  }

  if (!Ratio_F_M_possible) {
    if (n_female_adult > 0 && n_male_adult > 0) {
      Ratio_F_M <- n_female_adult / n_male_adult
      Ratio_F_M_possible <- TRUE
      note_Ratio_F_M <- "Computed from adult (15+) sexed counts because mortAAR did not return the ratio."
    } else {
      note_Ratio_F_M <- "Insufficient adult (15+) female and male counts for Ratio_F_M."
    }
  }

  if (!MMR1_possible) {
    note_MMR1 <- "MMR1 could not be computed from mortAAR output for this site."
  }

  if (!MMR2_possible) {
    note_MMR2 <- "MMR2 could not be computed from mortAAR output for this site."
  }

  any_metric_possible <- MI_possible | Ratio_F_M_possible | MMR1_possible | MMR2_possible

  tibble::tibble(
    n_burials = n_burials,
    n_female = n_female,
    n_male = n_male,
    n_indet = n_indet,
    n_female_adult = n_female_adult,
    n_male_adult = n_male_adult,
    MI = MI,
    Ratio_F_M = Ratio_F_M,
    MMR1 = MMR1,
    MMR2 = MMR2,
    MI_possible = MI_possible,
    Ratio_F_M_possible = Ratio_F_M_possible,
    MMR1_possible = MMR1_possible,
    MMR2_possible = MMR2_possible,
    any_metric_possible = any_metric_possible,
    note_MI = note_MI,
    note_Ratio_F_M = note_Ratio_F_M,
    note_MMR1 = note_MMR1,
    note_MMR2 = note_MMR2
  )
}

# -------------------------------------------------------------------------
# Filter helper
# -------------------------------------------------------------------------

#' Keep only rows where at least one mortAAR metric is available
#'
#' @param metrics_tbl Output from site-level mortAAR metric helpers.
#'
#' @return Filtered tibble.
filter_mortaar_any_possible <- function(metrics_tbl) {
  metrics_tbl %>%
    dplyr::filter(
      dplyr::if_any(
        dplyr::all_of(c("MI_possible", "Ratio_F_M_possible", "MMR1_possible", "MMR2_possible")),
        ~ .x
      )
    )
}

# -------------------------------------------------------------------------
# Whole-dataset site metrics
# -------------------------------------------------------------------------

#' Compute mortAAR site metrics for the full dataset
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#'
#' @return Tibble with one row per site.
compute_mortaar_site_metrics_all <- function(final, prepared) {
  mort_tbl <- build_mortaar_burial_table(final, prepared)

  mort_tbl %>%
    dplyr::group_by(site_id, site_name) %>%
    dplyr::group_modify(~ compute_site_mortaar_metrics(.x)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(site_name)
}

# -------------------------------------------------------------------------
# Sliding windows
# -------------------------------------------------------------------------

#' Create sliding time windows on the median-bin scale
#'
#' @param mort_tbl MortAAR-ready burial table.
#' @param window_width Width of each window in bin units.
#' @param step_size Step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#'
#' @return Tibble with `window_id`, `window_start`, `window_end`, `window_mid`.
make_mortaar_windows <- function(mort_tbl,
                                 window_width,
                                 step_size,
                                 start_bin = NULL,
                                 end_bin = NULL) {
  if (!is.numeric(window_width) || length(window_width) != 1L || is.na(window_width) || window_width <= 0) {
    stop("`window_width` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(step_size) || length(step_size) != 1L || is.na(step_size) || step_size <= 0) {
    stop("`step_size` must be a single positive number.", call. = FALSE)
  }

  mb <- sort(unique(stats::na.omit(mort_tbl$median_bin)))
  if (length(mb) == 0) {
    stop("No valid `median_bin` values available.", call. = FALSE)
  }

  bin_step <- if (length(mb) > 1) min(diff(mb)) else step_size

  if (is.null(start_bin)) {
    start_bin <- floor(min(mb, na.rm = TRUE) / bin_step) * bin_step
  }
  if (is.null(end_bin)) {
    end_bin <- ceiling(max(mb, na.rm = TRUE) / bin_step) * bin_step
  }

  starts <- seq(from = start_bin, to = end_bin - window_width, by = step_size)
  if (length(starts) == 0) {
    starts <- start_bin
  }

  tibble::tibble(
    window_id = seq_along(starts),
    window_start = starts,
    window_end = starts + window_width,
    window_mid = starts + window_width / 2
  )
}

#' Filter mortAAR burial table by one sliding window
#'
#' @param mort_tbl MortAAR-ready burial table.
#' @param window_start Lower inclusive bound.
#' @param window_end Upper inclusive bound.
#'
#' @return Filtered tibble.
filter_mortaar_by_window <- function(mort_tbl, window_start, window_end) {
  mort_tbl %>%
    dplyr::filter(
      !is.na(median_bin),
      median_bin >= window_start,
      median_bin <= window_end
    )
}

# -------------------------------------------------------------------------
# Site metrics by sliding windows
# -------------------------------------------------------------------------

#' Compute mortAAR site metrics across sliding windows
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#' @param window_width Width of each window in bin units.
#' @param step_size Step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#'
#' @return Tibble with one row per site per window.
compute_mortaar_site_metrics_windows <- function(final,
                                                 prepared,
                                                 window_width,
                                                 step_size,
                                                 start_bin = NULL,
                                                 end_bin = NULL) {
  mort_tbl <- build_mortaar_burial_table(final, prepared)
  windows <- make_mortaar_windows(
    mort_tbl = mort_tbl,
    window_width = window_width,
    step_size = step_size,
    start_bin = start_bin,
    end_bin = end_bin
  )

  purrr::map_dfr(seq_len(nrow(windows)), function(i) {
    w <- windows[i, ]

    mort_tbl_w <- filter_mortaar_by_window(
      mort_tbl,
      window_start = w$window_start,
      window_end = w$window_end
    )

    if (nrow(mort_tbl_w) == 0) {
      return(tibble::tibble(
        window_id = w$window_id,
        window_start = w$window_start,
        window_end = w$window_end,
        window_mid = w$window_mid,
        site_id = NA_character_,
        site_name = NA_character_,
        n_burials = 0,
        n_female = 0,
        n_male = 0,
        n_indet = 0,
        n_female_adult = 0,
        n_male_adult = 0,
        MI = NA_real_,
        Ratio_F_M = NA_real_,
        MMR2 = NA_real_,
        MI_possible = FALSE,
        Ratio_F_M_possible = FALSE,
        MMR2_possible = FALSE,
        any_metric_possible = FALSE,
        note_MI = "No burials in window.",
        note_Ratio_F_M = "No burials in window.",
        note_MMR2 = "No burials in window."
      ))
    }

    mort_tbl_w %>%
      dplyr::group_by(site_id, site_name) %>%
      dplyr::group_modify(~ compute_site_mortaar_metrics(.x)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        window_id = w$window_id,
        window_start = w$window_start,
        window_end = w$window_end,
        window_mid = w$window_mid
      ) %>%
      dplyr::relocate(window_id, window_start, window_end, window_mid, .before = 1)
  })
}

# -------------------------------------------------------------------------
# Life tables by sliding windows
# -------------------------------------------------------------------------

#' Compute mortAAR life tables across sliding windows
#'
#' @param final Finalized ensemble object.
#' @param prepared Prepared input object.
#' @param window_width Width of each window in bin units.
#' @param step_size Step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#' @param group Optional grouping column for mortAAR, e.g. "site_name".
#'   Use NULL or "none" for the whole dataset.
#' @param method Method passed to `mortAAR::prep.life.table()`.
#'
#' @return Named list of mortAAR life-table objects.
compute_life_tables_windows <- function(final,
                                        prepared,
                                        window_width,
                                        step_size,
                                        start_bin = NULL,
                                        end_bin = NULL,
                                        group = NULL,
                                        method = "Standard") {
  .assert_package_mortAAR()

  mort_tbl <- build_mortaar_burial_table(final, prepared)
  windows <- make_mortaar_windows(
    mort_tbl = mort_tbl,
    window_width = window_width,
    step_size = step_size,
    start_bin = start_bin,
    end_bin = end_bin
  )

  out <- purrr::map(seq_len(nrow(windows)), function(i) {
    w <- windows[i, ]

    mort_tbl_w <- filter_mortaar_by_window(
      mort_tbl,
      window_start = w$window_start,
      window_end = w$window_end
    )

    if (nrow(mort_tbl_w) == 0) {
      return(NULL)
    }

    prep_obj <- .run_prep_life_table(
      x = mort_tbl_w,
      dec = NA,
      agebeg = "from",
      ageend = "to",
      group = group,
      method = method,
      agerange = "included"
    )

    mortAAR::life.table(prep_obj)
  })

  names(out) <- paste0(windows$window_start, "_", windows$window_end)
  out
}

# -------------------------------------------------------------------------
# Main wrapper
# -------------------------------------------------------------------------

#' Run the mortAAR bridge workflow
#'
#' @description
#' Convenience wrapper that prepares burial-level mortAAR input from the
#' ensemble output and computes site-level sex-relation metrics for the full
#' dataset and, optionally, across sliding time windows.
#'
#' @param final Finalized ensemble object from `finalize_ensemble()`.
#' @param prepared Prepared input object from `assemble_prepared_inputs()`.
#' @param window_width Optional width of sliding windows in bin units.
#' @param step_size Optional step size between successive windows.
#' @param start_bin Optional lower bound for the first window start.
#' @param end_bin Optional upper bound for the final window end.
#' @param return_life_tables If TRUE, also compute mortAAR life tables for each window.
#' @param life_table_group Optional grouping variable passed to
#'   `prep.life.table()`, e.g. "site_name". Use NULL or "none" for the whole dataset.
#' @param life_table_method Method passed to mortAAR `prep.life.table()`.
#'
#' @return A list with:
#' \itemize{
#'   \item `mortaar_burial_table`
#'   \item `site_metrics_all`
#'   \item `site_metrics_all_any_possible`
#'   \item `site_metrics_windows`
#'   \item `site_metrics_windows_any_possible`
#'   \item `life_tables_windows`
#'   \item `windows`
#' }
run_mortAAR_bridge <- function(final,
                               prepared,
                               window_width = NULL,
                               step_size = NULL,
                               start_bin = NULL,
                               end_bin = NULL,
                               return_life_tables = FALSE,
                               life_table_group = NULL,
                               life_table_method = "Standard") {
  if (!is.list(final)) {
    rlang::abort("`final` must be a finalized ensemble object.")
  }

  if (identical(life_table_group, "none")) {
    life_table_group <- NULL
  }

  mortaar_burial_table <- build_mortaar_burial_table(final, prepared)
  site_metrics_all <- compute_mortaar_site_metrics_all(final, prepared)

  windows <- NULL
  site_metrics_windows <- NULL
  pooled_metrics_windows <- NULL
  life_tables_windows <- NULL

  if (!is.null(window_width) && !is.null(step_size)) {
    windows <- make_mortaar_windows(
      mort_tbl = mortaar_burial_table,
      window_width = window_width,
      step_size = step_size,
      start_bin = start_bin,
      end_bin = end_bin
    )

    site_metrics_windows <- compute_mortaar_site_metrics_windows(
      final = final,
      prepared = prepared,
      window_width = window_width,
      step_size = step_size,
      start_bin = start_bin,
      end_bin = end_bin
    )

    pooled_metrics_windows <- compute_mortaar_pooled_metrics_windows(
      final = final,
      prepared = prepared,
      window_width = window_width,
      step_size = step_size,
      start_bin = start_bin,
      end_bin = end_bin
    )

    if (isTRUE(return_life_tables)) {
      life_tables_windows <- compute_life_tables_windows(
        final = final,
        prepared = prepared,
        window_width = window_width,
        step_size = step_size,
        start_bin = start_bin,
        end_bin = end_bin,
        group = life_table_group,
        method = life_table_method
      )
    }
  }

  list(
    mortaar_burial_table = mortaar_burial_table,
    site_metrics_all = site_metrics_all,
    site_metrics_all_any_possible = filter_mortaar_any_possible(site_metrics_all),
    site_metrics_windows = site_metrics_windows,
    site_metrics_windows_any_possible = if (is.null(site_metrics_windows)) NULL else filter_mortaar_any_possible(site_metrics_windows),
    pooled_metrics_windows = pooled_metrics_windows,
    pooled_metrics_windows_any_possible = if (is.null(pooled_metrics_windows)) NULL else filter_mortaar_any_possible(pooled_metrics_windows),
    life_tables_windows = life_tables_windows,
    windows = windows
  )
}

# -------------------------------------------------------------------------
# Optional export helper
# -------------------------------------------------------------------------

#' Export mortAAR bridge outputs
#'
#' @param mortaar_run Result from `run_mortAAR_bridge()`.
#' @param out_dir Output directory.
#'
#' @return Output directory invisibly.
export_mortAAR_bridge <- function(mortaar_run, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    mortaar_run$mortaar_burial_table,
    file.path(out_dir, "mortaar_burial_table.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    mortaar_run$site_metrics_all,
    file.path(out_dir, "mortaar_site_metrics_all.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    mortaar_run$site_metrics_all_any_possible,
    file.path(out_dir, "mortaar_site_metrics_all_any_possible.csv"),
    row.names = FALSE
  )

  if (!is.null(mortaar_run$site_metrics_windows)) {
    utils::write.csv(
      mortaar_run$site_metrics_windows,
      file.path(out_dir, "mortaar_site_metrics_windows.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(mortaar_run$site_metrics_windows_any_possible)) {
    utils::write.csv(
      mortaar_run$site_metrics_windows_any_possible,
      file.path(out_dir, "mortaar_site_metrics_windows_any_possible.csv"),
      row.names = FALSE
    )
  }

  invisible(out_dir)
}

# -------------------------------------------------------------------------
# Sliding-window trajectory plotting
# -------------------------------------------------------------------------

#' Plot one mortAAR metric through sliding windows
#'
#' @description
#' Plots a selected mortAAR metric across sliding temporal windows using the
#' window midpoint as the x-axis.
#'
#' @param metrics_windows Output from `compute_mortaar_site_metrics_windows()`
#'   or `run_mortAAR_bridge()$site_metrics_windows`.
#' @param metric One of `"MI"`, `"Ratio_F_M"`, or `"MMR2"`.
#' @param facet_by_site If TRUE, facet by site. If FALSE, draw all selected
#'   sites in one panel.
#' @param site_names Optional character vector of site names to include.
#' @param only_possible If TRUE, keep only rows where the selected metric is
#'   actually available.
#' @param add_points If TRUE, add points to the trajectory lines.
#' @param add_overall If TRUE, add a pooled mean trajectory across included
#'   sites. This is descriptive only.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_mortaar_metric_trajectory <- function(metrics_windows,
                                           metric = c("MI", "Ratio_F_M", "MMR1", "MMR2"),
                                           facet_by_site = TRUE,
                                           site_names = NULL,
                                           only_possible = TRUE,
                                           add_points = TRUE,
                                           add_overall = FALSE,
                                           aggregate_sites = FALSE,
                                           save_path = NULL) {
  metric <- match.arg(metric)
  dat <- tibble::as_tibble(metrics_windows)

  validate_required_fields(
    dat,
    c("window_mid", metric),
    "metrics_windows"
  )

  if (!"site_name" %in% names(dat)) {
    aggregate_sites <- TRUE
  }

  possible_col <- paste0(metric, "_possible")
  if (isTRUE(only_possible) && possible_col %in% names(dat)) {
    dat <- dat %>%
      dplyr::filter(.data[[possible_col]])
  }

  if (!aggregate_sites) {
    validate_required_fields(dat, "site_name", "metrics_windows")

    if (!is.null(site_names)) {
      dat <- dat %>%
        dplyr::filter(.data$site_name %in% site_names)
    }
  }

  # General finite-value filter
  dat <- dat %>%
    dplyr::filter(is.finite(.data[[metric]]))

  # Additional strict filter for maternal mortality measures
  if (metric %in% c("MMR1", "MMR2")) {
    dat <- dat %>%
      dplyr::filter(!is.na(.data[[metric]]), is.finite(.data[[metric]]), .data[[metric]] >= 0)
  }

  if (nrow(dat) == 0) {
    stop("No rows available for plotting after filtering.", call. = FALSE)
  }

  y_lab <- dplyr::case_when(
    metric == "MI" ~ "Masculinity Index (MI)",
    metric == "Ratio_F_M" ~ "Ratio of Females to Males",
    metric == "MMR1" ~ "Maternal Mortality Rate (MMR1)",
    metric == "MMR2" ~ "Maternal Mortality Rate (MMR2)",
    TRUE ~ metric
  )

  if (aggregate_sites) {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = window_mid, y = .data[[metric]])
    ) +
      ggplot2::geom_line(linewidth = 0.9, colour = "#2C3E50") +
      ggplot2::labs(
        x = "Window midpoint",
        y = y_lab,
        title = paste0(metric, " through sliding windows (all sites pooled)")
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )

    if (isTRUE(add_points)) {
      p <- p + ggplot2::geom_point(size = 1.8, colour = "#2C3E50")
    }
  } else if (isTRUE(facet_by_site)) {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = window_mid, y = .data[[metric]], group = site_name)
    ) +
      ggplot2::geom_line(linewidth = 0.7, colour = "#2C3E50") +
      ggplot2::labs(
        x = "Window midpoint",
        y = y_lab,
        title = paste0(metric, " through sliding windows")
      ) +
      ggplot2::facet_wrap(~ site_name, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold")
      )

    if (isTRUE(add_points)) {
      p <- p + ggplot2::geom_point(size = 1.6, colour = "#2C3E50")
    }

    if (isTRUE(add_overall)) {
      overall <- dat %>%
        dplyr::group_by(window_mid) %>%
        dplyr::summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

      p <- p +
        ggplot2::geom_line(
          data = overall,
          ggplot2::aes(x = window_mid, y = value),
          inherit.aes = FALSE,
          linewidth = 1.1,
          colour = "#C0392B"
        )
    }
  } else {
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = window_mid, y = .data[[metric]], colour = site_name, group = site_name)
    ) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::labs(
        x = "Window midpoint",
        y = y_lab,
        colour = "Site",
        title = paste0(metric, " through sliding windows")
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )

    if (isTRUE(add_points)) {
      p <- p + ggplot2::geom_point(size = 1.4)
    }

    if (isTRUE(add_overall)) {
      overall <- dat %>%
        dplyr::group_by(window_mid) %>%
        dplyr::summarise(value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

      p <- p +
        ggplot2::geom_line(
          data = overall,
          ggplot2::aes(x = window_mid, y = value),
          inherit.aes = FALSE,
          linewidth = 1.2,
          colour = "#C0392B"
        )
    }
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, plot = p, width = 11, height = 7, dpi = 300, bg = "white")
  }

  p
}

#' Plot all core mortAAR trajectories
#'
#' @description
#' Convenience wrapper returning plots for MI, Ratio_F_M, and MMR2.
#'
#' @param metrics_windows Output from `compute_mortaar_site_metrics_windows()`
#'   or `run_mortAAR_bridge()$site_metrics_windows`.
#' @param facet_by_site If TRUE, facet by site.
#' @param site_names Optional character vector of site names to include.
#' @param only_possible If TRUE, keep only rows where the respective metric is available.
#' @param add_points If TRUE, add points to the trajectory lines.
#' @param add_overall If TRUE, add a pooled mean trajectory across included sites.
#' @param save_dir Optional directory for saving all plots.
#'
#' @return Named list of ggplot objects.
plot_mortaar_all_trajectories <- function(metrics_windows,
                                          facet_by_site = TRUE,
                                          site_names = NULL,
                                          only_possible = TRUE,
                                          add_points = TRUE,
                                          add_overall = FALSE,
                                          save_dir = NULL) {
  metrics <- c("MI", "Ratio_F_M", "MMR1", "MMR2")

  plots <- lapply(metrics, function(m) {
    save_path <- NULL
    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      save_path <- file.path(save_dir, paste0("mortaar_trajectory_", m, ".png"))
    }

    plot_mortaar_metric_trajectory(
      metrics_windows = metrics_windows,
      metric = m,
      facet_by_site = facet_by_site,
      site_names = site_names,
      only_possible = only_possible,
      add_points = add_points,
      add_overall = add_overall,
      save_path = save_path
    )
  })

  names(plots) <- metrics
  plots
}

#' Plot one site's mortAAR trajectories through time
#'
#' @description
#' Produces a long-format plot for one site showing MI, Ratio_F_M, and MMR2
#' as separate facets through sliding windows.
#'
#' @param metrics_windows Output from `compute_mortaar_site_metrics_windows()`
#'   or `run_mortAAR_bridge()$site_metrics_windows`.
#' @param site_name_value Site name to plot.
#' @param only_possible If TRUE, keep only rows where the respective metric is available.
#' @param add_points If TRUE, add points.
#' @param save_path Optional file path for saving.
#'
#' @return ggplot object.
plot_mortaar_site_trajectory_panel <- function(metrics_windows,
                                               site_name_value,
                                               only_possible = TRUE,
                                               add_points = TRUE,
                                               save_path = NULL) {
  dat <- tibble::as_tibble(metrics_windows) %>%
    dplyr::filter(site_name == site_name_value)

  if (nrow(dat) == 0) {
    stop("No rows found for the requested site.", call. = FALSE)
  }

  long <- dat %>%
    dplyr::transmute(
      window_mid,
      MI,
      Ratio_F_M,
      MMR1,
      MMR2,
      MI_possible,
      Ratio_F_M_possible,
      MMR1_possible,
      MMR2_possible
    ) %>%
    tidyr::pivot_longer(
      cols = c("MI", "Ratio_F_M", "MMR1", "MMR2"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      possible = dplyr::case_when(
        metric == "MI" ~ MI_possible,
        metric == "Ratio_F_M" ~ Ratio_F_M_possible,
        metric == "MMR1" ~ MMR1_possible,
        metric == "MMR2" ~ MMR2_possible,
        TRUE ~ FALSE
      )
    )

  if (isTRUE(only_possible)) {
    long <- long %>%
      dplyr::filter(possible)
  }

  long <- long %>%
    dplyr::filter(is.finite(value))

  long <- long %>%
    dplyr::filter(!(metric %in% c("MMR1", "MMR2")) | value >= 0)

  if (nrow(long) == 0) {
    stop("No plottable rows remain for the requested site.", call. = FALSE)
  }

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = window_mid, y = value, group = metric)
  ) +
    ggplot2::geom_line(linewidth = 0.8, colour = "#2C3E50") +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = paste("mortAAR trajectories:", site_name_value),
      x = "Window midpoint",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(add_points)) {
    p <- p + ggplot2::geom_point(size = 1.6, colour = "#2C3E50")
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, plot = p, width = 8, height = 9, dpi = 300, bg = "white")
  }

  p
}

# -------------------------------------------------------------------------
# Life-table cluster comparison helpers
# -------------------------------------------------------------------------

#' Parse mortAAR age labels from life-table column `x`
#'
#' @param x Character vector like "0--0", "1--4", ...
#'
#' @return Tibble with age, age_start, age_end, age_mid.
.parse_mortaar_x_age <- function(x) {
  x_chr <- as.character(x)

  age_start <- suppressWarnings(as.numeric(sub("--.*", "", x_chr)))
  age_end   <- suppressWarnings(as.numeric(sub(".*--", "", x_chr)))

  tibble::tibble(
    age = paste0(age_start, "-", age_end),
    age_start = age_start,
    age_end = age_end,
    age_mid = (age_start + age_end) / 2
  )
}

#' Extract one life-table measure across windows
#'
#' @param life_tables_windows Named list of mortAAR life-table objects.
#' @param measure One of "qx", "dx", "lx", "ex", "rel_popx".
#'
#' @return Tibble with window, age, age_start, age_end, age_mid, value.
extract_life_table_measure_windows <- function(life_tables_windows,
                                               measure = c("qx", "dx", "lx", "ex", "rel_popx")) {
  measure <- match.arg(measure)

  if (is.null(life_tables_windows) || length(life_tables_windows) == 0) {
    stop("`life_tables_windows` is empty.", call. = FALSE)
  }

  out <- purrr::map_dfr(names(life_tables_windows), function(nm) {
    lt <- life_tables_windows[[nm]]
    if (is.null(lt)) return(NULL)

    df <- tibble::as_tibble(as.data.frame(lt))

    if (!measure %in% names(df)) return(NULL)
    if (!"x" %in% names(df)) {
      stop("Life table object does not contain column `x` with age labels.", call. = FALSE)
    }

    age_tbl <- .parse_mortaar_x_age(df$x)

    parts <- strsplit(nm, "_")[[1]]
    window_start <- suppressWarnings(as.numeric(parts[1]))
    window_end   <- suppressWarnings(as.numeric(parts[2]))
    window_mid   <- if (is.finite(window_start) && is.finite(window_end)) {
      (window_start + window_end) / 2
    } else {
      NA_real_
    }

    tibble::tibble(
      window = nm,
      window_label = nm,
      window_start = window_start,
      window_end = window_end,
      window_mid = window_mid,
      age = age_tbl$age,
      age_start = age_tbl$age_start,
      age_end = age_tbl$age_end,
      age_mid = age_tbl$age_mid,
      value = suppressWarnings(as.numeric(df[[measure]]))
    )
  })

  out %>%
    dplyr::filter(is.finite(.data$value), is.finite(.data$age_mid))
}

#' Label extracted window table by cluster
#'
#' @param x Extracted life-table measure table.
#' @param cluster_A Character vector of window labels.
#' @param cluster_B Character vector of window labels.
#' @param cluster_names Length-2 character vector naming the clusters.
#'
#' @return Tibble with `cluster` column added.
label_life_table_clusters <- function(x,
                                      cluster_A,
                                      cluster_B,
                                      cluster_names = c("A", "B")) {
  if (length(cluster_names) != 2L) {
    stop("`cluster_names` must have length 2.", call. = FALSE)
  }

  x %>%
    dplyr::mutate(
      cluster = dplyr::case_when(
        .data$window %in% cluster_A ~ cluster_names[1],
        .data$window %in% cluster_B ~ cluster_names[2],
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.data$cluster))
}

#' Compute mean age at death per window from dx
#'
#' @param dx_tbl Output from `label_life_table_clusters()` for measure "dx".
#' @param dx_scale Either "percent" or "proportion".
#'
#' @return Tibble with one row per window.
compute_mean_age_by_window <- function(dx_tbl,
                                       dx_scale = c("percent", "proportion")) {
  dx_scale <- match.arg(dx_scale)

  if (!all(c("window", "cluster", "age_mid", "value") %in% names(dx_tbl))) {
    stop("`dx_tbl` is missing required columns.", call. = FALSE)
  }

  out <- dx_tbl %>%
    dplyr::mutate(
      value_prop = if (dx_scale == "percent") .data$value / 100 else .data$value
    ) %>%
    dplyr::group_by(.data$window, .data$cluster) %>%
    dplyr::summarise(
      sum_dx = sum(.data$value_prop, na.rm = TRUE),
      mean_age = ifelse(sum_dx > 0, sum(.data$age_mid * .data$value_prop, na.rm = TRUE) / sum_dx, NA_real_),
      window_start = dplyr::first(.data$window_start),
      window_end = dplyr::first(.data$window_end),
      window_mid = dplyr::first(.data$window_mid),
      .groups = "drop"
    )

  out
}

#' Run age-specific Wilcoxon tests by cluster
#'
#' @param measure_tbl Output from `label_life_table_clusters()`.
#' @param p_adjust Method for `p.adjust()`.
#'
#' @return Tibble with p-values by age.
compute_age_class_cluster_tests <- function(measure_tbl,
                                            p_adjust = "BH") {
  measure_tbl %>%
    dplyr::group_by(.data$age, .data$age_start, .data$age_end) %>%
    dplyr::summarise(
      p_value = tryCatch(
        stats::wilcox.test(value ~ cluster, data = dplyr::cur_data(), exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(.data$p_value, method = p_adjust)
    ) %>%
    dplyr::arrange(.data$age_start, .data$age_end)
}

#' Run window-level cluster tests on mean age at death
#'
#' @param mean_age_df Output from `compute_mean_age_by_window()`.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#'
#' @return Named list of test outputs.
compute_mean_age_cluster_tests <- function(mean_age_df,
                                           n_perm = 5000,
                                           seed = 123) {
  if (!all(c("mean_age", "cluster") %in% names(mean_age_df))) {
    stop("`mean_age_df` is missing required columns.", call. = FALSE)
  }

  wt <- tryCatch(
    stats::wilcox.test(mean_age ~ cluster, data = mean_age_df, exact = FALSE),
    error = function(e) e
  )

  set.seed(seed)
  cl <- mean_age_df$cluster
  vals <- mean_age_df$mean_age

  obs_diff <- abs(diff(tapply(vals, cl, mean)))
  perm_diffs <- replicate(n_perm, {
    sh <- sample(cl)
    abs(diff(tapply(vals, sh, mean)))
  })
  perm_p <- mean(perm_diffs >= obs_diff)

  list(
    wilcox = wt,
    permutation = list(
      observed_difference = obs_diff,
      permutation_p = perm_p,
      perm_diffs = perm_diffs,
      n_perm = n_perm
    ),
    summary = mean_age_df %>%
      dplyr::group_by(.data$cluster) %>%
      dplyr::summarise(
        n = dplyr::n(),
        unique_vals = dplyr::n_distinct(.data$mean_age),
        min = min(.data$mean_age, na.rm = TRUE),
        max = max(.data$mean_age, na.rm = TRUE),
        mean = mean(.data$mean_age, na.rm = TRUE),
        median = stats::median(.data$mean_age, na.rm = TRUE),
        .groups = "drop"
      )
  )
}

#' Plot difference in one life-table measure between clusters by age
#'
#' @param measure_tbl Output from `label_life_table_clusters()`.
#' @param save_path Optional file path.
#'
#' @return ggplot object.
plot_cluster_measure_difference_by_age <- function(measure_tbl,
                                                   save_path = NULL) {
  cluster_levels <- unique(measure_tbl$cluster)
  if (length(cluster_levels) != 2L) {
    stop("Need exactly two clusters for difference plotting.", call. = FALSE)
  }

  diff_tbl <- measure_tbl %>%
    dplyr::group_by(.data$age, .data$age_start, .data$cluster) %>%
    dplyr::summarise(mean_val = mean(.data$value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = .data$cluster, values_from = .data$mean_val) %>%
    dplyr::mutate(diff = .data[[cluster_levels[2]]] - .data[[cluster_levels[1]]]) %>%
    dplyr::arrange(.data$age_start)

  p <- ggplot2::ggplot(
    diff_tbl,
    ggplot2::aes(x = factor(.data$age, levels = .data$age), y = .data$diff)
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = paste0("Difference in measure between clusters (", cluster_levels[2], " - ", cluster_levels[1], ")"),
      x = "Age range",
      y = "Difference"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  }

  p
}

#' Plot mean age at death by cluster
#'
#' @param mean_age_df Output from `compute_mean_age_by_window()`.
#' @param save_path Optional file path.
#'
#' @return ggplot object.
plot_mean_age_by_cluster <- function(mean_age_df,
                                     save_path = NULL) {
  p <- ggplot2::ggplot(mean_age_df, ggplot2::aes(x = .data$cluster, y = .data$mean_age)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.08, size = 2) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "Mean age at death by cluster",
      x = "Cluster",
      y = "Mean age at death"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 7, height = 5, dpi = 300, bg = "white")
  }

  p
}

#' Plot one life-table measure by cluster and age
#'
#' @param measure_tbl Output from `label_life_table_clusters()`.
#' @param save_path Optional file path.
#'
#' @return ggplot object.
plot_cluster_measure_profiles <- function(measure_tbl,
                                          save_path = NULL) {
  prof_tbl <- measure_tbl %>%
    dplyr::group_by(.data$age, .data$age_start, .data$cluster) %>%
    dplyr::summarise(mean_val = mean(.data$value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(.data$age_start)

  p <- ggplot2::ggplot(
    prof_tbl,
    ggplot2::aes(
      x = factor(.data$age, levels = unique(.data$age)),
      y = .data$mean_val,
      colour = .data$cluster,
      group = .data$cluster
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "Cluster life-table profiles by age",
      x = "Age range",
      y = "Mean value",
      colour = "Cluster"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  }

  p
}

#' Run full cluster comparison workflow for one life-table measure
#'
#' @param life_tables_windows Named list of mortAAR life tables.
#' @param measure One of "qx", "dx", "lx", "ex", "rel_popx".
#' @param cluster_A Character vector of window labels for cluster A.
#' @param cluster_B Character vector of window labels for cluster B.
#' @param cluster_names Length-2 vector naming clusters.
#' @param n_perm Number of permutations for mean-age test.
#' @param seed Random seed.
#' @param export_dir Optional directory to export tables and plots.
#'
#' @return Named list with extracted tables, tests, and plots.
run_life_table_cluster_comparison <- function(life_tables_windows,
                                              measure = c("qx", "dx", "lx", "ex", "rel_popx"),
                                              cluster_A,
                                              cluster_B,
                                              cluster_names = c("A", "B"),
                                              n_perm = 5000,
                                              seed = 123,
                                              export_dir = NULL) {
  measure <- match.arg(measure)

  measure_tbl <- extract_life_table_measure_windows(
    life_tables_windows = life_tables_windows,
    measure = measure
  )

  measure_tbl <- label_life_table_clusters(
    x = measure_tbl,
    cluster_A = cluster_A,
    cluster_B = cluster_B,
    cluster_names = cluster_names
  )

  age_tests <- compute_age_class_cluster_tests(measure_tbl)

  out <- list(
    measure_table = measure_tbl,
    age_tests = age_tests,
    mean_age = NULL,
    mean_age_tests = NULL,
    plots = list()
  )

  out$plots$profiles <- plot_cluster_measure_profiles(measure_tbl)
  out$plots$difference <- plot_cluster_measure_difference_by_age(measure_tbl)

  if (measure == "dx") {
    mean_age_df <- compute_mean_age_by_window(measure_tbl, dx_scale = "percent")
    mean_age_tests <- compute_mean_age_cluster_tests(mean_age_df, n_perm = n_perm, seed = seed)

    out$mean_age <- mean_age_df
    out$mean_age_tests <- mean_age_tests
    out$plots$mean_age <- plot_mean_age_by_cluster(mean_age_df)
  }

  if (!is.null(export_dir)) {
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

    utils::write.csv(
      out$measure_table,
      file.path(export_dir, paste0("cluster_measure_table_", measure, ".csv")),
      row.names = FALSE
    )

    utils::write.csv(
      out$age_tests,
      file.path(export_dir, paste0("cluster_age_tests_", measure, ".csv")),
      row.names = FALSE
    )

    if (!is.null(out$mean_age)) {
      utils::write.csv(
        out$mean_age,
        file.path(export_dir, "cluster_mean_age_dx.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        out$mean_age_tests$summary,
        file.path(export_dir, "cluster_mean_age_dx_summary.csv"),
        row.names = FALSE
      )
    }

    out$plots$profiles <- plot_cluster_measure_profiles(
      measure_tbl,
      save_path = file.path(export_dir, paste0("cluster_profiles_", measure, ".png"))
    )

    out$plots$difference <- plot_cluster_measure_difference_by_age(
      measure_tbl,
      save_path = file.path(export_dir, paste0("cluster_difference_", measure, ".png"))
    )

    if (!is.null(out$mean_age)) {
      out$plots$mean_age <- plot_mean_age_by_cluster(
        out$mean_age,
        save_path = file.path(export_dir, "cluster_mean_age_dx.png")
      )
    }
  }

  out
}

#' Run cluster comparisons for all standard life-table measures
#'
#' @param life_tables_windows Named list of mortAAR life tables.
#' @param cluster_A Character vector of window labels for cluster A.
#' @param cluster_B Character vector of window labels for cluster B.
#' @param cluster_names Length-2 vector naming clusters.
#' @param n_perm Number of permutations for dx mean-age test.
#' @param seed Random seed.
#' @param export_dir Optional base export directory.
#'
#' @return Named list of per-measure results.
run_all_life_table_cluster_comparisons <- function(life_tables_windows,
                                                   cluster_A,
                                                   cluster_B,
                                                   cluster_names = c("A", "B"),
                                                   n_perm = 5000,
                                                   seed = 123,
                                                   export_dir = NULL) {
  measures <- c("qx", "dx", "lx", "ex", "rel_popx")

  out <- lapply(measures, function(m) {
    subdir <- if (is.null(export_dir)) NULL else file.path(export_dir, m)

    run_life_table_cluster_comparison(
      life_tables_windows = life_tables_windows,
      measure = m,
      cluster_A = cluster_A,
      cluster_B = cluster_B,
      cluster_names = cluster_names,
      n_perm = n_perm,
      seed = seed,
      export_dir = subdir
    )
  })

  names(out) <- measures
  out
}


# -------------------------------------------------------------------------
# Sex-specific life-table cluster comparison helpers
# -------------------------------------------------------------------------

#' Extract one life-table measure across windows by sex_gender
#'
#' @param life_tables_windows_by_sex Named list of sex-specific mortAAR life tables.
#' @param measure One of "qx", "dx", "lx", "ex", "rel_popx".
#'
#' @return Tibble with window, sex_gender, age, age_start, age_end, age_mid, value.
extract_life_table_measure_windows_by_sex <- function(life_tables_windows_by_sex,
                                                      measure = c("qx", "dx", "lx", "ex", "rel_popx")) {
  measure <- match.arg(measure)

  if (is.null(life_tables_windows_by_sex) || length(life_tables_windows_by_sex) == 0) {
    stop("`life_tables_windows_by_sex` is empty.", call. = FALSE)
  }

  out <- purrr::map_dfr(names(life_tables_windows_by_sex), function(nm) {
    lt_sex <- life_tables_windows_by_sex[[nm]]
    if (is.null(lt_sex)) return(NULL)

    purrr::map_dfr(names(lt_sex), function(sex_gender) {
      if (sex_gender == "All") return(NULL)

      lt <- lt_sex[[sex_gender]]
      if (is.null(lt)) return(NULL)

      df <- tibble::as_tibble(as.data.frame(lt))

      if (!measure %in% names(df)) return(NULL)
      if (!"x" %in% names(df)) {
        stop("Life table missing `x` column.", call. = FALSE)
      }

      age_tbl <- .parse_mortaar_x_age(df$x)

      parts <- strsplit(nm, "_")[[1]]
      window_start <- suppressWarnings(as.numeric(parts[1]))
      window_end   <- suppressWarnings(as.numeric(parts[2]))
      window_mid   <- if (is.finite(window_start) && is.finite(window_end)) {
        (window_start + window_end) / 2
      } else {
        NA_real_
      }

      tibble::tibble(
        window = nm,
        window_start = window_start,
        window_end = window_end,
        window_mid = window_mid,
        sex_gender = sex_gender,
        age = age_tbl$age,
        age_start = age_tbl$age_start,
        age_end = age_tbl$age_end,
        age_mid = age_tbl$age_mid,
        value = suppressWarnings(as.numeric(df[[measure]]))
      )
    })
  })

  out %>%
    dplyr::filter(is.finite(.data$value), is.finite(.data$age_mid))
}

#' Compute mean age at death per window from dx by sex_gender
#'
#' @param dx_tbl Output from `label_life_table_clusters()` on sex-specific dx table.
#' @param dx_scale Either "percent" or "proportion".
#'
#' @return Tibble with one row per window and sex_gender.
compute_mean_age_by_window_by_sex <- function(dx_tbl,
                                              dx_scale = c("percent", "proportion")) {
  dx_scale <- match.arg(dx_scale)

  if (!all(c("window", "cluster", "sex_gender", "age_mid", "value") %in% names(dx_tbl))) {
    stop("`dx_tbl` is missing required columns.", call. = FALSE)
  }

  dx_tbl %>%
    dplyr::filter(.data$sex_gender %in% c("f", "m", "indet")) %>%
    dplyr::mutate(
      value_prop = if (dx_scale == "percent") .data$value / 100 else .data$value
    ) %>%
    dplyr::group_by(.data$window, .data$cluster, .data$sex_gender) %>%
    dplyr::summarise(
      sum_dx = sum(.data$value_prop, na.rm = TRUE),
      mean_age = ifelse(sum_dx > 0, sum(.data$age_mid * .data$value_prop, na.rm = TRUE) / sum_dx, NA_real_),
      window_start = dplyr::first(.data$window_start),
      window_end = dplyr::first(.data$window_end),
      window_mid = dplyr::first(.data$window_mid),
      .groups = "drop"
    )
}

#' Test interaction: sex × cluster (difference-in-differences)
#'
#' @param measure_tbl Output of extract_life_table_measure_by_sex() + clustering
#' @param n_perm Number of permutations
#' @param seed Random seed
#'
#' @return List with observed effect and permutation p-value
test_sex_cluster_interaction <- function(measure_tbl,
                                         n_perm = 5000,
                                         seed = 123) {

  required <- c("value", "cluster", "sex_gender")
  if (!all(required %in% names(measure_tbl))) {
    stop("Missing required columns.", call. = FALSE)
  }

  # remove indeterminate sex if desired
  df <- measure_tbl %>%
    dplyr::filter(sex_gender %in% c("f", "m"))

  # compute observed effect
  means <- df %>%
    dplyr::group_by(sex_gender, cluster) %>%
    dplyr::summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  delta <- means %>%
    tidyr::pivot_wider(names_from = cluster, values_from = mean_val) %>%
    dplyr::mutate(diff = .data[[2]] - .data[[1]])

  obs_effect <- diff(delta$diff)

  # permutation test
  set.seed(seed)

  perm_effects <- replicate(n_perm, {
    df_perm <- df %>%
      dplyr::mutate(cluster = sample(cluster))

    means_p <- df_perm %>%
      dplyr::group_by(sex_gender, cluster) %>%
      dplyr::summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

    delta_p <- means_p %>%
      tidyr::pivot_wider(names_from = cluster, values_from = mean_val) %>%
      dplyr::mutate(diff = .data[[2]] - .data[[1]])

    diff(delta_p$diff)
  })

  p_val <- mean(abs(perm_effects) >= abs(obs_effect))

  list(
    observed_effect = obs_effect,
    p_value = p_val,
    perm_distribution = perm_effects,
    summary = means
  )
}

#' Age-specific interaction tests by sex_gender and cluster
#'
#' @param measure_tbl Sex-specific extracted and clustered measure table.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#'
#' @return Tibble with age-specific interaction effects and p-values.
test_sex_cluster_interaction_by_age <- function(measure_tbl,
                                                n_perm = 2000,
                                                seed = 123) {
  ages <- unique(measure_tbl$age)

  purrr::map_dfr(ages, function(a) {
    df <- measure_tbl %>%
      dplyr::filter(age == a, sex_gender %in% c("f", "m"))

    if (nrow(df) < 4 || dplyr::n_distinct(df$cluster) < 2 || dplyr::n_distinct(df$sex_gender) < 2) {
      return(tibble::tibble(
        age = a,
        age_start = dplyr::first(df$age_start %||% NA_real_),
        age_end = dplyr::first(df$age_end %||% NA_real_),
        effect = NA_real_,
        p_value = NA_real_
      ))
    }

    res <- tryCatch(
      test_sex_cluster_interaction(df, n_perm = n_perm, seed = seed),
      error = function(e) NULL
    )

    if (is.null(res)) {
      return(tibble::tibble(
        age = a,
        age_start = dplyr::first(df$age_start),
        age_end = dplyr::first(df$age_end),
        effect = NA_real_,
        p_value = NA_real_
      ))
    }

    tibble::tibble(
      age = a,
      age_start = dplyr::first(df$age_start),
      age_end = dplyr::first(df$age_end),
      effect = res$observed_effect,
      p_value = res$p_value
    )
  }) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(.data$p_value, method = "BH")
    ) %>%
    dplyr::arrange(.data$age_start, .data$age_end)
}

#' Plot one life-table measure by sex_gender and cluster
#'
#' @param measure_tbl Sex-specific extracted and clustered measure table.
#' @param save_path Optional file path.
#'
#' @return ggplot object.
plot_measure_by_sex_and_cluster <- function(measure_tbl,
                                            save_path = NULL) {
  prof <- measure_tbl %>%
    dplyr::group_by(.data$age, .data$age_start, .data$cluster, .data$sex_gender) %>%
    dplyr::summarise(mean_val = mean(.data$value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(.data$age_start) %>%
    dplyr::mutate(
      sex_gender = factor(.data$sex_gender, levels = c("f", "m", "indet"))
    )

  age_levels <- prof %>%
    dplyr::distinct(.data$age, .data$age_start) %>%
    dplyr::arrange(.data$age_start) %>%
    dplyr::pull(.data$age)

  p <- ggplot2::ggplot(
    prof,
    ggplot2::aes(
      x = factor(.data$age, levels = age_levels),
      y = .data$mean_val,
      colour = .data$sex_gender,
      group = .data$sex_gender
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~ cluster, ncol = 1, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "Life-table profiles by sex_gender and cluster",
      x = "Age range",
      y = "Mean value",
      colour = "sex_gender"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 10, height = 8, dpi = 300, bg = "white")
  }

  p
}

#' Plot sex_gender × cluster interaction effects by age
#'
#' @param interaction_tbl Output from `test_sex_cluster_interaction_by_age()`.
#' @param save_path Optional file path.
#'
#' @return ggplot object.
plot_sex_cluster_interaction_by_age <- function(interaction_tbl,
                                                save_path = NULL) {
  plot_tbl <- interaction_tbl %>%
    dplyr::mutate(
      sig = .data$p_adj < 0.05
    )

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = factor(.data$age, levels = .data$age), y = .data$effect, fill = .data$sig)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("grey70", "red")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = "sex_gender × cluster interaction by age",
      x = "Age range",
      y = "Interaction effect (Δf - Δm)",
      fill = "Significant"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, p, width = 10, height = 6, dpi = 300, bg = "white")
  }

  p
}

#' Run one sex-specific life-table cluster comparison
#'
#' @param life_tables_windows_by_sex Named list of sex-specific mortAAR life tables.
#' @param measure One of "qx", "dx", "lx", "ex", "rel_popx".
#' @param cluster_A Character vector of window labels for cluster A.
#' @param cluster_B Character vector of window labels for cluster B.
#' @param cluster_names Length-2 vector naming clusters.
#' @param n_perm Number of permutations for interaction tests.
#' @param seed Random seed.
#' @param export_dir Optional directory to export tables and plots.
#'
#' @return Named list with extracted tables, tests, and plots.
run_life_table_cluster_comparison_by_sex <- function(life_tables_windows_by_sex,
                                                     measure = c("qx", "dx", "lx", "ex", "rel_popx"),
                                                     cluster_A,
                                                     cluster_B,
                                                     cluster_names = c("A", "B"),
                                                     n_perm = 2000,
                                                     seed = 123,
                                                     export_dir = NULL) {
  measure <- match.arg(measure)

  measure_tbl <- extract_life_table_measure_windows_by_sex(
    life_tables_windows_by_sex = life_tables_windows_by_sex,
    measure = measure
  )

  measure_tbl <- label_life_table_clusters(
    x = measure_tbl,
    cluster_A = cluster_A,
    cluster_B = cluster_B,
    cluster_names = cluster_names
  )

  by_sex <- split(measure_tbl, measure_tbl$sex_gender) %>%
    purrr::map(function(df_sex) {
      list(
        age_tests = compute_age_class_cluster_tests(df_sex),
        mean_age = if (measure == "dx") compute_mean_age_by_window(df_sex, dx_scale = "percent") else NULL
      )
    })

  sex_difference <- measure_tbl %>%
    dplyr::group_by(.data$age, .data$cluster) %>%
    dplyr::summarise(
      p_value = tryCatch(
        stats::wilcox.test(value ~ sex_gender, data = dplyr::cur_data(), exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(.data$p_value, method = "BH")
    )

  interaction <- test_sex_cluster_interaction(
    measure_tbl,
    n_perm = n_perm,
    seed = seed
  )

  interaction_by_age <- test_sex_cluster_interaction_by_age(
    measure_tbl,
    n_perm = n_perm,
    seed = seed
  )

  plots <- list(
    profiles = plot_measure_by_sex_and_cluster(measure_tbl),
    interaction = plot_sex_cluster_interaction_by_age(interaction_by_age)
  )

  mean_age_by_sex <- NULL
  if (measure == "dx") {
    mean_age_by_sex <- compute_mean_age_by_window_by_sex(measure_tbl, dx_scale = "percent")
    plots$mean_age_by_sex <- ggplot2::ggplot(
      mean_age_by_sex,
      ggplot2::aes(x = .data$cluster, y = .data$mean_age, colour = .data$sex_gender)
    ) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.06), size = 2) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.2) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = "Mean age at death by cluster and sex_gender",
        x = "Cluster",
        y = "Mean age at death",
        colour = "sex_gender"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")
      )
  }

  out <- list(
    measure_table = measure_tbl,
    by_sex = by_sex,
    sex_difference = sex_difference,
    interaction = interaction,
    interaction_by_age = interaction_by_age,
    mean_age_by_sex = mean_age_by_sex,
    plots = plots
  )

  if (!is.null(export_dir)) {
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

    utils::write.csv(
      out$measure_table,
      file.path(export_dir, paste0("cluster_measure_table_by_sex_", measure, ".csv")),
      row.names = FALSE
    )

    utils::write.csv(
      out$sex_difference,
      file.path(export_dir, paste0("cluster_sex_difference_", measure, ".csv")),
      row.names = FALSE
    )

    utils::write.csv(
      out$interaction_by_age,
      file.path(export_dir, paste0("cluster_interaction_by_age_", measure, ".csv")),
      row.names = FALSE
    )

    if (!is.null(out$mean_age_by_sex)) {
      utils::write.csv(
        out$mean_age_by_sex,
        file.path(export_dir, "cluster_mean_age_by_sex_dx.csv"),
        row.names = FALSE
      )
    }

    out$plots$profiles <- plot_measure_by_sex_and_cluster(
      measure_tbl,
      save_path = file.path(export_dir, paste0("cluster_profiles_by_sex_", measure, ".png"))
    )

    out$plots$interaction <- plot_sex_cluster_interaction_by_age(
      out$interaction_by_age,
      save_path = file.path(export_dir, paste0("cluster_interaction_by_age_", measure, ".png"))
    )

    if (!is.null(out$mean_age_by_sex)) {
      ggplot2::ggsave(
        file.path(export_dir, "cluster_mean_age_by_sex_dx.png"),
        out$plots$mean_age_by_sex,
        width = 8,
        height = 5,
        dpi = 300,
        bg = "white"
      )
    }
  }

  out
}

#' Run sex-specific life-table cluster comparisons for all standard measures
#'
#' @param life_tables_windows_by_sex Named list of sex-specific mortAAR life tables.
#' @param cluster_A Character vector of window labels for cluster A.
#' @param cluster_B Character vector of window labels for cluster B.
#' @param cluster_names Length-2 vector naming clusters.
#' @param n_perm Number of permutations for interaction tests.
#' @param seed Random seed.
#' @param export_dir Optional base export directory.
#'
#' @return Named list of per-measure results.
run_all_life_table_cluster_comparisons_by_sex <- function(life_tables_windows_by_sex,
                                                          cluster_A,
                                                          cluster_B,
                                                          cluster_names = c("A", "B"),
                                                          n_perm = 2000,
                                                          seed = 123,
                                                          export_dir = NULL) {
  measures <- c("qx", "dx", "lx", "ex", "rel_popx")

  out <- lapply(measures, function(m) {
    subdir <- if (is.null(export_dir)) NULL else file.path(export_dir, m)

    run_life_table_cluster_comparison_by_sex(
      life_tables_windows_by_sex = life_tables_windows_by_sex,
      measure = m,
      cluster_A = cluster_A,
      cluster_B = cluster_B,
      cluster_names = cluster_names,
      n_perm = n_perm,
      seed = seed,
      export_dir = subdir
    )
  })

  names(out) <- measures
  out
}
