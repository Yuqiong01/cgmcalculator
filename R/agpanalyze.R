# ================================ AGP analyze =================================
#' Build AGP-Style 24-Hour CGM Profiles and Run Optional Functional Regression
#'
#' @description
#' Read per-subject CGM CSV files and construct a complete 24-hour
#' time-of-day glucose profile on a fixed grid (e.g., every 5 minutes,
#' 00:00..23:55). Timestamps are converted to minutes-of-day and values
#' are aggregated across days for each subject (AGP-style profile),
#' using one of: (i) mean within floor-time bins, (ii) nearest observation
#' to each grid center, or (iii) linear interpolation evaluated at grid centers.
#'
#' Optionally merge subject-level covariates by id, create an hour-scale
#' matrix sampled at 00:30..23:30 (after optional rolling-mean smoothing
#' on the grid), and optionally fit a functional regression model via
#' refund::pffr on the hour-scale matrix.
#'
#' When a pffr model is fitted, smooth-term coefficient curves can be
#' extracted and plotted on the original y-index grid and, optionally,
#' linearly interpolated to a finer time grid for plotting or saving.
#'
#' @details
#' ## Input format and behavior
#' Each CGM CSV file in \code{inputdir} must contain columns
#' \code{timestamp} and \code{sensorglucose}, plus a subject id column.
#' If \code{cgm_id_col} is not present, the first column in the file is
#' assumed to be the id and will be renamed to \code{cgm_id_col}.
#'
#' Timestamps are parsed with \code{lubridate::parse_date_time()} using a
#' flexible set of common date-time orders. If timestamps cannot be parsed
#' in a given file, the function stops.
#'
#' A complete 24-hour grid is built as bin centers from \code{00:00} to
#' \code{(1440 - grid_mins)} minutes, i.e., labels corresponding to
#' \code{seq(0, 1440 - grid_mins, by = grid_mins)}.
#'
#' If the inferred raw sampling interval is coarser than \code{grid_mins},
#' the function always interpolates to grid centers using
#' \code{stats::approx(rule = upsample_rule)}, regardless of
#' \code{raw_to_grid_strategy}.
#'
#' ## Computational notes
#' Fitting \code{refund::pffr()} can be computationally expensive. Runtime and
#' memory usage generally increase with the number of subjects, the response-grid
#' length (number of y-index points), the number of model terms, and the basis
#' dimensions used for smoothing (as specified in \code{pffr_bs_yindex}).
#' Increasing the number of basis functions can therefore lead to substantially
#' longer fitting times.
#'
#' For faster prototyping, consider starting with fewer terms and/or a coarser
#' response grid (e.g., the default 24-point hour-scale), then refining settings
#' for the final model. Using \code{pffr_algorithm = "bam"} and
#' \code{pffr_discrete = TRUE} may also improve performance for larger datasets.
#'
#' @param inputdir Directory containing per-subject CGM CSV files.
#'   Each file must contain columns: \code{timestamp} and \code{sensorglucose},
#'   plus an id column.
#' @param outputdir Output directory for any files written by this function.
#'   Defaults to \code{tempdir()}. When \code{make_plot=TRUE} and
#'   \code{plot_save=TRUE} and/or \code{coef_save=TRUE}, any *bare file names*
#'   (e.g., \code{"pffr_coef_plot.pdf"}) will be saved under \code{outputdir}.
#'   If a full/relative path is provided (contains \code{/} or \code{\\}),
#'   the path is respected and \code{outputdir} is not prepended.
#' @param covariate_file Optional CSV path for covariates to merge by id.
#'
#' @param tz Timezone used to parse and format timestamps (e.g., "UTC").
#'   Default is "UTC".
#' @param grid_mins Target time-of-day grid resolution in minutes (e.g., 5).
#'   Must divide 1440 exactly. Default is 5.
#' @param cgm_id_col Column name for subject id in CGM files. If missing,
#'   the first column is assumed to be the id. Default is "subjectid".
#' @param cov_id_col Column name for subject id in the covariate file.
#'   Default is "subjectid".
#' @param keep_cov_order If TRUE, keep covariate row order and align CGM by id;
#'   otherwise keep only common ids and align. Default is TRUE.
#' @param raw_to_grid_strategy Strategy when raw interval is finer than grid_mins:
#'   \code{"mean"} (average within each floor-time bin),
#'   \code{"nearest"} (closest observation to each grid center),
#'   or \code{"linear"} (linearly interpolate and evaluate at grid centers).
#'   Default is \code{"mean"}.
#' @param upsample_rule Rule for interpolation when raw interval is coarser
#'   than \code{grid_mins} (or when \code{raw_to_grid_strategy = "linear"}),
#'   passed to \code{stats::approx(rule = ...)}. Use 1 to avoid extrapolation
#'   (outside range -> NA), 2 to extrapolate. Default is 1.
#'
#' @param make_famm If TRUE, create an hour-scale matrix sampled at centers
#'   00:30..23:30. Default is TRUE.
#' @param smooth_k Integer. Window length (in number of grid points) for
#'   rolling-mean smoothing on the time-of-day grid before sampling the
#'   hour-scale matrix. For example, when \code{grid_mins = 5}, \code{smooth_k = 12}
#'   corresponds to a 60-minute window. Default is 12.
#'
#' @param pffr_filter Optional. Subset rule applied only when \code{run_pffr=TRUE}.
#'   Provide a list in one of two modes:
#'   \itemize{
#'     \item \strong{Expression mode}: \code{list(expr = <condition>, drop_na = TRUE)}.
#'       The condition is evaluated within the pffr data (e.g.,
#'       \code{list(expr = rlang::expr(a > 5))} or
#'       \code{list(expr = quote(a > 5))}).
#'     \item \strong{Membership mode}: \code{list(var = "col", keep = <values>, drop_na = TRUE)},
#'       which keeps rows where \code{data[[var]] \%in\% keep}.
#'   }
#' @param run_pffr If TRUE, fit a functional regression model via \code{refund::pffr}
#'   on the hour-scale matrix. Requires \code{make_famm = TRUE}. Default is FALSE.
#' @param pffr_group Name of the exposure/group variable in covariates/data.
#'   Required when \code{run_pffr = TRUE}.
#' @param pffr_covars Optional character vector of additional covariate names
#'   to adjust for. If NULL, no additional covariates are included.
#'   Only used when \code{run_pffr = TRUE}. Default is NULL.
#' @param pffr_id_re Optional subject id variable name for adding a random-effect
#'   smooth \code{s(id, bs = "re")}. If NULL, no random-effect term is added.
#'   Only used when \code{run_pffr = TRUE}. Default is NULL.
#' @param pffr_yind Optional numeric vector passed to \code{refund::pffr}
#'   (length = ncol(hour-scale matrix)). If NULL, uses hour centers
#'   \code{seq(0.5, 23.5, length.out = ncol(MIMS_hour_mat))}.
#' @param pffr_bs_yindex Passed to \code{refund::pffr} as \code{bs.yindex}.
#'   Only used when \code{run_pffr = TRUE}. Default uses a cyclic P-spline over 24 points.
#' @param pffr_algorithm Passed to \code{refund::pffr} (e.g., "bam").
#'   Only used when \code{run_pffr = TRUE}. Default is "bam".
#' @param pffr_discrete Passed to \code{refund::pffr}.
#'   Only used when \code{run_pffr = TRUE}. Default is TRUE.
#'
#' @param make_plot If TRUE, extract smooth-term coefficient curves and
#'   build ggplot objects. Only applicable when \code{run_pffr = TRUE}.
#'   Default behavior: if user does not specify, then it follows \code{run_pffr}
#'   (\code{run_pffr=TRUE -> make_plot=TRUE; run_pffr=FALSE -> make_plot=FALSE}).
#' @param plot_terms Character vector of terms to plot (without "(yind)" suffix).
#'   If NULL, plot all eligible smooth terms. Default is NULL.
#' @param plot_drop_intercept Drop the \code{Intercept(yind)} term from plots.
#'   Default is FALSE.
#' @param plot_drop_re Drop random-effect smooth terms (e.g., \code{s(id)}).
#'   Default is TRUE.
#' @param plot_drop_yind_smooth Drop baseline y-index smooth terms such as \code{s(yind)}.
#'   Default is TRUE.
#' @param coef_grid_mins Optional; if not NULL, interpolate coefficient curves
#'   to this grid in minutes (e.g., 5). Must divide 1440 exactly.
#'   Default behavior: if user does not specify and \code{run_pffr=TRUE} & \code{make_plot=TRUE},
#'   then \code{coef_grid_mins} defaults to 5 to produce the interpolated page/table.
#' @param plot_save If TRUE, save plots to files. Only applicable when
#'   \code{make_plot = TRUE}. Default behavior: follows \code{make_plot}.
#' @param plot_file File path for saving the coefficient plot (PDF). When a bare
#'   file name is provided (no \code{/} or \code{\\}), it is saved under
#'   \code{outputdir}. Default is \code{"pffr_coef_plot.pdf"}.
#'
#' @param coef_save If TRUE, save coefficient tables as CSV (only when
#'   \code{make_plot = TRUE}). File names are generated automatically as
#'   \code{pffr_coef_df_<group_var>.csv}, where \code{<group_var>} is the name of
#'   \code{pffr_group}; if \code{coef_grid_mins} is set and interpolation succeeds,
#'   an additional file \code{pffr_coef_df_grid<mins>min_<group_var>.csv} is also created.
#'   Default behavior: follows \code{make_plot}.
#'
#' @param pffr_summary_save If TRUE, save \code{summary(res$pffr_fit)} to a text file.
#'   Only applicable when \code{run_pffr = TRUE}. Default behavior: follows \code{run_pffr}.
#' @param pffr_summary_file File path for saving the pffr summary. When a bare
#'   file name is provided, the group variable name (\code{pffr_group}) is appended
#'   to the base name (before the extension) and the file is saved under
#'   \code{outputdir}. Default is \code{"pffr_fit_summary.txt"}.
#'
#' @param verbose Logical; if \code{TRUE}, emit informative messages (via
#'   \code{message()}) about output paths and files written. Default is
#'   \code{FALSE} so the function is quiet by default and does not write to
#'   the console unless explicitly requested. This argument does not affect
#'   the returned results.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{cgm_df}: Subject-level tibble containing \code{MIMS_mat} and \code{MIMS_tf}
#'     (and covariates if provided).
#'   \item \code{cgm_famm_df}: Hour-scale tibble containing \code{MIMS_hour_mat} and hour-wide columns
#'     when \code{make_famm = TRUE}; otherwise \code{NULL}.
#'   \item \code{pffr_fit}: Fitted \code{refund::pffr} model when \code{run_pffr = TRUE}; otherwise \code{NULL}.
#'   \item \code{pffr_call}: Bookkeeping list describing the pffr call and output paths when applicable.
#'   \item \code{pffr_coef_df}: Coefficient table on the original y-index grid when \code{make_plot = TRUE};
#'     otherwise \code{NULL}.
#'   \item \code{pffr_coef_df_grid}: Interpolated coefficient table when \code{coef_grid_mins} is not \code{NULL}
#'     and interpolation succeeds; otherwise \code{NULL}.
#'   \item \code{pffr_plot}: ggplot object for coefficient curves on the original grid when available; otherwise \code{NULL}.
#'   \item \code{pffr_plot_grid}: ggplot object for coefficient curves on the interpolated grid when available; otherwise \code{NULL}.
#'   \item Additional fields for file/grid bookkeeping, including \code{outputdir}, \code{plot_file},
#'     \code{coef_file}, \code{coef_file_grid}, \code{pffr_summary_file}, \code{grid_ok}, and \code{grid_reason}.
#' }
#'
#' @section Output files:
#' When \code{run_pffr=TRUE} and the user provides \code{pffr_group}, the default behavior
#' is to produce and save:
#' \itemize{
#'   \item Two PDFs: one for all eligible terms and one for the group term only; each PDF contains
#'     the original coefficient curve page and (if interpolation succeeds) the interpolated page.
#'   \item Coefficient CSV tables on original and interpolated grids (if enabled / succeeds).
#'   \item A pffr summary text file.
#' }
#' When \code{plot_terms} is provided, the function saves plots for the selected terms only.
#' When \code{coef_grid_mins = NULL}, the interpolated page/table are disabled.
#'
#' @examples
#' \dontrun{
#' inputdir  <- system.file("exdata", "tidydata", package = "cgmcalculator")
#' outputdir <- system.file("exdata", package = "cgmcalculator")
#' covariate_file <- system.file("exdata", "covariates.csv", package = "cgmcalculator")
#'
#' # ---- run with minimal args  ----
#' res <- agpanalyze(
#'   inputdir       = inputdir,
#'   outputdir      = outputdir,
#'   covariate_file = covariate_file,
#'   run_pffr       = TRUE,
#'   pffr_bs_yindex = list(bs = "cp", k = 12, m = c(2, 1)),
#'   pffr_group     = "micro",
#'   pffr_covars    = c("age","male","education","dm_duration","egfr","tc","tg"),
#'   pffr_id_re     = "subjectid",
#'   verbose        = TRUE
#' )
#' }
#'
#' @seealso \link[refund:pffr]{refund::pffr}
#' @seealso \code{tf::tfd}
#' @export
agpanalyze <- function(
                        # 1) I/O
                        inputdir,
                        outputdir = tempdir(),
                        covariate_file = NULL,

                        # 2) parsing / grid
                        tz = "UTC",
                        grid_mins = 5,
                        cgm_id_col = "subjectid",
                        cov_id_col = "subjectid",
                        keep_cov_order = TRUE,
                        raw_to_grid_strategy = c("mean", "nearest", "linear"),
                        upsample_rule = 1,

                        # 3) hour-matrix (FAMM-style)
                        make_famm = TRUE,
                        smooth_k = 12,

                        # 4) pffr modeling
                        pffr_filter = NULL,
                        run_pffr = FALSE,
                        pffr_group = NULL,
                        pffr_covars = NULL,
                        pffr_id_re = NULL,
                        pffr_yind = NULL,
                        pffr_bs_yindex = list(bs = "cp", k = 24, m = c(2, 1)),
                        pffr_algorithm = "bam",
                        pffr_discrete = TRUE,

                        # 5) plotting / saving (AUTO defaults: follow run_pffr unless user overrides)
                        make_plot = NA,
                        plot_terms = NULL,
                        plot_drop_intercept = FALSE,
                        plot_drop_re = TRUE,
                        plot_drop_yind_smooth = TRUE,
                        coef_grid_mins = NA,
                        plot_save = NA,
                        plot_file = "pffr_coef_plot.pdf",
                        coef_save = NA,

                        # 6) summary output (AUTO defaults: follow run_pffr unless user overrides)
                        pffr_summary_save = NA,
                        pffr_summary_file = "pffr_fit_summary.txt",

                        # 7) misc
                        verbose = FALSE
                      ) {

  # -------------------- auto-default normalization --------------------

  # make_plot: NA means "auto"
  if (length(make_plot) == 1L && is.na(make_plot)) {
    make_plot <- isTRUE(run_pffr)
  }
  if (!is.logical(make_plot) || length(make_plot) != 1L) {
    stop("`make_plot` must be TRUE/FALSE (or NA for auto).")
  }

  # plot_save: NA means "auto follow make_plot"
  if (length(plot_save) == 1L && is.na(plot_save)) {
    plot_save <- isTRUE(make_plot)
  }
  if (!is.logical(plot_save) || length(plot_save) != 1L) {
    stop("`plot_save` must be TRUE/FALSE (or NA for auto).")
  }

  # coef_save: NA means "auto follow make_plot"
  if (length(coef_save) == 1L && is.na(coef_save)) {
    coef_save <- isTRUE(make_plot)
  }
  if (!is.logical(coef_save) || length(coef_save) != 1L) {
    stop("`coef_save` must be TRUE/FALSE (or NA for auto).")
  }

  # pffr_summary_save: NA means "auto follow run_pffr"
  if (length(pffr_summary_save) == 1L && is.na(pffr_summary_save)) {
    pffr_summary_save <- isTRUE(run_pffr)
  }
  if (!is.logical(pffr_summary_save) || length(pffr_summary_save) != 1L) {
    stop("`pffr_summary_save` must be TRUE/FALSE (or NA for auto).")
  }

  # coef_grid_mins: NA means "auto default to 5 when run_pffr & make_plot"
  # NULL means "explicitly disable interpolation page/table"
  if (length(coef_grid_mins) == 1L && is.na(coef_grid_mins)) {
    if (isTRUE(run_pffr) && isTRUE(make_plot))
      coef_grid_mins <- 5
    else
      coef_grid_mins <- NULL
  }

  # -------------------- state holders --------------------
  pffr_fit <- NULL
  pffr_call <- NULL
  pffr_coef_df <- NULL
  pffr_coef_df_grid <- NULL
  pffr_plot <- NULL
  pffr_plot_grid <- NULL

  grid_ok <- NA
  grid_reason <- NA_character_
  plot_file_grid_written <- FALSE
  coef_file_grid_written <- FALSE
  coef_file <- NULL
  coef_file_grid <- NULL
  pffr_summary_file_written <- NULL

  # -------------------- helpers --------------------
  .is_bare_filename <- function(x) {
    is.character(x) && length(x) == 1L && nzchar(x) && !grepl("(/|\\\\)", x)
  }

  .rollmean_extend <- function(x, k) {
    if (length(x) < k) return(rep(NA_real_, length(x)))
    w <- rep(1 / k, k)
    y <- as.numeric(stats::filter(x, filter = w, sides = 2))
    if (all(is.na(y))) return(y)
    first_ok <- which(!is.na(y))[1L]
    last_ok  <- utils::tail(which(!is.na(y)), 1L)
    if (first_ok > 1L) y[1L:(first_ok - 1L)] <- y[first_ok]
    if (last_ok  < length(y)) y[(last_ok + 1L):length(y)] <- y[last_ok]
    y
  }

  .hour_to_hhmm <- function(hour) {
    h <- as.numeric(hour)
    mins <- as.integer(round(h * 60))
    mins[mins < 0] <- 0L
    mins[mins > 1440] <- 1440L
    hh <- mins %/% 60L
    mm <- mins %% 60L
    sprintf("%02d:%02d", hh, mm)  # 1440 -> "24:00"
  }

  .default_plot_file_all <- function(group) {
    paste0("pffr_coef_plot_all_", group, ".pdf")
  }

  .default_plot_file_group <- function(group) {
    paste0("pffr_coef_plot_", group, ".pdf")
  }

  .default_coef_file <- function(group) {
    paste0("pffr_coef_df_", group, ".csv")
  }

  .default_coef_file_grid <- function(group, coef_grid_mins) {
    paste0("pffr_coef_df_grid", coef_grid_mins, "min_", group, ".csv")
  }

  .derive_plot_file_grid <- function(plot_file, coef_grid_mins) {
    if (is.null(coef_grid_mins)) return(NULL)
    ext <- tools::file_ext(plot_file)
    base <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", plot_file) else plot_file
    suf <- paste0("_grid", coef_grid_mins, "min")
    if (nzchar(ext)) paste0(base, suf, ".", ext) else paste0(base, suf)
  }

  .draw_ggplot_to_device <- function(p) {
    if (is.null(p)) return(invisible(FALSE))
    if (!inherits(p, "ggplot")) stop("Expected a ggplot object.")
    grid::grid.newpage()
    grid::grid.draw(ggplot2::ggplotGrob(p))
    invisible(TRUE)
  }

  .save_two_pages_pdf <- function(file, plot1, plot2 = NULL,
                                  width = 10, height = 5, useDingbats = FALSE) {
    if (is.null(plot1) && is.null(plot2)) return(invisible(FALSE))

    grDevices::pdf(file = file, onefile = TRUE,
                   width = width, height = height, useDingbats = useDingbats)
    on.exit(grDevices::dev.off(), add = TRUE)

    if (!is.null(plot1)) .draw_ggplot_to_device(plot1)
    if (!is.null(plot2)) .draw_ggplot_to_device(plot2)

    invisible(TRUE)
  }

  # For interpolation stability only: close each term at 0/24
  .close_ends <- function(df, x = "hour") {
    x <- rlang::as_name(rlang::ensym(x))

    out <- df |>
      dplyr::group_by(.data$term) |>
      dplyr::arrange(.data[[x]]) |>
      dplyr::summarise(
        x_tmp = c(0, .data[[x]], 24),
        coef  = c(dplyr::first(.data$coef), .data$coef, dplyr::last(.data$coef)),
        lb    = c(dplyr::first(.data$lb),   .data$lb,   dplyr::last(.data$lb)),
        ub    = c(dplyr::first(.data$ub),   .data$ub,   dplyr::last(.data$ub)),
        .groups = "drop"
      )

    names(out)[names(out) == "x_tmp"] <- x
    out
  }

  .plot_termcurve_df <- function(df, title = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for plotting.")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required for plotting.")
    }
    if (!requireNamespace("rlang", quietly = TRUE)) {
      stop("Package 'rlang' is required for plotting.")
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
      stop("Package 'grid' is required for plotting.")
    }

    # ---- coerce types ----
    if ("term" %in% names(df)) df$term <- as.character(df$term)
    for (nm in intersect(c("yind", "hour", "coef", "lb", "ub"), names(df))) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }

    # ---- ensure hour exists ----
    if (!("hour" %in% names(df))) {
      if (!("yind" %in% names(df))) stop("plot df must have `hour` or `yind`.")
      # full-chain convention: yind is hour centers (0.5..23.5)
      df$hour <- as.numeric(df$yind)
    }

    # ---- ensure time label exists ----
    if (!("time" %in% names(df))) {
      .hour_to_hhmm_local <- function(hour) {
        h <- as.numeric(hour)
        mins <- as.integer(round(h * 60))
        mins[mins < 0] <- 0L
        mins[mins > 1440] <- 1440L
        hh <- mins %/% 60L
        mm <- mins %% 60L
        sprintf("%02d:%02d", hh, mm)  # 1440 -> "24:00"
      }
      df$time <- .hour_to_hhmm_local(df$hour)
    }

    # ---- stable facet/order ----
    if (!is.factor(df$term)) {
      df$term <- factor(as.character(df$term), levels = unique(as.character(df$term)))
    }

    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$hour,
        y = .data$coef,
        color = .data$term,
        fill = .data$term
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.35,
        colour = "grey35"
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lb, ymax = .data$ub),
        alpha = 0.22,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.65, show.legend = TRUE) +
      ggplot2::facet_wrap(~ .data$term, scales = "free_y") +
      ggplot2::scale_x_continuous(
        breaks = seq(0, 24, by = 6),
        minor_breaks = seq(0, 24, by = 1),
        labels = function(x) sprintf("%02d:00", as.integer(x)),
        limits = c(0, 24),
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        title = title,
        x = "Time of day (24 hour clock)",
        y = "Difference in glucose (mg/dL)"
      ) +
      # NOTE: these are ggplot2's viridis scales (no viridis::scale_*_viridis_d needed)
      ggplot2::scale_colour_viridis_d(option = "D", end = 1) +
      ggplot2::scale_fill_viridis_d(option = "D", end = 1) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 14,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 8)
        ),
        strip.text = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        legend.position = "right",
        legend.title = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "grey35", fill = NA, linewidth = 0.4),
        axis.line = ggplot2::element_line(linewidth = 0.35, colour = "grey35"),
        axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "grey35"),
        axis.ticks.length = grid::unit(0.1, "cm"),
        axis.text = ggplot2::element_text(colour = "grey20"),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
        panel.grid.minor = ggplot2::element_line(linewidth = 0.18, colour = "grey95")
      )
  }

  .upsample_termcurve_grid <- function(df, grid_mins) {
    if (is.null(grid_mins)) return(NULL)
    grid_h <- seq(0, 24, by = grid_mins / 60)

    # ---- preserve term order from incoming df ----
    term_levels <- if (is.factor(df$term)) levels(df$term) else unique(as.character(df$term))

    df$term <- if (is.factor(df$term)) df$term else factor(as.character(df$term), levels = term_levels)
    df$hour <- as.numeric(df$hour)
    df$coef <- as.numeric(df$coef)
    df$lb   <- as.numeric(df$lb)
    df$ub   <- as.numeric(df$ub)

    df2 <- .close_ends(df, x = "hour")
    # keep the same factor levels after close_ends()
    df2$term <- factor(as.character(df2$term), levels = term_levels)

    # ---- iterate in the exact same order as term_levels ----
    out_list <- lapply(term_levels, function(tt) {
      d <- df2[df2$term == tt, , drop = FALSE]
      d <- d[order(d$hour), , drop = FALSE]

      # ensure x is strictly increasing for approx()
      if (any(duplicated(d$hour))) {
        d <- dplyr::as_tibble(d) |>
          dplyr::group_by(.data$hour) |>
          dplyr::summarise(
            term = dplyr::first(.data$term),
            coef = mean(.data$coef, na.rm = TRUE),
            lb   = mean(.data$lb,   na.rm = TRUE),
            ub   = mean(.data$ub,   na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(.data$hour) |>
          as.data.frame()
      }

      data.frame(
        term = tt,
        hour = grid_h,
        time = .hour_to_hhmm(grid_h),
        coef = stats::approx(d$hour, d$coef, xout = grid_h, method = "linear", rule = 2)$y,
        lb   = stats::approx(d$hour, d$lb,   xout = grid_h, method = "linear", rule = 2)$y,
        ub   = stats::approx(d$hour, d$ub,   xout = grid_h, method = "linear", rule = 2)$y,
        stringsAsFactors = FALSE
      )
    })

    df_grid <- dplyr::bind_rows(out_list)
    df_grid$term <- factor(df_grid$term, levels = term_levels)
    df_grid
  }

  .reorder_coef_df_raw <- function(df) {
    if (is.null(df)) return(NULL)
    need <- c("term","yind","hour","time","coef","se","lb","ub")
    miss <- setdiff(need, names(df))
    if (length(miss)) stop("pffr_coef_df missing columns: ", paste(miss, collapse = ", "))
    df[, need, drop = FALSE]
  }

  .reorder_coef_df_grid <- function(df) {
    if (is.null(df)) return(NULL)
    need <- c("term","hour","time","coef","lb","ub")
    miss <- setdiff(need, names(df))
    if (length(miss)) stop("pffr_coef_df_grid missing columns: ", paste(miss, collapse = ", "))
    df[, need, drop = FALSE]
  }

  .termcurve_from_pffr_smterms <- function(pffr_fit,
                                           terms = NULL,
                                           drop_intercept = TRUE,
                                           drop_re = TRUE,
                                           drop_yind_smooth = TRUE,
                                           clean_term_label = TRUE) {
    sm <- stats::coef(pffr_fit)$smterms
    if (is.null(sm) || !length(sm)) stop("pffr_fit has no smterms in coef().")

    term_all <- names(sm)
    suffix <- if (any(grepl("\\(yindex\\)$", term_all))) "yindex" else "yind"

    if (isTRUE(drop_re)) {
      term_all <- term_all[!grepl("^s\\(.+\\)$", term_all)]
    }
    if (isTRUE(drop_intercept)) {
      term_all <- term_all[!grepl(paste0("^Intercept\\(", suffix, "\\)$"), term_all)]
    }
    if (isTRUE(drop_yind_smooth)) {
      term_all <- term_all[!grepl(paste0("^s\\(", suffix, "\\)$"), term_all)]
      term_all <- term_all[!grepl(paste0("^s\\(", suffix, ",.*\\)$"), term_all)]
    }

    normalize_terms <- function(t0) {
      t0 <- as.character(t0)
      t0[t0 %in% c("(Intercept)", "Intercept")] <- "Intercept"
      t0[t0 %in% c("exposure", "group")] <- "pffr_group_internal"
      if (!is.null(pffr_group) && is.character(pffr_group) && length(pffr_group) == 1L) {
        t0[t0 == pffr_group] <- "pffr_group_internal"
      }
      # remove (yind)/(yindex) if user already included
      t0 <- sub("\\((yind|yindex)\\)$", "", t0)
      paste0(t0, "(", suffix, ")")
    }

    if (!is.null(terms) && length(terms) >= 1L) {
      t_full <- normalize_terms(terms)
      term_use <- intersect(term_all, t_full)

      if (!length(term_use)) {
        stop(
          "None of requested terms found in smterms.\n",
          "Requested (after normalize): ", paste(t_full, collapse = ", "), "\n",
          "Available: ", paste(term_all, collapse = ", ")
        )
      }
    } else {
      term_use <- term_all
    }

    out <- dplyr::bind_rows(lapply(term_use, function(tt) {
      tbl <- sm[[tt]]$coef
      df <- tibble::as_tibble(tbl)

      x_candidates <- intersect(names(df), c("yind.vec", "yindex.vec", "yind", "yindex", "epoch", "arg"))
      if (!length(x_candidates)) {
        stop("Cannot find y-index column for term: ", tt,
             "\nAvailable columns: ", paste(names(df), collapse = ", "))
      }
      xcol <- x_candidates[1L]

      if (!("value" %in% names(df))) stop("No 'value' column for term: ", tt)
      if (!("se" %in% names(df))) stop("No 'se' column for term: ", tt)

      dplyr::transmute(
        df,
        term = tt,
        yind = as.numeric(.data[[xcol]]),
        coef = as.numeric(.data$value),
        se   = as.numeric(.data$se),
        lb   = as.numeric(.data$value) - 1.96 * as.numeric(.data$se),
        ub   = as.numeric(.data$value) + 1.96 * as.numeric(.data$se)
      )
    }))

    out <- out |>
      dplyr::mutate(
        hour = as.numeric(.data$yind),
        time = .hour_to_hhmm(.data$hour)
      )

    if (isTRUE(clean_term_label)) {
      out <- out |>
        dplyr::mutate(
          term = sub("\\((yind|yindex)\\)$", "", .data$term),
          term = dplyr::if_else(.data$term == "pffr_group_internal", "group", .data$term)
        )
    }

    out
  }

  # --- build a complete time-of-day grid (bin centers) ---
  .grid_centers_minutes <- function(grid_mins) {
    seq(0, 1440 - grid_mins, by = grid_mins)
  }

  .minutes_to_hhmm <- function(mins) {
    mins <- as.integer(round(mins))
    mins[mins < 0] <- 0L
    mins[mins > 1440] <- 1440L
    hh <- mins %/% 60L
    mm <- mins %% 60L
    sprintf("%02d:%02d", hh, mm)
  }

  # --- robust per-subject aggregation onto full grid ---
  .aggregate_subject_to_grid <- function(df,
                                         tz,
                                         grid_mins,
                                         raw_to_grid_strategy = c("mean","nearest","linear"),
                                         upsample_rule = 1) {
    raw_to_grid_strategy <- match.arg(raw_to_grid_strategy)

    if (!("timestamp" %in% names(df))) stop("Input CGM data must contain 'timestamp'.")
    if (!("sensorglucose" %in% names(df))) stop("Input CGM data must contain 'sensorglucose'.")

    df <- df[!is.na(df$timestamp) & !is.na(df$sensorglucose), , drop = FALSE]
    if (!nrow(df)) return(NULL)

    # compute minutes-of-day
    lt <- as.POSIXlt(df$timestamp, tz = tz)
    mod_min <- lt$hour * 60 + lt$min + lt$sec / 60
    mod_min <- mod_min %% 1440

    o <- order(mod_min)
    mod_min_o <- mod_min[o]
    y_o <- as.numeric(df$sensorglucose[o])

    # Pre-aggregate duplicates for interpolation modes:
    if (raw_to_grid_strategy %in% c("linear","nearest")) {
      tmp <- data.frame(mod_min = mod_min_o, y = y_o)
      tmp <- tmp |>
        dplyr::group_by(.data$mod_min) |>
        dplyr::summarise(y = mean(.data$y, na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(.data$mod_min)
      mod_min_o <- tmp$mod_min
      y_o <- tmp$y
    }

    dmins <- diff(mod_min_o)
    dmins <- dmins[is.finite(dmins) & dmins > 0]
    step_min <- if (length(dmins)) stats::median(dmins) else NA_real_

    grid_centers <- .grid_centers_minutes(grid_mins)
    grid_labels <- .minutes_to_hhmm(grid_centers)

    # CASE 1: raw coarser than grid => interpolate to grid centers
    if (!is.na(step_min) && step_min > grid_mins) {
      y_at_centers <- stats::approx(
        x = mod_min_o,
        y = y_o,
        xout = grid_centers,
        method = "linear",
        rule = upsample_rule
      )$y
      return(data.frame(time = grid_labels, value = y_at_centers, stringsAsFactors = FALSE))
    }

    # CASE 2: raw finer/equal to grid => explicit strategy
    if (raw_to_grid_strategy == "linear") {
      y_at_centers <- stats::approx(
        x = mod_min_o,
        y = y_o,
        xout = grid_centers,
        method = "linear",
        rule = upsample_rule
      )$y
      return(data.frame(time = grid_labels, value = y_at_centers, stringsAsFactors = FALSE))
    }

    if (raw_to_grid_strategy == "nearest") {
      y_near <- vapply(grid_centers, function(c0) {
        j <- which.min(abs(mod_min - c0))
        as.numeric(df$sensorglucose[j])
      }, numeric(1L))
      return(data.frame(time = grid_labels, value = y_near, stringsAsFactors = FALSE))
    }

    # raw_to_grid_strategy == "mean":
    bin_left <- floor(mod_min / grid_mins) * grid_mins
    bin_left[bin_left >= 1440] <- 0

    agg <- data.frame(bin = bin_left, y = as.numeric(df$sensorglucose), stringsAsFactors = FALSE) |>
      dplyr::group_by(.data$bin) |>
      dplyr::summarise(value = mean(.data$y, na.rm = TRUE), .groups = "drop")

    full <- data.frame(bin = grid_centers, time = grid_labels, stringsAsFactors = FALSE) |>
      dplyr::left_join(agg, by = c("bin" = "bin"))

    data.frame(time = full$time, value = full$value, stringsAsFactors = FALSE)
  }

  # ---------------------------------- checks ----------------------------------
  if (!dir.exists(inputdir)) stop("inputdir does not exist: ", inputdir)
  if (!is.numeric(grid_mins) || length(grid_mins) != 1L || grid_mins <= 0) {
    stop("`grid_mins` must be a positive numeric scalar (minutes).")
  }
  if (1440 %% grid_mins != 0) {
    stop("`grid_mins` must divide 1440 exactly (e.g., 1,2,3,4,5,6,8,10,12,15,20,30,60).")
  }
  if (!is.logical(keep_cov_order) || length(keep_cov_order) != 1L) stop("`keep_cov_order` must be TRUE/FALSE.")
  if (!is.logical(make_famm) || length(make_famm) != 1L) stop("`make_famm` must be TRUE/FALSE.")
  if (!is.numeric(smooth_k) || length(smooth_k) != 1L || smooth_k <= 0 || smooth_k %% 1 != 0) {
    stop("`smooth_k` must be a positive integer.")
  }
  if (!is.logical(run_pffr) || length(run_pffr) != 1L) stop("`run_pffr` must be TRUE/FALSE.")
  raw_to_grid_strategy <- match.arg(raw_to_grid_strategy)

  if (!is.numeric(upsample_rule) || length(upsample_rule) != 1L || !(upsample_rule %in% c(1,2))) {
    stop("`upsample_rule` must be 1 (no extrapolation) or 2 (extrapolate).")
  }

  if (!is.null(coef_grid_mins)) {
    if (!is.numeric(coef_grid_mins) || length(coef_grid_mins) != 1L || coef_grid_mins <= 0) {
      stop("`coef_grid_mins` must be NULL or a positive numeric scalar (minutes).")
    }
    if (1440 %% coef_grid_mins != 0) stop("`coef_grid_mins` must divide 1440 exactly.")
  }

  file_paths <- list.files(inputdir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(file_paths)) stop("No .csv files found in inputdir: ", inputdir)

  if (!is.character(outputdir) || length(outputdir) != 1L || !nzchar(outputdir)) {
    stop("`outputdir` must be a non-empty character scalar (a directory path).")
  }
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
  }
  outputdir <- normalizePath(outputdir, winslash = "/", mustWork = TRUE)

  # -------------------- read & aggregate per subject (robust grid) --------------------
  all_rows <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    dat <- utils::read.csv(file_path, stringsAsFactors = FALSE)

    if (!(cgm_id_col %in% names(dat))) {
      names(dat)[1L] <- cgm_id_col
    }
    if (!("timestamp" %in% names(dat))) stop("File ", basename(file_path), " does not contain 'timestamp'.")
    if (!("sensorglucose" %in% names(dat))) stop("File ", basename(file_path), " does not contain 'sensorglucose'.")

    ts_parsed <- suppressWarnings(
      lubridate::parse_date_time(
        as.character(dat$timestamp),
        orders = c("ymd HMS","ymd HM","ymd","Ymd HMS","Ymd HM","Ymd",
                   "mdy HMS","mdy HM","mdy",
                   "dmy HMS","dmy HM","dmy"),
        tz = tz
      )
    )
    if (all(is.na(ts_parsed))) stop("Unable to parse 'timestamp' in file: ", basename(file_path))
    dat$timestamp <- ts_parsed

    id_vec <- unique(dat[[cgm_id_col]])
    if (!length(id_vec)) {
      warning("File ", basename(file_path), " has no id; skipping.")
      next
    }
    if (length(id_vec) > 1L) warning("File ", basename(file_path), " contains multiple IDs; using the first.")
    id_first <- as.character(id_vec[1L])

    g <- .aggregate_subject_to_grid(
      df = dat,
      tz = tz,
      grid_mins = grid_mins,
      raw_to_grid_strategy = raw_to_grid_strategy,
      upsample_rule = upsample_rule
    )
    if (is.null(g)) {
      warning("File ", basename(file_path), " has no valid rows after cleaning; skipping.")
      next
    }

    row_wide <- as.data.frame(t(g$value))
    colnames(row_wide) <- as.character(g$time)
    row_wide[[cgm_id_col]] <- id_first
    row_wide <- row_wide[, c(cgm_id_col, setdiff(names(row_wide), cgm_id_col)), drop = FALSE]
    all_rows[[i]] <- row_wide
  }

  all_rows <- all_rows[!vapply(all_rows, is.null, logical(1L))]
  if (!length(all_rows)) stop("No valid rows could be constructed from inputdir: ", inputdir)

  cgm_wide_hhmm <- dplyr::bind_rows(all_rows)

  # -------------------- HH:MM fixed grid -> minutes matrix (complete, consistent) --------------------
  numeric_cols_hhmm <- setdiff(names(cgm_wide_hhmm), cgm_id_col)
  grid_centers <- .grid_centers_minutes(grid_mins)
  grid_labels  <- .minutes_to_hhmm(grid_centers)

  missing_cols <- setdiff(grid_labels, numeric_cols_hhmm)
  if (length(missing_cols)) {
    for (nm in missing_cols) cgm_wide_hhmm[[nm]] <- NA_real_
    numeric_cols_hhmm <- setdiff(names(cgm_wide_hhmm), cgm_id_col)
  }
  cgm_wide_hhmm <- cgm_wide_hhmm[, c(cgm_id_col, grid_labels), drop = FALSE]

  cgm_matrix_hhmm <- as.matrix(cgm_wide_hhmm[, grid_labels, drop = FALSE])
  rownames(cgm_matrix_hhmm) <- cgm_wide_hhmm[[cgm_id_col]]

  MINS_mat <- cgm_matrix_hhmm
  colnames(MINS_mat) <- as.character(grid_centers)

  # -------------------- covariates --------------------
  covariates <- NULL
  if (!is.null(covariate_file)) {
    covariates <- utils::read.csv(covariate_file, stringsAsFactors = FALSE)
    if (!(cov_id_col %in% names(covariates))) stop("covariate_file must contain id column: ", cov_id_col)
  }

  # -------------------- build cgm_df_out --------------------
  arg_hours_out <- grid_centers / 60

  if (!is.null(covariates)) {
    cgm_df_out <- tibble::as_tibble(covariates)

    if (isTRUE(keep_cov_order)) {
      idx <- match(cgm_df_out[[cov_id_col]], rownames(MINS_mat))
      MIMS_mat_out <- MINS_mat[idx, , drop = FALSE]
    } else {
      common <- intersect(cgm_df_out[[cov_id_col]], rownames(MINS_mat))
      cgm_df_out <- cgm_df_out[match(common, cgm_df_out[[cov_id_col]]), , drop = FALSE]
      MIMS_mat_out <- MINS_mat[match(common, rownames(MINS_mat)), , drop = FALSE]
    }

    cgm_df_out$MIMS_mat <- MIMS_mat_out
    cgm_df_out$MIMS_tf  <- tf::tfd(MIMS_mat_out, arg = arg_hours_out)

    if (cov_id_col != cgm_id_col && !(cgm_id_col %in% names(cgm_df_out))) {
      cgm_df_out[[cgm_id_col]] <- cgm_df_out[[cov_id_col]]
    }
  } else {
    cgm_df_out <- tibble::as_tibble(stats::setNames(
      data.frame(rownames(MINS_mat), stringsAsFactors = FALSE),
      cgm_id_col
    ))
    cgm_df_out$MIMS_mat <- MINS_mat
    cgm_df_out$MIMS_tf  <- tf::tfd(MINS_mat, arg = arg_hours_out)
  }

  # -------------------- cgm_famm_df (hour-scale) --------------------
  cgm_famm_df <- NULL
  MIMS_hour_mat <- NULL

  if (isTRUE(make_famm)) {
    hour_arg <- seq(0.5, 23.5, by = 1)  # centers
    arg_5h <- arg_hours_out

    M5 <- cgm_df_out$MIMS_mat
    M5_smooth <- t(apply(M5, 1, .rollmean_extend, k = smooth_k))

    MIMS_hour_mat <- t(apply(M5_smooth, 1, function(y) {
      stats::approx(
        x = arg_5h,
        y = as.numeric(y),
        xout = hour_arg,
        method = "linear",
        rule = 2
      )$y
    }))

    MIMS_hour_mat <- as.matrix(MIMS_hour_mat)
    colnames(MIMS_hour_mat) <- as.character(hour_arg)

    cgm_famm_df <- cgm_df_out
    cgm_famm_df$MIMS_hour_mat <- I(MIMS_hour_mat)
    hour_wide <- tibble::as_tibble(as.data.frame(MIMS_hour_mat, check.names = FALSE))
    cgm_famm_df <- dplyr::bind_cols(cgm_famm_df, hour_wide)
  }

  # --------------------------------- run pffr ---------------------------------
  if (isTRUE(run_pffr)) {
    if (!isTRUE(make_famm) || is.null(cgm_famm_df) || is.null(MIMS_hour_mat)) {
      stop("run_pffr=TRUE requires make_famm=TRUE (need hour-scale matrix).")
    }
    if (is.null(pffr_group) || !is.character(pffr_group) || length(pffr_group) != 1L) {
      stop("Please provide `pffr_group` as a single column name.")
    }
    if (!is.null(pffr_covars) && (!is.character(pffr_covars) || length(pffr_covars) < 1L)) {
      stop("`pffr_covars` must be NULL or a non-empty character vector.")
    }
    if (!is.null(pffr_id_re) && (!is.character(pffr_id_re) || length(pffr_id_re) != 1L)) {
      stop("`pffr_id_re` must be NULL or a single column name.")
    }

    # ---- optional filtering for pffr subset (single-arg, two modes) ----
    if (!is.null(pffr_filter)) {
      if (!is.list(pffr_filter)) stop("`pffr_filter` must be NULL or a list().")

      drop_na <- if (!is.null(pffr_filter$drop_na)) pffr_filter$drop_na else TRUE
      if (!is.logical(drop_na) || length(drop_na) != 1L) {
        stop("`pffr_filter$drop_na` must be TRUE/FALSE.")
      }

      n_before <- nrow(cgm_famm_df)

      # helper: safely compute keep index
      .eval_filter_expr <- function(expr_obj, data) {
        if (!requireNamespace("rlang", quietly = TRUE)) {
          stop("Package 'rlang' is required when using `pffr_filter$expr`.")
        }

        # expr_obj can be:
        #   - a quoted expression: quote(a > 0)
        #   - rlang expr/quosure: rlang::expr(a > 0) / quo(a > 0)
        #   - a value that user wrote without quoting (e.g., a > 5),
        #     which might have been pre-evaluated or errored elsewhere; we still try to recover.
        if (rlang::is_quosure(expr_obj)) {
          q <- expr_obj
        } else if (rlang::is_call(expr_obj) || rlang::is_symbol(expr_obj) || is.language(expr_obj)) {
          q <- rlang::new_quosure(expr_obj)
        } else {
          # last resort: treat as constant logical / numeric vector
          return(expr_obj)
        }

        rlang::eval_tidy(q, data = data)
      }

      # Mode 1: expr = <expression evaluated in cgm_famm_df>
      if (!is.null(pffr_filter$expr)) {

        expr <- pffr_filter$expr

        # If user accidentally passed an unquoted expression via list(),
        # expr might already be evaluated outside; to be robust, try to capture it.
        # Example: list(expr = a > 5)  (bad user-side)
        # We attempt to recover the original expression from the call.
        if (!is.language(expr) && !requireNamespace("rlang", quietly = TRUE)) {
          # rlang missing handled inside helper
        }

        keep_idx <- NULL
        keep_idx <- tryCatch(
          .eval_filter_expr(expr, data = cgm_famm_df),
          error = function(e) {
            # Recover from common misuse: list(expr = a > 5)
            # by capturing the unevaluated expression from pffr_filter call
            expr_sub <- tryCatch(substitute(pffr_filter$expr), error = function(e2) NULL)
            if (!is.null(expr_sub) && (is.language(expr_sub) || is.symbol(expr_sub))) {
              .eval_filter_expr(expr_sub, data = cgm_famm_df)
            } else {
              stop("Failed to evaluate `pffr_filter$expr`: ", conditionMessage(e))
            }
          }
        )

        # allow constant logical vector returned directly
        if (is.numeric(keep_idx)) keep_idx <- keep_idx != 0

        if (!is.logical(keep_idx) || length(keep_idx) != n_before) {
          stop("`pffr_filter$expr` must evaluate to a logical vector of length nrow(data).")
        }
        if (isTRUE(drop_na)) keep_idx[is.na(keep_idx)] <- FALSE

      } else {
        # Mode 2: var + keep (membership filter)
        var <- pffr_filter$var
        keep <- pffr_filter$keep

        if (!is.character(var) || length(var) != 1L || !nzchar(var)) {
          stop("Provide either `pffr_filter$expr` or (`pffr_filter$var` + `pffr_filter$keep`).")
        }
        if (!(var %in% names(cgm_famm_df))) stop("`pffr_filter$var` not found: ", var)
        if (is.null(keep) || length(keep) < 1L) stop("`pffr_filter$keep` must be non-empty.")

        v <- cgm_famm_df[[var]]
        keep_idx <- v %in% keep
        if (isTRUE(drop_na)) keep_idx <- keep_idx & !is.na(v)
      }

      # apply filter + keep matrices consistent
      keep_idx <- as.logical(keep_idx)
      if (isTRUE(drop_na)) keep_idx[is.na(keep_idx)] <- FALSE

      cgm_famm_df <- cgm_famm_df[keep_idx, , drop = FALSE]

      # IMPORTANT: keep MIMS_hour_mat consistent with filtered rows
      if (!is.null(MIMS_hour_mat)) {
        MIMS_hour_mat <- MIMS_hour_mat[keep_idx, , drop = FALSE]
        cgm_famm_df$MIMS_hour_mat <- I(MIMS_hour_mat)
      }

      n_after <- nrow(cgm_famm_df)
      if (n_after < 2L) {
        stop("pffr_filter left too few rows: ", n_after, " (before: ", n_before, ").")
      }

      if (isTRUE(verbose)) {
        base::message("pffr_filter applied: n=", n_before, " -> ", n_after)
      }
    }

    needed <- c(pffr_group, "MIMS_hour_mat")
    if (!is.null(pffr_covars)) needed <- c(needed, pffr_covars)
    if (!is.null(pffr_id_re))  needed <- c(needed, pffr_id_re)

    missing_cols <- setdiff(needed, names(cgm_famm_df))
    if (length(missing_cols)) stop("Missing columns in cgm_famm_df: ", paste(missing_cols, collapse = ", "))

    if (!requireNamespace("refund", quietly = TRUE)) stop("Package 'refund' is required for run_pffr=TRUE.")
    if (isTRUE(make_plot) && !requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for make_plot=TRUE.")

    if (!is.null(pffr_id_re)) {
      cgm_famm_df[[pffr_id_re]] <- as.factor(cgm_famm_df[[pffr_id_re]])
    }

    g0 <- cgm_famm_df[[pffr_group]]

    # ---- FORCE group into TRUE numeric 0/1 (and fix common 1/2 coding) ----
    .coerce_binary01 <- function(x) {
      x_raw <- x

      if (is.factor(x)) x <- as.character(x)
      if (is.logical(x)) return(ifelse(is.na(x), NA_real_, as.numeric(x)))

      if (is.character(x)) {
        xx <- trimws(tolower(x))
        xx[xx %in% c("false","f","no","n")]  <- "0"
        xx[xx %in% c("true","t","yes","y")] <- "1"
        x <- suppressWarnings(as.numeric(xx))
      }

      x <- suppressWarnings(as.numeric(x))
      u <- sort(unique(x[!is.na(x)]))

      if (!length(u)) {
        u0 <- sort(unique(stats::na.omit(as.character(x_raw))))
        stop("`pffr_group` cannot be coerced to numeric. Values are: ",
             paste(utils::head(u0, 20), collapse = ", "),
             if (length(u0) > 20) " ..." else "",
             ". Please recode to 0/1 (or 1/2).")
      }

      if (length(u) <= 2 && all(u %in% c(0, 1))) return(x)
      if (length(u) <= 2 && all(u %in% c(1, 2))) return(x - 1)

      stop("`pffr_group` must be binary coded as 0/1 (or 1/2). Got: ",
           paste(u, collapse = ", "),
           ". Please recode your group before running.")
    }

    cgm_famm_df[["pffr_group_internal"]] <- .coerce_binary01(g0)

    u_int <- sort(unique(cgm_famm_df[["pffr_group_internal"]][!is.na(cgm_famm_df[["pffr_group_internal"]])]))
    if (!all(u_int %in% c(0,1))) stop("Internal group coding failed: ", paste(u_int, collapse=", "))

    cgm_famm_df$MIMS_hour_mat <- I(cgm_famm_df$MIMS_hour_mat)
    nyindex <- ncol(cgm_famm_df$MIMS_hour_mat)

    yind <- pffr_yind
    if (is.null(yind)) {
      yind <- seq(0.5, 23.5, length.out = nyindex)
    } else {
      if (!is.numeric(yind) || length(yind) != nyindex) {
        stop("`pffr_yind` length must equal ncol(MIMS_hour_mat) = ", nyindex)
      }
    }

    rhs_terms <- c("pffr_group_internal")
    if (!is.null(pffr_covars) && length(pffr_covars)) rhs_terms <- c(rhs_terms, pffr_covars)
    if (!is.null(pffr_id_re) && nzchar(pffr_id_re)) rhs_terms <- c(rhs_terms, sprintf("s(%s, bs = \"re\")", pffr_id_re))

    fml <- stats::as.formula(paste("MIMS_hour_mat ~", paste(rhs_terms, collapse = " + ")))

    pffr_fit <- refund::pffr(
      formula   = fml,
      yind      = yind,
      data      = cgm_famm_df,
      algorithm = pffr_algorithm,
      discrete  = pffr_discrete,
      bs.yindex = pffr_bs_yindex
    )

    # ---- save pffr summary to txt (optional; default follows run_pffr) ----
    if (isTRUE(pffr_summary_save)) {

      sum_file <- pffr_summary_file

      if (.is_bare_filename(sum_file)) {
        ext  <- tools::file_ext(sum_file)
        base <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", sum_file) else sum_file
        base <- paste0(base, "_", pffr_group)

        sum_file <- if (nzchar(ext)) {
          file.path(outputdir, paste0(base, ".", ext))
        } else {
          file.path(outputdir, base)
        }
      }

      txt <- utils::capture.output(summary(pffr_fit))
      writeLines(txt, con = sum_file, useBytes = TRUE)
      pffr_summary_file_written <- sum_file
    }

    pffr_call <- list(
      formula    = fml,
      yind       = yind,
      group      = pffr_group,
      covars     = pffr_covars,
      id_re      = pffr_id_re,
      algorithm  = pffr_algorithm,
      discrete   = pffr_discrete,
      bs.yindex  = pffr_bs_yindex,
      coef_grid_mins = coef_grid_mins,
      plot_file = plot_file,
      coef_save = coef_save,
      raw_to_grid_strategy = raw_to_grid_strategy,
      upsample_rule = upsample_rule,
      outputdir = outputdir
    )

    # ---------------- make_plot block ----------------
    if (isTRUE(make_plot)) {

      group_tag <- if (!is.null(pffr_group) && nzchar(pffr_group)) pffr_group else "group"

      # if user did not customize plot_file and kept default name, swap to "ALL" name
      if (.is_bare_filename(plot_file) && identical(plot_file, "pffr_coef_plot.pdf")) {
        plot_file <- .default_plot_file_all(group_tag)
      }

      coef_file <- .default_coef_file(group_tag)
      coef_file_grid <- if (!is.null(coef_grid_mins)) .default_coef_file_grid(group_tag, coef_grid_mins) else NULL

      # ---- base coef df + base plot ----
      pffr_coef_df <- .termcurve_from_pffr_smterms(
        pffr_fit         = pffr_fit,
        terms            = plot_terms,
        drop_intercept   = plot_drop_intercept,
        drop_re          = plot_drop_re,
        drop_yind_smooth = plot_drop_yind_smooth,
        clean_term_label = TRUE
      )
      pffr_coef_df <- .reorder_coef_df_raw(pffr_coef_df)
      pffr_plot <- .plot_termcurve_df(
        pffr_coef_df,
        title = "Adjusted time-varying coefficient functions (original grid)"
      )

      # ---- grid df + grid plot (soft-fail) ----
      pffr_coef_df_grid <- NULL
      pffr_plot_grid <- NULL
      grid_ok <- NA
      grid_reason <- NA_character_

      if (!is.null(coef_grid_mins)) {
        grid_try <- tryCatch(
          {
            df_grid <- .upsample_termcurve_grid(pffr_coef_df, grid_mins = coef_grid_mins)
            need_cols <- c("term", "hour", "time", "coef", "lb", "ub")
            if (is.null(df_grid) || !is.data.frame(df_grid) || nrow(df_grid) == 0L) stop("grid coef df is empty.")
            if (!all(need_cols %in% names(df_grid))) {
              stop("grid coef df missing columns: ",
                   paste(setdiff(need_cols, names(df_grid)), collapse = ", "))
            }
            plt_grid <- .plot_termcurve_df(
              df_grid,
              title = paste0("Adjusted time-varying coefficient functions (", coef_grid_mins, "-min grid)")
            )
            list(ok = TRUE, df = df_grid, plot = plt_grid, reason = NA_character_)
          },
          error = function(e) list(ok = FALSE, df = NULL, plot = NULL, reason = conditionMessage(e))
        )

        grid_ok <- isTRUE(grid_try$ok)
        grid_reason <- grid_try$reason

        if (grid_ok) {
          pffr_coef_df_grid <- .reorder_coef_df_grid(grid_try$df)
          pffr_plot_grid <- grid_try$plot
        } else {
          pffr_coef_df_grid <- NULL
          pffr_plot_grid <- NULL
          coef_file_grid <- NULL
        }
      }

      # ---- save plots (DEFAULT: when run_pffr=TRUE, plot_save follows make_plot => TRUE) ----
      if (isTRUE(plot_save)) {

        if (.is_bare_filename(plot_file)) plot_file <- file.path(outputdir, plot_file)

        if (!is.null(plot_terms) && length(plot_terms) >= 1L) {
          ggplot2::ggsave(plot_file, pffr_plot, width = 10, height = 5, dpi = 300)

          if (!is.null(pffr_plot_grid)) {
            plot_file_grid_local <- .derive_plot_file_grid(plot_file, coef_grid_mins)
            ggplot2::ggsave(plot_file_grid_local, pffr_plot_grid, width = 10, height = 5, dpi = 300)
            plot_file_grid_written <- TRUE
            if (!is.null(pffr_call)) pffr_call$plot_file_grid <- plot_file_grid_local
          }

        } else {

          # DEFAULT MODE:
          #   - write ALL terms (original + optional grid) into ONE pdf
          #   - write GROUP term (original + optional grid) into ANOTHER pdf

          # 2.1 ALL TERMS
          file_all <- plot_file
          if (.is_bare_filename(file_all) && identical(file_all, "pffr_coef_plot.pdf")) {
            file_all <- .default_plot_file_all(group_tag)
          }
          if (.is_bare_filename(file_all)) file_all <- file.path(outputdir, file_all)

          .save_two_pages_pdf(
            file  = file_all,
            plot1 = pffr_plot,
            plot2 = pffr_plot_grid
          )

          plot_file <- file_all

          # 2.2 GROUP TERM ONLY
          file_group <- .default_plot_file_group(group_tag)
          if (.is_bare_filename(file_group)) file_group <- file.path(outputdir, file_group)

          nd0 <- cgm_famm_df
          nd1 <- cgm_famm_df
          nd0[["pffr_group_internal"]] <- 0
          nd1[["pffr_group_internal"]] <- 1

          p0 <- stats::predict(pffr_fit, newdata = nd0, type = "response", se.fit = TRUE)
          p1 <- stats::predict(pffr_fit, newdata = nd1, type = "response", se.fit = TRUE)

          fit0 <- p0$fit; fit1 <- p1$fit
          se0  <- p0$se.fit; se1 <- p1$se.fit

          nyindex2 <- length(pffr_call$yind)
          yind2 <- pffr_call$yind

          if (!is.null(dim(fit0)) && ncol(fit0) == nyindex2) {
            diff_fit <- colMeans(fit1 - fit0, na.rm = TRUE)
            diff_se  <- sqrt(colMeans(se1^2 + se0^2, na.rm = TRUE))
          } else if (!is.null(dim(fit0)) && nrow(fit0) == nyindex2) {
            diff_fit <- rowMeans(fit1 - fit0, na.rm = TRUE)
            diff_se  <- sqrt(rowMeans(se1^2 + se0^2, na.rm = TRUE))
          } else {
            stop("Unexpected predict() output shape: ", paste(dim(fit0), collapse="x"))
          }

          df_group_raw <- data.frame(
            term = "group",
            yind = yind2,
            hour = as.numeric(yind2),
            time = .hour_to_hhmm(as.numeric(yind2)),
            coef = as.numeric(diff_fit),
            se   = as.numeric(diff_se),
            lb   = as.numeric(diff_fit - 1.96 * diff_se),
            ub   = as.numeric(diff_fit + 1.96 * diff_se),
            stringsAsFactors = FALSE
          )
          df_group_raw <- .reorder_coef_df_raw(df_group_raw)

          plt_group <- .plot_termcurve_df(
            df_group_raw,
            title = paste0("Adjusted time-varying effect of ", pffr_group, " (1 vs 0; original grid)")
          )

          plt_group_grid <- NULL
          if (!is.null(coef_grid_mins)) {
            df_g_grid <- .upsample_termcurve_grid(df_group_raw, grid_mins = coef_grid_mins)
            plt_group_grid <- .plot_termcurve_df(
              df_g_grid,
              title = paste0("Adjusted time-varying effect of ", pffr_group,
                             " (1 vs 0; ", coef_grid_mins, "-min grid)")
            )
          }

          .save_two_pages_pdf(
            file  = file_group,
            plot1 = plt_group,
            plot2 = plt_group_grid
          )

          if (!is.null(pffr_call)) {
            pffr_call$plot_file_all   <- file_all
            pffr_call$plot_file_group <- file_group
            pffr_call$plot_file       <- file_all
          }
        }
      }

      # ---- save coef tables (DEFAULT: follows make_plot => TRUE) ----
      if (isTRUE(coef_save)) {

        if (.is_bare_filename(coef_file)) coef_file <- file.path(outputdir, coef_file)
        utils::write.csv(pffr_coef_df, coef_file, row.names = FALSE)

        if (!is.null(pffr_coef_df_grid) && !is.null(coef_file_grid)) {
          if (.is_bare_filename(coef_file_grid)) coef_file_grid <- file.path(outputdir, coef_file_grid)
          utils::write.csv(pffr_coef_df_grid, coef_file_grid, row.names = FALSE)
          coef_file_grid_written <- TRUE
        }
      }

      # ---- write back to pffr_call ----
      if (!is.null(pffr_call)) {
        pffr_call$outputdir       <- outputdir
        pffr_call$plot_file       <- plot_file
        pffr_call$coef_file       <- coef_file
        pffr_call$coef_file_grid  <- coef_file_grid
        pffr_call$grid_ok         <- grid_ok
        pffr_call$grid_reason     <- grid_reason
        pffr_call$plot_file_grid_written <- plot_file_grid_written
        pffr_call$coef_file_grid_written <- coef_file_grid_written
      }
    } # end make_plot
  } # end run_pffr

  # -------------------- return (ALWAYS return, regardless of run_pffr) --------------------
  out <- list(
    cgm_df            = cgm_df_out,
    cgm_famm_df       = cgm_famm_df,

    pffr_fit          = pffr_fit,
    pffr_call         = pffr_call,

    pffr_coef_df      = pffr_coef_df,
    pffr_coef_df_grid = pffr_coef_df_grid,

    pffr_plot         = pffr_plot,
    pffr_plot_grid    = pffr_plot_grid,

    coef_grid_mins     = coef_grid_mins,
    plot_file          = plot_file,
    coef_save          = coef_save,
    coef_file          = coef_file,
    coef_file_grid     = coef_file_grid,
    pffr_summary_file  = pffr_summary_file_written,
    outputdir          = outputdir,

    raw_to_grid_strategy = raw_to_grid_strategy,
    upsample_rule = upsample_rule,

    grid_ok = grid_ok,
    grid_reason = grid_reason,
    plot_file_grid_written = plot_file_grid_written,
    coef_file_grid_written = coef_file_grid_written
  )

  if (isTRUE(verbose) &&
      (isTRUE(plot_save) || isTRUE(coef_save) || isTRUE(pffr_summary_save))) {

    base::message("Outputs saved under: ", outputdir)

    if (isTRUE(plot_save) && !is.null(plot_file) && nzchar(plot_file))
      base::message("Plot (all terms): ", plot_file)

    if (!is.null(pffr_call$plot_file_group) && nzchar(pffr_call$plot_file_group))
      base::message("Plot (group): ", pffr_call$plot_file_group)

    if (isTRUE(plot_file_grid_written) &&
        !is.null(pffr_call$plot_file_grid) && nzchar(pffr_call$plot_file_grid)) {
      base::message("Plot (grid): ", pffr_call$plot_file_grid)
    }

    if (isTRUE(coef_save) && !is.null(coef_file) && nzchar(coef_file))
      base::message("Coef: ", coef_file)

    if (isTRUE(coef_file_grid_written) && !is.null(coef_file_grid) && nzchar(coef_file_grid))
      base::message("Coef (grid): ", coef_file_grid)

    if (!is.null(pffr_summary_file_written) && nzchar(pffr_summary_file_written))
      base::message("pffr summary: ", pffr_summary_file_written)
  }

  return(out)
}
