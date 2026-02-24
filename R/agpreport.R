if (requireNamespace("utils", quietly = TRUE)) {
  utils::globalVariables(c("xmin", "xmax", "ymin", "ymax"))
}

# ============================= Generate AGP Report ============================
#' Generate Ambulatory Glucose Profile Report
#'
#' @description
#' This function takes a directory of cleaned CGM data and generates aggregate AGPs
#' for all of the files combined, using either raw bin-wise quantiles, Tukey 3R
#' smoothing, or a LOESS-smoothed median curve. It also produces an optional set of
#' per-subject AGP reports (each subject in a separate PDF), including:
#' \enumerate{
#'   \item subject-level AGP (original / Tukey / LOESS),
#'   \item multi-day glucose profile overlays, and
#'   \item daily glucose profile cards.
#' }
#'
#' Both the subject-level 24-hour overlay plot (labeled
#' \dQuote{Subject-level 24-hour Glucose Profiles}) and the multi-day glucose overlay
#' plot follow the value of \code{agg_smooth}. Specifically, when
#' \code{agg_smooth = "original"}, curves are drawn from raw observations without
#' smoothing; when \code{agg_smooth = "tukey"}, curves are smoothed using Tukey 3R
#' (applied separately within each subject or each day, respectively); and when
#' \code{agg_smooth = "loess"}, curves are LOESS-smoothed accordingly.
#'
#' @param inputdir The directory containing all cleaned, formatted CGM data to be analyzed.
#'
#' @param outputdir The directory where plot PDF files should be written.
#'
#' @param tz The time zone in which the data were recorded.
#'
#' @param yaxis Numeric vector of length 2 giving the y-axis limits (in mg/dL) used for all plots.
#'
#' @param agg_smooth Character; aggregate AGP smoothing method. One of
#'   \code{"original"}, \code{"tukey"}, or \code{"loess"}:
#'   \itemize{
#'     \item \code{"original"}: no smoothing, use raw bin-wise quantiles based on the native sampling interval.
#'     \item \code{"tukey"}: apply Tukey 3R smoothing to the bin-wise quantiles.
#'     \item \code{"loess"}: LOESS-smooth the aggregate median curve.
#'   }
#'
#' @param ptype Character; participant type used to choose default horizontal reference lines.
#'   One of: \code{"dm"}, \code{"gdm"}, \code{"normal"}, or \code{"other"}.
#'
#' @param hlines Numeric vector; custom horizontal reference lines (mg/dL).
#'   If provided, these values are used as the horizontal reference lines. If \code{NULL},
#'   defaults are chosen based on \code{ptype}. When \code{ptype = "other"} and \code{hlines = NULL},
#'   the default is \code{c(70, 180)}.
#'
#' @param mapping_file Optional; path to an external mapping CSV file. If provided, this file is used
#'   to map \code{subjectid} to \code{group_var}. The file must contain at least columns
#'   \code{subjectid} and the specified \code{group_var}. If \code{mapping_file = NULL}, the function will
#'   look for the internal \code{inst/exdata/mapping.csv} shipped with the package. If neither is available,
#'   grouping/filtering via mapping is skipped.
#'
#' @param filter_var Optional; filtering rule evaluated on columns in the mapping file.
#'   Supports a lightweight expression syntax with parentheses and AND/OR:
#'   \itemize{
#'     \item AND: \code{&}, OR: \code{|}, parentheses: \code{( ... )}
#'     \item Comparisons: \code{==, !=, >, >=, <, <=}
#'     \item Ranges: \code{var between a and b} or \code{var in a:b}
#'     \item Membership: \code{site==1,2,3} (comma-separated values)
#'   }
#'   Example: \code{"male==1 & (age>=60 | egfr<90)"}.
#'
#' @param group_var Optional; column name in the mapping file used to group subjects in aggregate AGP plots.
#'   If \code{NULL} (default), no grouping is applied. The mapping file must contain \code{subjectid} and the
#'   supplied \code{group_var}.
#'
#' @param bin_mins Numeric; width of the time-of-day bins (in minutes) used to compute within-bin quantiles
#'   for the aggregate AGP. Defaults to 5 minutes.
#'
#' @param per_subject Logical; whether to generate individual-level AGP PDF reports for each subject.
#'   \itemize{
#'     \item \code{FALSE} (default): only the population-level AGP report \code{AGP_Report_All.pdf} is produced.
#'     \item \code{TRUE}: for each subject, a separate PDF is saved containing:
#'       \enumerate{
#'         \item the subject-specific AGP (according to \code{agg_smooth}),
#'         \item multi-day glucose overlays, and
#'         \item daily glucose profile cards with zone-colored segments.
#'       }
#'   }
#'
#' @param subject_overlay Logical; whether to compute and draw the subject-level 24-hour overlay plot
#'   (scatter + per-subject smoothed curves). Defaults to \code{FALSE} because it can be computationally
#'   expensive when the sample size is large.
#'
#' @param page_width Width of output PDF (in inches).
#' @param page_height Height of output PDF (in inches).
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item \code{data}: cleaned and combined CGM records.
#'     \item \code{quartiles}: aggregate AGP quantiles (bin-wise).
#'     \item \code{quartiles_grouped}: grouped quantiles (if grouping used).
#'     \item \code{quartiles_grouped_original}: grouped quantiles for the original (native-interval) AGP.
#'     \item \code{aggAGP_original}, \code{aggAGPtukey}, \code{aggAGPloess}: aggregate AGPs.
#'     \item \code{aggAGP_original_grouped}, \code{aggAGPtukey_grouped}, \code{aggAGPloess_grouped}:
#'       grouped aggregate AGPs (if grouping used).
#'     \item \code{AGPloess}: subject-level 24-hour glucose profile overlay
#'       (only when \code{subject_overlay = TRUE}; otherwise \code{NULL}). The smoothing method follows
#'       \code{agg_smooth}: raw (no smoothing), Tukey 3R, or LOESS smoothing.
#'   }
#'
#' @examples
#' \donttest{
#' ## Minimal runnable example
#' agpreport(
#'   inputdir  = system.file("exdata", "tidydata", package = "cgmcalculator"),
#'   outputdir = tempdir()
#' )
#' }
#'
#' \dontrun{
#' ## 2) LOESS smoothing
#' agpreport(
#'   inputdir   = "path/to/tidydata",
#'   outputdir  = "results/",
#'   agg_smooth = "loess"
#' )
#'
#' ## 3) Per-subject PDFs
#' agpreport(
#'   inputdir    = "path/to/tidydata",
#'   outputdir   = "results/",
#'   per_subject = TRUE
#' )
#'
#' ## 4) Group-specific AGP
#' agpreport(
#'   inputdir  = "path/to/tidydata",
#'   outputdir = "results/",
#'   group_var = "group"
#' )
#' }
#'
#' @importFrom rlang .data
#' @importFrom grid grid.draw
#' @importFrom ggplot2 ggplotGrob
#' @importFrom tidyr complete
#' @seealso viridis::scale_color_viridis_d
#' @seealso viridis::scale_fill_viridis_d
#' @export
agpreport <- function(inputdir,
                      outputdir = tempdir(),
                      tz = "UTC",
                      yaxis = c(0, 400),
                      agg_smooth = c("original", "loess", "tukey"),
                      ptype = c("dm", "gdm", "normal", "other"),
                      hlines = NULL,
                      mapping_file = NULL,
                      group_var = NULL,
                      filter_var = NULL,
                      bin_mins = 5,
                      per_subject = FALSE,
                      subject_overlay = FALSE,
                      page_width = 11,
                      page_height = 8.5) {

  agg_smooth <- match.arg(agg_smooth)
  ptype      <- match.arg(ptype)

  # ------------------ 1. Basic checks ---------------------------------------
  if (!dir.exists(inputdir)) {
    stop("Input directory '", inputdir, "' does not exist.", call. = FALSE)
  }

  files <- base::list.files(path = inputdir, full.names = TRUE, pattern = "\\.[Cc][Ss][Vv]$")
  if (length(files) == 0L) {
    stop("No files found in 'inputdir'. Please check the directory path.",
         call. = FALSE)
  }

  if (!is.numeric(yaxis) || length(yaxis) != 2L || yaxis[1] >= yaxis[2]) {
    stop("'yaxis' must be a numeric vector of length 2 with yaxis[1] < yaxis[2].",
         call. = FALSE)
  }

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- bin_mins: validate + coerce to integer minutes (safer for lubridate::round_date) ----
  if (length(bin_mins) != 1L) {
    stop("'bin_mins' must be a positive numeric scalar.", call. = FALSE)
  }
  bin_mins <- as.numeric(bin_mins)
  if (!is.finite(bin_mins) || bin_mins <= 0) {
    stop("'bin_mins' must be a positive numeric scalar.", call. = FALSE)
  }

  bin_mins_int <- as.integer(round(bin_mins))
  if (is.na(bin_mins_int) || bin_mins_int <= 0) {
    stop("'bin_mins' must be a positive numeric scalar.", call. = FALSE)
  }
  if (bin_mins_int != bin_mins) {
    warning("'bin_mins' was not an integer; using rounded value: ", bin_mins_int)
  }
  bin_mins <- bin_mins_int

  if ((1440 %% bin_mins) != 0L) {
    warning(
      "'bin_mins' (", bin_mins, ") does not divide evenly into 1440.\n",
      "The last bin may be slightly shorter or longer."
    )
  }

  ## ---- 1.1 Decide reference lines from ptype / hlines ----------------------
  if (is.null(hlines)) {
    ref_lines <- switch(
      ptype,
      "dm"     = c(70, 180),
      "gdm"    = c(63, 140),
      "normal" = c(70, 140),
      "other"  = c(70, 180)  # default if ptype = "other" and hlines not specified
    )
  } else {
    ref_lines <- as.numeric(hlines)
    ref_lines <- ref_lines[is.finite(ref_lines)]
  }

  # Infer low / high thresholds from the reference lines (for example 70 / 180)
  if (length(ref_lines) >= 2L) {
    low_thr  <- min(ref_lines)
    high_thr <- max(ref_lines)
  } else if (length(ref_lines) == 1L) {
    # Single reference line: treat it as upper threshold, use a default lower threshold
    low_thr  <- 70
    high_thr <- ref_lines[1]
  } else {
    # No reference lines at all: fall back to 70-180
    low_thr  <- 70
    high_thr <- 180
  }

  # ------------------ 2. Read and combine all files -------------------------
  datalist <- lapply(files, function(f) {
    cgmdata <- utils::read.csv(
      f,
      stringsAsFactors = FALSE,
      header = TRUE,
      skipNul = TRUE
    )

    required_cols <- c("subjectid", "timestamp", "sensorglucose")
    missing_cols  <- setdiff(required_cols, names(cgmdata))
    if (length(missing_cols) > 0L) {
      stop("File '", basename(f), "' is missing required column(s): ",
           paste(missing_cols, collapse = ", "),
           call. = FALSE)
    }

    # Ensure subjectid is constant per file (same behaviour as original code)
    id <- cgmdata$subjectid[1]
    if (any(cgmdata$subjectid != id, na.rm = TRUE)) {
      warning("File '", basename(f), "' contains multiple subjectid values; using the first one: ", id)
    }
    cgmdata$subjectid <- id

    cgmdata
  })

  aggregateAGPdata <- do.call(rbind, datalist)

  # ------------------ 3. Clean numeric + timestamps -------------------------
  # Convert glucose to numeric (Libre / Dexcom "Low"/"High" should ideally be handled in tidy step)
  aggregateAGPdata$sensorglucose <- suppressWarnings(
    as.numeric(aggregateAGPdata$sensorglucose)
  )

  # Parse timestamp into the recording timezone
  # First coerce to character to avoid weird internal types
  ts_chr <- as.character(aggregateAGPdata$timestamp)

  # Treat empty strings as NA
  ts_chr[!nzchar(ts_chr)] <- NA_character_

  # Use lubridate::parse_date_time to flexibly parse multiple formats
  aggregateAGPdata$timestamp <- suppressWarnings(
    lubridate::parse_date_time(
      ts_chr,
      orders = c(
        "ymd HMS",
        "ymd HM",
        "ymd",
        "Ymd HMS",
        "Ymd HM",
        "Ymd",
        "mdy HMS",
        "mdy HM",
        "mdy",
        "dmy HMS",
        "dmy HM",
        "dmy"
      ),
      tz = tz
    )
  )

  # Do not error on unparsed timestamps; just warn with examples
  bad_idx <- is.na(aggregateAGPdata$timestamp) & !is.na(ts_chr)
  if (any(bad_idx)) {
    warning(
      "Some timestamps could not be parsed. Examples: ",
      paste(utils::head(unique(ts_chr[bad_idx]), 5), collapse = ", ")
    )
  }

  # Keep only rows with complete key fields
  keep <- stats::complete.cases(
    aggregateAGPdata[, c("subjectid", "timestamp", "sensorglucose")]
  )
  aggregateAGPdata <- aggregateAGPdata[keep, ]

  if (!nrow(aggregateAGPdata)) {
    stop("No valid CGM records after cleaning 'subjectid', 'timestamp', and ",
         "'sensorglucose'. Please check the input files.", call. = FALSE)
  }

  # ------------------ 4. Time-of-day (minutes) ------------------------------
  # Work purely in "minutes since midnight" to avoid POSIXct / timezone quirks.
  ts_local <- aggregateAGPdata$timestamp

  # Continuous time-of-day in minutes (for LOESS overlays)
  tod_minutes <- lubridate::hour(ts_local) * 60 +
    lubridate::minute(ts_local) +
    lubridate::second(ts_local) / 60

  # Rounded bins in X minutes (for quantiles / Tukey AGP)
  ts_bin <- lubridate::round_date(ts_local, unit = paste(bin_mins, "mins"))
  hourmin_minutes <- lubridate::hour(ts_bin) * 60 +
    lubridate::minute(ts_bin)

  aggregateAGPdata$time_minutes    <- tod_minutes
  aggregateAGPdata$hourmin_minutes <- hourmin_minutes

  # Drop rows where time-of-day construction failed (extremely rare)
  aggregateAGPdata <- aggregateAGPdata[
    stats::complete.cases(aggregateAGPdata[, c("time_minutes", "hourmin_minutes")]),
    ,
    drop = FALSE
  ]

  ## ---- 4.1 Native sampling interval (for all "original" AGPs) --------------
  dt_all <- diff(sort(unique(aggregateAGPdata$timestamp)))
  dt_mins <- as.numeric(dt_all, units = "mins")
  native_mins <- if (length(dt_mins)) stats::median(dt_mins[is.finite(dt_mins) & dt_mins > 0]) else NA_real_

  # Fallback: if native interval cannot be determined, revert to bin_mins
  if (is.na(native_mins) || !is.finite(native_mins) || native_mins <= 0) {
    native_mins <- bin_mins
  }

  bin_mins_original <- max(1L, round(native_mins))

  ts_bin_orig <- lubridate::round_date(
    ts_local,
    unit = paste(bin_mins_original, "mins")
  )
  hourmin_minutes_orig <- lubridate::hour(ts_bin_orig) * 60 +
    lubridate::minute(ts_bin_orig)

  aggregateAGPdata$hourmin_minutes_orig <- hourmin_minutes_orig

  if (!nrow(aggregateAGPdata)) {
    stop("No valid records after constructing time-of-day variables. ",
         "Check timestamp formats and time zone.", call. = FALSE)
  }

  ## ------------------ 4.2 Optional: mapping (grouping and/or filtering) -------
  need_map <- !is.null(group_var) || !is.null(filter_var)

  if (isTRUE(need_map)) {

    # ---- helper: read + normalize mapping ----
    .read_mapping <- function(path, group_var = NULL, filter_var = NULL) {

      df <- utils::read.csv(
        path,
        stringsAsFactors = FALSE,
        check.names      = FALSE,
        fileEncoding     = "UTF-8-BOM"
      )

      # 1) normalize column names: trim/BOM/lower
      nm <- names(df)
      nm <- trimws(nm)
      nm <- gsub("\ufeff", "", nm, fixed = TRUE)
      nm <- tolower(nm)
      names(df) <- nm

      # 2) normalize requested variable names to match
      if (!is.null(group_var))  group_var  <- tolower(trimws(group_var))
      if (!is.null(filter_var)) filter_var <- tolower(trimws(filter_var))

      # 3) hard check subjectid
      if (!"subjectid" %in% names(df)) {
        stop(
          "Mapping file must contain a column named 'subjectid'.\n",
          "After cleaning, available columns: ",
          paste(names(df), collapse = ", "),
          call. = FALSE
        )
      }

      # 4) join key (character)
      df$subjectid_chr <- as.character(df$subjectid)

      list(df = df, group_var = group_var, filter_var = filter_var)
    }

    mapping_df   <- NULL
    mapping_path <- NULL

    ### -------------------- 4.2.1 decide mapping path --------------------
    if (!is.null(mapping_file)) {

      if (!file.exists(mapping_file)) {
        warning(
          "A mapping file was provided, but the file does not exist: ",
          mapping_file,
          "\nProceeding without mapping-based grouping/filtering."
        )
        mapping_path <- NULL
      } else {
        mapping_path <- mapping_file
      }

    } else {

      mapping_path0 <- system.file("exdata", "mapping.csv", package = "cgmcalculator")
      if (!nzchar(mapping_path0) || !file.exists(mapping_path0)) {
        warning(
          "Mapping is needed (group_var and/or filter_var specified) but no mapping_file was provided and\n",
          "'mapping.csv' was not found in system.file('exdata','mapping.csv', package='cgmcalculator').\n",
          "Proceeding without mapping-based grouping/filtering."
        )
        mapping_path <- NULL
      } else {
        mapping_path <- mapping_path0
      }
    }

    ### ------------ 4.2.2 read mapping if we have a valid path ------------
    if (!is.null(mapping_path)) {

      tmp <- .read_mapping(mapping_path, group_var = group_var, filter_var = filter_var)
      mapping_df <- tmp$df
      group_var  <- tmp$group_var
      filter_var <- tmp$filter_var

      # join key for CGM data
      if (is.null(aggregateAGPdata$subjectid_chr)) {
        aggregateAGPdata$subjectid_chr <- as.character(aggregateAGPdata$subjectid)
      }
    }

    ## ------------------ 4.3 Apply mapping: grouping and/or filtering -----------
    if (!is.null(mapping_df)) {

      ### ------------------ 4.3.1 Grouping (only if group_var provided) ----------
      if (!is.null(group_var)) {

        if (!group_var %in% names(mapping_df)) {
          warning(
            "Mapping file does not contain group_var column: ", group_var,
            ". Proceeding without grouping."
          )
          group_var <- NULL

        } else {

          idx <- match(aggregateAGPdata$subjectid_chr, mapping_df$subjectid_chr)

          # add group_var column and standardized 'group' column
          aggregateAGPdata[[group_var]] <- mapping_df[[group_var]][idx]
          aggregateAGPdata$group <- as.character(aggregateAGPdata[[group_var]])
          aggregateAGPdata$group[!nzchar(aggregateAGPdata$group)] <- NA_character_

          if (all(is.na(aggregateAGPdata$group))) {
            warning(
              "group_var was specified, but no matching subjectid was found in the mapping file.\n",
              "Proceeding without grouping."
            )
            group_var <- NULL
            aggregateAGPdata$group <- NULL
          }
        }
      }

      ### ------------------ 4.3.2 Filtering (independent of group_var) ----------
      if (!is.null(filter_var)) {

        if (!is.character(filter_var) || length(filter_var) != 1L || !nzchar(filter_var)) {
          stop("`filter_var` must be a non-empty single string like 'male==1 & (age>60 | egfr<90)'.", call. = FALSE)
        }

        rule_all <- trimws(filter_var)

        # ---- atomic rule parser: returns list(var, op, rhs) ----
        .parse_atomic_rule <- function(rule) {

          rule <- trimws(rule)
          if (!nzchar(rule)) stop("Empty atomic rule found.", call. = FALSE)

          # between: var between a and b
          m1 <- regexec("^\\s*([A-Za-z0-9_.]+)\\s+between\\s+(.+?)\\s+and\\s+(.+?)\\s*$",
                        rule, ignore.case = TRUE)
          mm1 <- regmatches(rule, m1)[[1L]]
          if (length(mm1) == 4L) {
            return(list(var = mm1[2L], op = "between", rhs = c(mm1[3L], mm1[4L])))
          }

          # in a:b  (range)
          m2 <- regexec("^\\s*([A-Za-z0-9_.]+)\\s+in\\s+(.+?)\\s*:\\s*(.+?)\\s*$",
                        rule, ignore.case = TRUE)
          mm2 <- regmatches(rule, m2)[[1L]]
          if (length(mm2) == 4L) {
            return(list(var = mm2[2L], op = "between", rhs = c(mm2[3L], mm2[4L])))
          }

          # comparison / equality: var == rhs / != / >= / <= / > / <
          m <- regexec("^\\s*([A-Za-z0-9_.]+)\\s*(==|=|!=|>=|<=|>|<)\\s*(.+?)\\s*$", rule)
          mm <- regmatches(rule, m)[[1L]]
          if (length(mm) == 4L) {
            op <- mm[3L]
            if (op == "=") op <- "=="
            return(list(var = mm[2L], op = op, rhs = mm[4L]))
          }

          stop(
            "Atomic rule not recognized: '", rule, "'. Examples:\n",
            "  male==1\n",
            "  group!=A\n",
            "  site==1,2,3\n",
            "  age>60\n",
            "  egfr<=90\n",
            "  age between 40 and 60\n",
            call. = FALSE
          )
        }

        # ---- evaluate ONE atomic rule on mapping_df -> logical vector ----
        .eval_atomic_rule <- function(rule, mapping_df) {
          pr <- .parse_atomic_rule(rule)

          fv <- tolower(trimws(pr$var))  # mapping_df names already lowercased
          if (!fv %in% names(mapping_df)) {
            stop(
              "filter column '", fv, "' not found in mapping file.\n",
              "Available columns: ", paste(names(mapping_df), collapse = ", "),
              call. = FALSE
            )
          }

          v_raw <- mapping_df[[fv]]
          keep_idx <- rep(FALSE, nrow(mapping_df))

          if (identical(pr$op, "between")) {
            lo <- suppressWarnings(as.numeric(trimws(pr$rhs[1L])))
            hi <- suppressWarnings(as.numeric(trimws(pr$rhs[2L])))
            if (!is.finite(lo) || !is.finite(hi)) {
              stop("`between` bounds must be numeric: ", rule, call. = FALSE)
            }
            if (lo > hi) { tmp <- lo; lo <- hi; hi <- tmp }

            v_num <- suppressWarnings(as.numeric(as.character(v_raw)))
            keep_idx <- is.finite(v_num) & (v_num >= lo) & (v_num <= hi)

          } else if (pr$op %in% c(">", ">=", "<", "<=")) {
            thr <- suppressWarnings(as.numeric(trimws(pr$rhs)))
            if (!is.finite(thr)) stop("Threshold must be numeric for rule: ", rule, call. = FALSE)

            v_num <- suppressWarnings(as.numeric(as.character(v_raw)))
            ok <- is.finite(v_num)

            keep_idx <- switch(
              pr$op,
              ">"  = ok & (v_num >  thr),
              ">=" = ok & (v_num >= thr),
              "<"  = ok & (v_num <  thr),
              "<=" = ok & (v_num <= thr)
            )

          } else if (pr$op %in% c("==", "!=")) {
            rhs <- trimws(pr$rhs)
            rhs <- gsub("^['\"]|['\"]$", "", rhs)  # strip surrounding quotes

            keep_chr <- trimws(unlist(strsplit(rhs, "\\s*,\\s*")))
            keep_chr <- keep_chr[nzchar(keep_chr)]
            if (!length(keep_chr)) stop("No values parsed from rule: ", rule, call. = FALSE)

            v_chr <- as.character(v_raw)
            hit <- v_chr %in% keep_chr

            keep_idx <- if (pr$op == "==") hit else !hit
            keep_idx[is.na(v_chr)] <- FALSE

          } else {
            stop("Unsupported operator in atomic rule: ", pr$op, call. = FALSE)
          }

          as.logical(keep_idx)
        }

        # =============================================================================
        #  Parentheses + AND/OR parser (Shunting-yard -> RPN -> evaluate)
        # =============================================================================

        # Tokenize expression into: ATOM, &, |, (, )
        .tokenize_filter <- function(expr) {
          s <- trimws(expr)
          if (!nzchar(s)) stop("Empty filter expression.", call. = FALSE)

          tokens <- character()
          i <- 1L
          n <- nchar(s)

          while (i <= n) {
            ch <- substr(s, i, i)

            # skip whitespace
            if (grepl("\\s", ch)) {
              i <- i + 1L
              next
            }

            # operators / parentheses
            if (ch %in% c("&", "|", "(", ")")) {
              tokens <- c(tokens, ch)
              i <- i + 1L
              next
            }

            # otherwise read an ATOM until next top-level operator/paren
            j <- i
            while (j <= n) {
              cj <- substr(s, j, j)
              if (cj %in% c("&", "|", "(", ")")) break
              j <- j + 1L
            }

            atom <- trimws(substr(s, i, j - 1L))
            if (!nzchar(atom)) stop("Empty atomic token near position ", i, call. = FALSE)
            tokens <- c(tokens, atom)
            i <- j
          }

          tokens
        }

        # Convert tokens -> RPN using shunting-yard
        .to_rpn <- function(tokens) {
          # precedence: & > |
          prec <- function(op) if (op == "&") 2L else if (op == "|") 1L else 0L

          output <- character()
          stack  <- character()

          for (tk in tokens) {
            if (tk %in% c("&", "|")) {
              while (length(stack) > 0L) {
                top <- stack[length(stack)]
                if (top %in% c("&", "|") && prec(top) >= prec(tk)) {
                  output <- c(output, top)
                  stack  <- stack[-length(stack)]
                } else break
              }
              stack <- c(stack, tk)

            } else if (tk == "(") {
              stack <- c(stack, tk)

            } else if (tk == ")") {
              found <- FALSE
              while (length(stack) > 0L) {
                top <- stack[length(stack)]
                stack <- stack[-length(stack)]
                if (top == "(") { found <- TRUE; break }
                output <- c(output, top)
              }
              if (!found) stop("Mismatched ')' in filter expression.", call. = FALSE)

            } else {
              # ATOM
              output <- c(output, tk)
            }
          }

          # pop remaining
          while (length(stack) > 0L) {
            top <- stack[length(stack)]
            stack <- stack[-length(stack)]
            if (top %in% c("(", ")")) stop("Mismatched '(' in filter expression.", call. = FALSE)
            output <- c(output, top)
          }

          output
        }

        # Evaluate RPN -> logical vector
        .eval_rpn <- function(rpn, mapping_df) {
          st <- list()

          for (tk in rpn) {
            if (tk %in% c("&", "|")) {
              if (length(st) < 2L) stop("Invalid expression: operator '", tk, "' lacks operands.", call. = FALSE)
              b <- st[[length(st)]]; st <- st[-length(st)]
              a <- st[[length(st)]]; st <- st[-length(st)]

              if (!is.logical(a) || !is.logical(b) ||
                  length(a) != nrow(mapping_df) || length(b) != nrow(mapping_df)) {
                stop("Internal error: operand is not a logical vector of mapping_df length.", call. = FALSE)
              }

              st[[length(st) + 1L]] <- if (tk == "&") (a & b) else (a | b)

            } else {
              st[[length(st) + 1L]] <- .eval_atomic_rule(tk, mapping_df)
            }
          }

          if (length(st) != 1L) stop("Invalid filter expression: leftover tokens after evaluation.", call. = FALSE)
          st[[1L]]
        }

        # ---- run ----
        tokens <- .tokenize_filter(rule_all)
        rpn    <- .to_rpn(tokens)
        keep_idx_map <- .eval_rpn(rpn, mapping_df)

        keep_ids <- unique(mapping_df$subjectid_chr[keep_idx_map])

        if (!length(keep_ids)) {
          stop(
            "No subjects matched filter rule in mapping file: ", filter_var, "\n",
            "Tip: check mapping values (e.g., '1' vs 1, blanks/NA, or thresholds).",
            call. = FALSE
          )
        }

        before_n <- length(unique(aggregateAGPdata$subjectid_chr))
        aggregateAGPdata <- aggregateAGPdata[aggregateAGPdata$subjectid_chr %in% keep_ids, , drop = FALSE]
        after_n  <- length(unique(aggregateAGPdata$subjectid_chr))

        if (!nrow(aggregateAGPdata)) {
          stop("No CGM records left after filtering subjects.", call. = FALSE)
        }

        message("Filtered subjects: ", after_n, " / ", before_n, " (rule: ", filter_var, ")")
      }
    } # end if (!is.null(mapping_df))

  } # end if (need_map)

  # refresh local timestamp vector after mapping/filtering
  ts_local <- aggregateAGPdata$timestamp

  # ------------------------------- 5. Quantiles ------------------------------

  ## ---------- 5.1 Quantiles per bin_mins (for Tukey / LOESS) ----------------
  split_list <- split(
    aggregateAGPdata$sensorglucose,
    aggregateAGPdata$hourmin_minutes
  )

  bin_minutes <- suppressWarnings(as.numeric(names(split_list)))
  ord         <- order(bin_minutes)

  split_list  <- split_list[ord]
  bin_minutes <- bin_minutes[ord]

  qmat <- t(vapply(
    split_list,
    FUN = function(x) {
      stats::quantile(
        x,
        probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
        na.rm = TRUE,
        names = FALSE
      )
    },
    FUN.VALUE = numeric(5L)
  ))

  quartiles <- data.frame(
    hourmin_minutes      = bin_minutes,
    sensorglucose5perc   = qmat[, 1],
    sensorglucoseqone    = qmat[, 2],
    sensorglucosemedian  = qmat[, 3],
    sensorglucoseqthree  = qmat[, 4],
    sensorglucose95perc  = qmat[, 5]
  )

  ## ---------- 5.2 Quantiles per native bin (for all "original" AGPs) --------
  split_list_o <- split(
    aggregateAGPdata$sensorglucose,
    aggregateAGPdata$hourmin_minutes_orig
  )

  bin_minutes_o <- suppressWarnings(as.numeric(names(split_list_o)))
  ord_o         <- order(bin_minutes_o)

  split_list_o  <- split_list_o[ord_o]
  bin_minutes_o <- bin_minutes_o[ord_o]

  qmat_o <- t(vapply(
    split_list_o,
    FUN = function(x) {
      stats::quantile(
        x,
        probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
        na.rm = TRUE,
        names = FALSE
      )
    },
    FUN.VALUE = numeric(5L)
  ))

  quartiles_original <- data.frame(
    hourmin_minutes      = bin_minutes_o,
    sensorglucose5perc   = qmat_o[, 1],
    sensorglucoseqone    = qmat_o[, 2],
    sensorglucosemedian  = qmat_o[, 3],
    sensorglucoseqthree  = qmat_o[, 4],
    sensorglucose95perc  = qmat_o[, 5]
  )

  # For "original" AGP also wrap to 24:00 (1440 minutes)
  if (!any(quartiles_original$hourmin_minutes == 1440)) {
    first_row_o <- quartiles_original[which.min(quartiles_original$hourmin_minutes), ]
    wrap_row_o  <- first_row_o
    wrap_row_o$hourmin_minutes <- 1440

    quartiles_original <- rbind(quartiles_original, wrap_row_o)
    quartiles_original <- quartiles_original[order(quartiles_original$hourmin_minutes), ]
  }

  # ------------------ 6. Tukey 3R smoothing ---------------------------------
  smooth_safe <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 3L) return(x)

    # ---- handle NA: fill by linear interpolation, then carry ends ----
    if (anyNA(x)) {
      idx <- which(!is.na(x))
      if (length(idx) < 2L) return(x)  # not enough info to fill

      # linear interpolation on internal NAs
      x_filled <- x
      x_filled <- stats::approx(
        x = idx,
        y = x[idx],
        xout = seq_along(x),
        method = "linear",
        rule = 2  # carry ends
      )$y
      x <- x_filled
    }

    # Prefer stats::smooth(kind="3R") if available; otherwise fall back safely
    if (exists("smooth", where = "package:stats", inherits = TRUE) ||
        exists("smooth", envir = as.environment("package:stats"), inherits = TRUE)) {
      return(as.numeric(stats::smooth(x, kind = "3R", twiceit = TRUE)))
    }

    # Fallback: two-pass running median (stable + no extra deps)
    x1 <- stats::runmed(x, k = 3, endrule = "median")
    x2 <- stats::runmed(x1, k = 3, endrule = "median")
    x2
  }

  quartiles$smooth5perc  <- smooth_safe(quartiles$sensorglucose5perc)
  quartiles$smoothqone   <- smooth_safe(quartiles$sensorglucoseqone)
  quartiles$smoothmed    <- smooth_safe(quartiles$sensorglucosemedian)
  quartiles$smoothqthree <- smooth_safe(quartiles$sensorglucoseqthree)
  quartiles$smooth95perc <- smooth_safe(quartiles$sensorglucose95perc)

  # Wrap to 24:00 (duplicate the first row at 1440 minutes)
  if (!any(quartiles$hourmin_minutes == 1440)) {
    # Take the earliest row (usually 0 minutes)
    first_row <- quartiles[which.min(quartiles$hourmin_minutes), ]
    wrap_row  <- first_row
    wrap_row$hourmin_minutes <- 1440

    quartiles <- rbind(quartiles, wrap_row)
    quartiles <- quartiles[order(quartiles$hourmin_minutes), ]
  }

  ## ------------------ 6.1 Group-specific quartiles (optional) ---------------
  quartiles_grouped          <- NULL  # for Tukey / LOESS (by bin_mins)
  quartiles_grouped_original <- NULL  # for original (by native sampling interval)

  if (!is.null(group_var) && "group" %in% names(aggregateAGPdata)) {
    aggregateAGPdata <- aggregateAGPdata[!is.na(aggregateAGPdata$group), , drop = FALSE]
    split_by_group <- split(aggregateAGPdata, aggregateAGPdata$group)

    ## 6.1a) Quantiles by bin_mins for each group (Tukey / LOESS)
    quartiles_grouped_list <- lapply(names(split_by_group), function(g) {
      dat_g <- split_by_group[[g]]

      split_list_g <- split(dat_g$sensorglucose, dat_g$hourmin_minutes)
      if (!length(split_list_g)) return(NULL)

      bin_minutes_g <- suppressWarnings(as.numeric(names(split_list_g)))
      ord_g         <- order(bin_minutes_g)

      split_list_g  <- split_list_g[ord_g]
      bin_minutes_g <- bin_minutes_g[ord_g]

      qmat_g <- t(vapply(
        split_list_g,
        FUN = function(x) {
          stats::quantile(
            x,
            probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
            na.rm = TRUE,
            names = FALSE
          )
        },
        FUN.VALUE = numeric(5L)
      ))

      quart_g <- data.frame(
        hourmin_minutes      = bin_minutes_g,
        sensorglucose5perc   = qmat_g[, 1],
        sensorglucoseqone    = qmat_g[, 2],
        sensorglucosemedian  = qmat_g[, 3],
        sensorglucoseqthree  = qmat_g[, 4],
        sensorglucose95perc  = qmat_g[, 5]
      )

      quart_g$smooth5perc  <- smooth_safe(quart_g$sensorglucose5perc)
      quart_g$smoothqone   <- smooth_safe(quart_g$sensorglucoseqone)
      quart_g$smoothmed    <- smooth_safe(quart_g$sensorglucosemedian)
      quart_g$smoothqthree <- smooth_safe(quart_g$sensorglucoseqthree)
      quart_g$smooth95perc <- smooth_safe(quart_g$sensorglucose95perc)

      if (!any(quart_g$hourmin_minutes == 1440)) {
        first_row_g <- quart_g[which.min(quart_g$hourmin_minutes), ]
        wrap_row_g  <- first_row_g
        wrap_row_g$hourmin_minutes <- 1440
        quart_g <- rbind(quart_g, wrap_row_g)
        quart_g <- quart_g[order(quart_g$hourmin_minutes), ]
      }

      quart_g$group <- g
      quart_g
    })

    if (length(quartiles_grouped_list)) {
      quartiles_grouped <- do.call(
        rbind,
        quartiles_grouped_list[
          !vapply(quartiles_grouped_list, is.null, logical(1L))
        ]
      )
    }

    ## ------------ 6.2 Quantiles by native sampling interval for each group --------------
    quartiles_grouped_orig_list <- lapply(names(split_by_group), function(g) {
      dat_g <- split_by_group[[g]]

      split_list_go <- split(dat_g$sensorglucose, dat_g$hourmin_minutes_orig)
      if (!length(split_list_go)) return(NULL)

      bin_minutes_go <- suppressWarnings(as.numeric(names(split_list_go)))
      ord_go         <- order(bin_minutes_go)

      split_list_go  <- split_list_go[ord_go]
      bin_minutes_go <- bin_minutes_go[ord_go]

      qmat_go <- t(vapply(
        split_list_go,
        FUN = function(x) {
          stats::quantile(
            x,
            probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
            na.rm = TRUE,
            names = FALSE
          )
        },
        FUN.VALUE = numeric(5L)
      ))

      quart_go <- data.frame(
        hourmin_minutes      = bin_minutes_go,
        sensorglucose5perc   = qmat_go[, 1],
        sensorglucoseqone    = qmat_go[, 2],
        sensorglucosemedian  = qmat_go[, 3],
        sensorglucoseqthree  = qmat_go[, 4],
        sensorglucose95perc  = qmat_go[, 5]
      )

      if (!any(quart_go$hourmin_minutes == 1440)) {
        first_row_go <- quart_go[which.min(quart_go$hourmin_minutes), ]
        wrap_row_go  <- first_row_go
        wrap_row_go$hourmin_minutes <- 1440
        quart_go <- rbind(quart_go, wrap_row_go)
        quart_go <- quart_go[order(quart_go$hourmin_minutes), ]
      }

      quart_go$group <- g
      quart_go
    })

    if (length(quartiles_grouped_orig_list)) {
      quartiles_grouped_original <- do.call(
        rbind,
        quartiles_grouped_orig_list[
          !vapply(quartiles_grouped_orig_list, is.null, logical(1L))
        ]
      )
    }
  }

  # ------------------------------ 7. AGP Plots ------------------------------

  ## ---- 7.0 Local scale helpers (no global side effects) --------------------
  # If viridis is installed, use it for discrete color/fill scales.
  # If not installed, fall back to ggplot2 defaults silently.
  .maybe_scale_color_viridis_d <- function(limits = NULL, option = "D", end = 1, ...) {

    # ---- viridis (newer): has scale_color_viridis_d ----
    if (requireNamespace("viridis", quietly = TRUE)) {
      exports <- tryCatch(getNamespaceExports("viridis"), error = function(e) character())
      if ("scale_color_viridis_d" %in% exports) {
        f <- getExportedValue("viridis", "scale_color_viridis_d")
        return(f(limits = limits, option = option, end = end, ...))
      }

      # ---- viridis (older): use scale_color_viridis(discrete=TRUE) ----
      if ("scale_color_viridis" %in% exports) {
        f <- getExportedValue("viridis", "scale_color_viridis")
        return(f(discrete = TRUE, limits = limits, option = option, end = end, ...))
      }
    }

    # ---- viridisLite fallback ----
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      if (!is.null(limits)) {
        cols <- viridisLite::viridis(length(limits), option = option, end = end, ...)
        names(cols) <- limits
        return(ggplot2::scale_color_manual(values = cols, limits = limits))
      }
      return(
        ggplot2::discrete_scale(
          aesthetics = "colour",
          scale_name = "viridisLite",
          palette = function(n) viridisLite::viridis(n, option = option, end = end, ...)
        )
      )
    }

    NULL
  }

  .maybe_scale_fill_viridis_d <- function(limits = NULL, option = "D", end = 1, ...) {

    # ---- viridis (newer): has scale_fill_viridis_d ----
    if (requireNamespace("viridis", quietly = TRUE)) {
      exports <- tryCatch(getNamespaceExports("viridis"), error = function(e) character())
      if ("scale_fill_viridis_d" %in% exports) {
        f <- getExportedValue("viridis", "scale_fill_viridis_d")
        return(f(limits = limits, option = option, end = end, ...))
      }

      # ---- viridis (older): use scale_fill_viridis(discrete=TRUE) ----
      if ("scale_fill_viridis" %in% exports) {
        f <- getExportedValue("viridis", "scale_fill_viridis")
        return(f(discrete = TRUE, limits = limits, option = option, end = end, ...))
      }
    }

    # ---- viridisLite fallback ----
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      if (!is.null(limits)) {
        cols <- viridisLite::viridis(length(limits), option = option, end = end, ...)
        names(cols) <- limits
        return(ggplot2::scale_fill_manual(values = cols, limits = limits))
      }
      return(
        ggplot2::discrete_scale(
          aesthetics = "fill",
          scale_name = "viridisLite",
          palette = function(n) viridisLite::viridis(n, option = option, end = end, ...)
        )
      )
    }

    NULL
  }

  ## ------------------ 7.1 AGP figure formatting helpers --------------------

  build_agg_agp <- function(df,
                            lower95, lower50,
                            med, upper50, upper95,
                            title_main, median_label,
                            median_size = 1.2,
                            yaxis, ref_lines) {

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[["hourmin_minutes"]])
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[[lower95]],
          ymax = .data[[upper95]],
          fill = "IQR90"
        ),
        alpha = 0.15
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[[lower50]],
          ymax = .data[[upper50]],
          fill = "IQR50"
        ),
        alpha = 0.30
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[[med]], color = median_label),
        linewidth = median_size
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1440),
        breaks = seq(0, 1440, by = 180),
        labels = c(
          "0:00", "03:00", "06:00", "09:00",
          "12:00", "15:00", "18:00", "21:00", "24:00"
        ),
        expand = c(0.005, 0.005)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(yaxis[1], yaxis[2], by = 50),
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = yaxis) +
      ggplot2::scale_color_manual(
        name   = "",
        values = stats::setNames("#3D7359", median_label)
      ) +
      ggplot2::scale_fill_manual(
        name   = "",
        values = c(
          "IQR90" = "#A0D8B3",
          "IQR50" = "#76B798"
        )
      ) +
      ggplot2::labs(
        title = title_main,
        x     = "Time of day (24-hour clock)",
        y     = "Glucose (mg/dL)"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.7),
        axis.title      = ggplot2::element_text(size = 14, face = "bold"),
        axis.text       = ggplot2::element_text(size = 11),
        legend.position = "right",
        legend.title    = ggplot2::element_blank(),
        legend.text     = ggplot2::element_text(size = 11),
        plot.title      = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      )

    if (length(ref_lines)) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = ref_lines,
          linetype   = "dashed",
          color      = "#084F6A",
          linewidth  = 0.7
        )
    }
    p
  }

  # Helper for group-wise aggregate AGP
  build_agg_agp_grouped <- function(df,
                                    lower95, lower50,
                                    med, upper50, upper95,
                                    title_main,
                                    median_size = 1.2,
                                    yaxis, ref_lines,
                                    group_levels = NULL) {

    if (is.null(df) || !nrow(df)) return(NULL)

    df$group <- as.character(df$group)

    # Use a global, fixed group order if provided
    if (is.null(group_levels)) {
      lv <- sort(unique(df$group))
    } else {
      lv <- as.character(group_levels)
    }

    # Force consistent factor levels across plots
    df$group <- factor(df$group, levels = lv)

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[["hourmin_minutes"]])
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[[lower95]],
          ymax = .data[[upper95]],
          fill = .data[["group"]]
        ),
        alpha = 0.15
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[[lower50]],
          ymax = .data[[upper50]],
          fill = .data[["group"]]
        ),
        alpha = 0.30
      ) +
      ggplot2::geom_line(
        ggplot2::aes(
          y = .data[[med]],
          color = .data[["group"]]
        ),
        linewidth = median_size
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1440),
        breaks = seq(0, 1440, by = 180),
        labels = c(
          "0:00", "03:00", "06:00", "09:00",
          "12:00", "15:00", "18:00", "21:00", "24:00"
        ),
        expand = c(0.005, 0.005)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(yaxis[1], yaxis[2], by = 50),
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(ylim = yaxis) +
      ggplot2::labs(
        title = title_main,
        x     = "Time of day (24-hour clock)",
        y     = "Glucose (mg/dL)"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.7),
        axis.title       = ggplot2::element_text(size = 14, face = "bold"),
        axis.text        = ggplot2::element_text(size = 11),
        legend.position  = "right",
        legend.title     = ggplot2::element_blank(),
        legend.text      = ggplot2::element_text(size = 11),
        plot.title       = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::guides(fill = "none")

    sc_col <- .maybe_scale_color_viridis_d(limits = lv)
    if (!is.null(sc_col)) p <- p + sc_col

    sc_fill <- .maybe_scale_fill_viridis_d(limits = lv)
    if (!is.null(sc_fill)) p <- p + sc_fill

    if (length(ref_lines)) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = ref_lines,
          linetype   = "dashed",
          color      = "#084F6A",
          linewidth  = 0.7
        )
    }

    p
  }

  ## ------------------ 7.2 Aggregate AGP (overall) ---------------------------

  # 1) Original
  aggAGP_original <- build_agg_agp(
    df         = quartiles,
    lower95    = "sensorglucose5perc",
    lower50    = "sensorglucoseqone",
    med        = "sensorglucosemedian",
    upper50    = "sensorglucoseqthree",
    upper95    = "sensorglucose95perc",
    title_main = "Aggregate Ambulatory Glucose Profile (Original, no smoothing)",
    median_label = "Median",
    median_size  = 1.0,
    yaxis        = yaxis,
    ref_lines    = ref_lines
  )

  # 2) Tukey
  aggAGPtukey <- build_agg_agp(
    df         = quartiles,
    lower95    = "smooth5perc",
    lower50    = "smoothqone",
    med        = "smoothmed",
    upper50    = "smoothqthree",
    upper95    = "smooth95perc",
    title_main = "Aggregate Ambulatory Glucose Profile (Tukey 3R smoothing)",
    median_label = "Median",
    median_size  = 1.3,
    yaxis        = yaxis,
    ref_lines    = ref_lines
  )

  # 3) Loess
  loess_med <- try(
    stats::loess(
      sensorglucosemedian ~ hourmin_minutes,
      data = quartiles,
      span = 0.3
    ),
    silent = TRUE
  )

  if (!inherits(loess_med, "try-error")) {
    quartiles$loess_med <- stats::predict(loess_med, quartiles$hourmin_minutes)
  } else {
    quartiles$loess_med <- quartiles$sensorglucosemedian
  }

  aggAGPloess <- build_agg_agp(
    df         = quartiles,
    lower95    = "sensorglucose5perc",
    lower50    = "sensorglucoseqone",
    med        = "loess_med",
    upper50    = "sensorglucoseqthree",
    upper95    = "sensorglucose95perc",
    title_main = "Aggregate Ambulatory Glucose Profile (LOESS smoothing)",
    median_label = "Median",
    median_size  = 1.2,
    yaxis        = yaxis,
    ref_lines    = ref_lines
  )

  ## --------------- 7.3 Grouped aggregate AGP --------------------------------
  aggAGP_original_grouped <- NULL
  aggAGPtukey_grouped     <- NULL
  aggAGPloess_grouped     <- NULL

  # Global group order for consistent colors across grouped plots
  group_levels <- NULL
  if (!is.null(group_var) && "group" %in% names(aggregateAGPdata)) {
    group_levels <- sort(unique(stats::na.omit(as.character(aggregateAGPdata$group))))
  }

  ## title
  smooth_label <- switch(
    agg_smooth,
    "loess"    = "LOESS smoothing",
    "tukey"    = "Tukey 3R smoothing",
    "original" = "no smoothing",
    # default
    agg_smooth
  )

  AGPloess_title <- paste0(
    "Subject-level 24-hour Glucose Profiles (", smooth_label, ")"
  )

  method_label <- switch(
    agg_smooth,
    "original" = "Original, no smoothing",
    "tukey"    = "Tukey 3R smoothing",
    "loess"    = "LOESS smoothing"
  )

  title_main <- paste0(
    "Aggregate Ambulatory Glucose Profile by Group (", method_label, ")"
  )

  # Decouple original vs tukey/loess availability
  if (!is.null(group_var) &&
      !is.null(group_levels) &&
      length(group_levels) > 0L) {

    # 1) Original grouped (native interval)
    if (!is.null(quartiles_grouped_original) &&
        nrow(quartiles_grouped_original) > 0L) {

      # aggAGP_original_grouped <- build_agg_agp_grouped(
      #   df         = quartiles_grouped_original,
      #   lower95    = "sensorglucose5perc",
      #   lower50    = "sensorglucoseqone",
      #   med        = "sensorglucosemedian",
      #   upper50    = "sensorglucoseqthree",
      #   upper95    = "sensorglucose95perc",
      #   title_main = title_main,
      #   median_size  = 1.0,
      #   yaxis        = yaxis,
      #   ref_lines    = ref_lines,
      #   group_levels = group_levels
      # )

      aggAGP_original_grouped <- build_agg_agp_grouped(
        df         = quartiles_grouped,
        lower95    = "sensorglucose5perc",
        lower50    = "sensorglucoseqone",
        med        = "sensorglucosemedian",
        upper50    = "sensorglucoseqthree",
        upper95    = "sensorglucose95perc",
        title_main = title_main,
        median_size  = 1.0,
        yaxis        = yaxis,
        ref_lines    = ref_lines,
        group_levels = group_levels
      )
    }

    # 2) Tukey + LOESS grouped (bin_mins)
    if (!is.null(quartiles_grouped) &&
        nrow(quartiles_grouped) > 0L) {

      aggAGPtukey_grouped <- build_agg_agp_grouped(
        df         = quartiles_grouped,
        lower95    = "smooth5perc",
        lower50    = "smoothqone",
        med        = "smoothmed",
        upper50    = "smoothqthree",
        upper95    = "smooth95perc",
        title_main = title_main,
        median_size  = 1.3,
        yaxis        = yaxis,
        ref_lines    = ref_lines,
        group_levels = group_levels
      )

      # LOESS for each group separately
      quartiles_grouped$loess_med <- NA_real_
      for (g in stats::na.omit(unique(quartiles_grouped$group))) {
        idx_g <- quartiles_grouped$group == g
        df_g  <- quartiles_grouped[idx_g, , drop = FALSE]

        loess_med_g <- try(
          stats::loess(
            sensorglucosemedian ~ hourmin_minutes,
            data = df_g,
            span = 0.3
          ),
          silent = TRUE
        )

        if (!inherits(loess_med_g, "try-error")) {
          quartiles_grouped$loess_med[idx_g] <-
            stats::predict(loess_med_g, df_g$hourmin_minutes)
        } else {
          quartiles_grouped$loess_med[idx_g] <- df_g$sensorglucosemedian
        }
      }

      aggAGPloess_grouped <- build_agg_agp_grouped(
        df         = quartiles_grouped,
        lower95    = "sensorglucose5perc",
        lower50    = "sensorglucoseqone",
        med        = "loess_med",
        upper50    = "sensorglucoseqthree",
        upper95    = "sensorglucose95perc",
        title_main = title_main,
        median_size  = 1.2,
        yaxis        = yaxis,
        ref_lines    = ref_lines,
        group_levels = group_levels
      )
    }
  }

  ## ---------------- 7.4 Subject-level 24h overlay (optional) -----------------
  # We draw one curve per subject. To make "original/tukey/loess" comparable,
  # we first bin into bin_mins and compute per-subject median in each bin,
  # then apply the chosen smoothing to the median series.

  AGPloess <- NULL

  if (isTRUE(subject_overlay)) {

  # refresh local timestamp vector after mapping/filtering
  ts_local <- aggregateAGPdata$timestamp

  ### ----------- 7.4.1 Build per-subject median-by-bin (bin_mins) --------------
  # Use floor_date for stable left-edge binning (recommended for AGP)
  ts_bin_overlay <- lubridate::floor_date(ts_local, unit = paste(bin_mins, "mins"))
  hourmin_overlay <- lubridate::hour(ts_bin_overlay) * 60 + lubridate::minute(ts_bin_overlay)

  overlay_df0 <- aggregateAGPdata
  overlay_df0$hourmin_minutes_bin <- hourmin_overlay

  overlay_med <- overlay_df0 |>
    dplyr::filter(stats::complete.cases(.data$subjectid, .data$hourmin_minutes_bin, .data$sensorglucose)) |>
    dplyr::group_by(.data$subjectid, .data$hourmin_minutes_bin) |>
    dplyr::summarise(med = stats::median(.data$sensorglucose, na.rm = TRUE), .groups = "drop")

  # complete a full 24h grid for each subject so curves align (0..1440)
  grid_bins <- seq(0, 1440, by = bin_mins)
  overlay_med <- overlay_med |>
    dplyr::group_by(.data$subjectid) |>
    tidyr::complete(hourmin_minutes_bin = grid_bins) |>
    dplyr::arrange(.data$hourmin_minutes_bin, .by_group = TRUE) |>
    dplyr::ungroup()

  ### ------- 7.4.2 Apply smoothing to each subject's median series depending on agg_smooth -------
  # smooth_safe() is already defined above (Tukey 3R fallback-safe)
  overlay_med <- overlay_med |>
    dplyr::group_by(.data$subjectid) |>
    dplyr::arrange(.data$hourmin_minutes_bin, .by_group = TRUE) |>
    dplyr::mutate(
      curve = dplyr::case_when(
        agg_smooth == "original" ~ .data$med,

        agg_smooth == "tukey" ~ {
          # smooth_safe() MUST be NA-safe (fill NA before stats::smooth)
          smooth_safe(.data$med)
        },

        agg_smooth == "loess" ~ {
          d <- dplyr::cur_data()

          # fit LOESS only on non-missing med
          d_fit <- d[!is.na(d$med) & !is.na(d$hourmin_minutes_bin), , drop = FALSE]

          if (nrow(d_fit) < 3L) {
            d$med
          } else {
            fit <- try(
              stats::loess(med ~ hourmin_minutes_bin, data = d_fit, span = 0.3),
              silent = TRUE
            )

            if (!inherits(fit, "try-error")) {
              as.numeric(stats::predict(
                fit,
                newdata = data.frame(hourmin_minutes_bin = d$hourmin_minutes_bin)
              ))
            } else {
              d$med
            }
          }
        },

        TRUE ~ .data$med
      )
    ) |>
    dplyr::ungroup()

  ### ------------------- Plot (points + smooth curves) ------------------------
  points_df <- aggregateAGPdata |>
    dplyr::filter(stats::complete.cases(.data$subjectid, .data$time_minutes, .data$sensorglucose)) |>
    dplyr::transmute(
      subjectid = .data$subjectid,
      x = .data$time_minutes,          # continuous minutes since midnight
      y = .data$sensorglucose
    )

  AGPloess <- ggplot2::ggplot() +
    # 1) raw points (scatter)
    ggplot2::geom_point(
      data = points_df,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$subjectid),
      alpha = 0.2,
      size  = 0.4,
      stroke = 0,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    # 2) smoothed curves (overlay_med result)
    ggplot2::geom_line(
      data = overlay_med,
      ggplot2::aes(
        x     = .data$hourmin_minutes_bin,
        y     = .data$curve,
        color = .data$subjectid
      ),
      linewidth = 0.9,
      alpha = 0.95,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    .maybe_scale_color_viridis_d() +
    ggplot2::scale_x_continuous(
      limits = c(0, 1440),
      breaks = seq(0, 1440, by = 180),
      labels = c("0:00","03:00","06:00","09:00","12:00","15:00","18:00","21:00","24:00"),
      expand = c(0.005, 0.005)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(yaxis[1], yaxis[2], by = 50),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = yaxis) +
    ggplot2::labs(
      title = AGPloess_title,
      x     = "Time of day (24-hour clock)",
      y     = "Glucose (mg/dL)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title       = ggplot2::element_text(size = 12, face = "bold"),
      axis.text        = ggplot2::element_text(size = 9),
      panel.spacing    = grid::unit(1.5, "lines"),
      legend.position  = "right",
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 9),
      plot.title       = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(ref_lines)) {
    AGPloess <- AGPloess +
      ggplot2::geom_hline(
        yintercept = ref_lines,
        linetype   = "dashed",
        color      = "#084F6A",
        linewidth  = 0.7
      )
    }
  }

  ## ------------------ 7.5 Daily glucose profile per subject ---------------

  ### ----------- 7.5.1 Multi−day Ambulatory Glucose Profile -------------

  # Version that outputs a selected AGP type (Original, Tukey, or LOESS)
  plot_daily_subject_cards <- function(aggregateAGPdata, tz, yaxis,
                                       low_thr, high_thr, ref_lines,
                                       outputdir, page_width, page_height,
                                       agg_smooth, bin_mins) {

    # Helper: build multi-day overlay (per subject)
    build_multiday_overlay <- function(df_sid, agg_smooth, yaxis, sid) {

      df_overlay <- df_sid[stats::complete.cases(
        df_sid[, c("time_minutes", "sensorglucose", "date_label")]
      ), ]
      if (!nrow(df_overlay)) return(NULL)

      # Sort by date and time within day
      df_overlay <- df_overlay[order(df_overlay$date_label,
                                     df_overlay$time_minutes), ]

      if (agg_smooth == "original") {
        # No smoothing: plain per-day connected lines
        p_multi <- ggplot2::ggplot(
          df_overlay,
          ggplot2::aes(
            x     = .data[["time_minutes"]],
            y     = .data[["sensorglucose"]],
            color = .data[["date_label"]]
          )
        ) +
          ggplot2::geom_line(linewidth = 0.8, alpha = 0.9)

      } else if (agg_smooth == "tukey") {
        # Tukey 3R: smooth each day separately before drawing lines
        df_sm_list <- lapply(
          split(df_overlay, df_overlay$date_label),
          function(d_day) {
            d_day <- d_day[order(d_day$time_minutes), ]
            if (nrow(d_day) < 3L) {
              # Too few points: keep original series
              d_day$sensorglucose_smooth <- d_day$sensorglucose
            } else {
              d_day$sensorglucose_smooth <- smooth_safe(d_day$sensorglucose)
            }
            d_day
          }
        )

        df_sm <- do.call(
          rbind,
          df_sm_list[!vapply(df_sm_list, is.null, logical(1L))]
        )

        p_multi <- ggplot2::ggplot(
          df_sm,
          ggplot2::aes(
            x     = .data[["time_minutes"]],
            y     = .data[["sensorglucose_smooth"]],
            color = .data[["date_label"]]
          )
        ) +
          ggplot2::geom_line(linewidth = 0.8, alpha = 0.95)

      } else { # agg_smooth == "loess"
        # LOESS: smooth each day using geom_smooth grouped by date_label
        p_multi <- ggplot2::ggplot(
          df_overlay,
          ggplot2::aes(
            x     = .data[["time_minutes"]],
            y     = .data[["sensorglucose"]],
            color = .data[["date_label"]]
          )
        ) +
          ggplot2::geom_smooth(se = FALSE, linewidth = 0.8)
      }

      method_label <- switch(
        agg_smooth,
        "original" = "Original, no smoothing",
        "tukey"    = "Tukey 3R smoothing",
        "loess"    = "LOESS smoothing"
      )

      title_main <- paste0(
        "Multi-day Ambulatory Glucose Profile (", method_label, ")"
      )

      p_multi +
        .maybe_scale_color_viridis_d() +
        ggplot2::scale_x_continuous(
          limits = c(0, 1440),
          breaks = seq(0, 1440, by = 180),
          labels = c(
            "0:00", "03:00", "06:00", "09:00",
            "12:00", "15:00", "18:00", "21:00", "24:00"
          ),
          expand = c(0.005, 0.005)
        ) +
        ggplot2::scale_y_continuous(
          breaks = seq(yaxis[1], yaxis[2], by = 50),
          expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(ylim = yaxis) +
        ggplot2::labs(
          title = title_main,
          x     = "Time of day (24-hour clock)",
          y     = "Glucose (mg/dL)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4),
          axis.title      = ggplot2::element_text(size = 12, face = "bold"),
          axis.text       = ggplot2::element_text(size = 9),
          legend.position = "right",
          legend.title    = ggplot2::element_blank(),
          legend.text     = ggplot2::element_text(size = 8),
          plot.title      = ggplot2::element_text(
            size = 16, face = "bold", hjust = 0.5
          ),
          panel.grid.minor = ggplot2::element_blank()
        )
    }

    ### ----------- 7.5.2 Aggregate Ambulatory Glucose Profile -------------

    # Main loop: output one PDF per subject
    unique_ids <- unique(aggregateAGPdata$subjectid)

    for (sid in unique_ids) {
      df_sid <- aggregateAGPdata[aggregateAGPdata$subjectid == sid, ]
      if (!nrow(df_sid)) next

      # A. Subject-level AGP quantiles (original / tukey / loess)
      df_agp <- df_sid[stats::complete.cases(
        df_sid[, c("sensorglucose", "hourmin_minutes")]
      ), ]
      subjAGP_original <- subjAGPtukey <- subjAGPloess <- NULL

      if (nrow(df_agp) >= 5L) {

        # 1) Subject-specific native sampling interval -> bin width
        dt_sid <- diff(sort(unique(df_sid$timestamp)))
        dt_sid_mins <- as.numeric(dt_sid, units = "mins")
        native_sid <- stats::median(
          dt_sid_mins[is.finite(dt_sid_mins) & dt_sid_mins > 0]
        )

        if (is.na(native_sid) || !is.finite(native_sid) || native_sid <= 0) {
          native_sid <- bin_mins  # fallback to global bin width
        }
        bin_sid <- max(1L, round(native_sid))

        ts_bin_sid <- lubridate::round_date(
          df_agp$timestamp,
          unit = paste(bin_sid, "mins")
        )
        hourmin_sid <- lubridate::hour(ts_bin_sid) * 60 +
          lubridate::minute(ts_bin_sid)

        # 2) Subject-level quantiles at 0.05/0.25/0.50/0.75/0.95
        split_list <- split(df_agp$sensorglucose, hourmin_sid)
        bin_minutes <- suppressWarnings(as.numeric(names(split_list)))
        ord <- order(bin_minutes)
        split_list  <- split_list[ord]
        bin_minutes <- bin_minutes[ord]

        qmat <- t(vapply(
          split_list,
          FUN = function(x) {
            stats::quantile(
              x,
              probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
              na.rm = TRUE,
              names = FALSE
            )
          },
          FUN.VALUE = numeric(5L)
        ))

        quart_sid <- data.frame(
          hourmin_minutes      = bin_minutes,
          sensorglucose5perc   = qmat[, 1],
          sensorglucoseqone    = qmat[, 2],
          sensorglucosemedian  = qmat[, 3],
          sensorglucoseqthree  = qmat[, 4],
          sensorglucose95perc  = qmat[, 5]
        )

        # 3) Tukey 3R smoothing at subject level
        quart_sid$smooth5perc  <- smooth_safe(quart_sid$sensorglucose5perc)
        quart_sid$smoothqone   <- smooth_safe(quart_sid$sensorglucoseqone)
        quart_sid$smoothmed    <- smooth_safe(quart_sid$sensorglucosemedian)
        quart_sid$smoothqthree <- smooth_safe(quart_sid$sensorglucoseqthree)
        quart_sid$smooth95perc <- smooth_safe(quart_sid$sensorglucose95perc)

        # 4) Wrap subject-level AGP to 24:00 (1440 minutes)
        if (!any(quart_sid$hourmin_minutes == 1440)) {
          first_row <- quart_sid[which.min(quart_sid$hourmin_minutes), ]
          wrap_row  <- first_row
          wrap_row$hourmin_minutes <- 1440
          quart_sid <- rbind(quart_sid, wrap_row)
          quart_sid <- quart_sid[order(quart_sid$hourmin_minutes), ]
        }

        # 5) LOESS median curve at subject level
        loess_med_sid <- try(
          stats::loess(
            sensorglucosemedian ~ hourmin_minutes,
            data = quart_sid,
            span = 0.3
          ),
          silent = TRUE
        )
        if (!inherits(loess_med_sid, "try-error")) {
          quart_sid$loess_med <- stats::predict(
            loess_med_sid,
            quart_sid$hourmin_minutes
          )
        } else {
          quart_sid$loess_med <- quart_sid$sensorglucosemedian
        }

        # # 6) Build subject-level AGPs using the same helper as for aggregate plots
        subjAGP_original <- build_agg_agp(
          df         = quart_sid,
          lower95    = "sensorglucose5perc",
          lower50    = "sensorglucoseqone",
          med        = "sensorglucosemedian",
          upper50    = "sensorglucoseqthree",
          upper95    = "sensorglucose95perc",
          title_main = "Aggregate Ambulatory Glucose Profile (Original, no smoothing)",
          median_label = "Median",
          median_size  = 1.0,
          yaxis        = yaxis,
          ref_lines    = ref_lines
        )

        subjAGPtukey <- build_agg_agp(
          df         = quart_sid,
          lower95    = "smooth5perc",
          lower50    = "smoothqone",
          med        = "smoothmed",
          upper50    = "smoothqthree",
          upper95    = "smooth95perc",
          title_main = "Aggregate Ambulatory Glucose Profile (Tukey 3R smoothing)",
          median_label = "Median",
          median_size  = 1.3,
          yaxis        = yaxis,
          ref_lines    = ref_lines
        )

        subjAGPloess <- build_agg_agp(
          df         = quart_sid,
          lower95    = "sensorglucose5perc",
          lower50    = "sensorglucoseqone",
          med        = "loess_med",
          upper50    = "sensorglucoseqthree",
          upper95    = "sensorglucose95perc",
          title_main = "Aggregate Ambulatory Glucose Profile (LOESS smoothing)",
          median_label = "Median",
          median_size  = 1.2,
          yaxis        = yaxis,
          ref_lines    = ref_lines
        )
      }

      # Select one subject-level AGP according to agg_smooth
      subjAGP_selected <- switch(
        agg_smooth,
        "original" = subjAGP_original,
        "tukey"    = subjAGPtukey,
        "loess"    = subjAGPloess
      )

      # B. Daily profile cards (faceted by date with zone coloring)

      # Ensure date and date_label exist
      df_sid$date       <- as.Date(df_sid$timestamp, tz = tz)
      df_sid$date_label <- format(df_sid$date, "%Y-%m-%d")

      df_sid <- df_sid[stats::complete.cases(
        df_sid[, c("time_minutes", "sensorglucose", "date_label")]
      ), ]
      if (!nrow(df_sid)) next

      # Multi-day overlay using the same agg_smooth
      p_multiday <- build_multiday_overlay(
        df_sid     = df_sid,
        agg_smooth = agg_smooth,
        yaxis      = yaxis,
        sid        = sid
      )

      ### ----------- 7.5.3 Daily Ambulatory Glucose Profile -------------

      # Helper: for one day's data, cut the trajectory into segments
      # according to thresholds and assign a zone to each segment
      build_segments_one_day <- function(d_day, low_thr, high_thr) {
        d_day <- d_day[order(d_day$time_minutes), ]
        n <- nrow(d_day)
        if (n < 2L) return(NULL)

        segs <- vector("list", 0L)

        for (i in seq_len(n - 1L)) {
          x1 <- d_day$time_minutes[i]
          y1 <- d_day$sensorglucose[i]
          x2 <- d_day$time_minutes[i + 1L]
          y2 <- d_day$sensorglucose[i + 1L]

          xs <- c(x1)
          ys <- c(y1)

          # Intersection with lower threshold
          if ((y1 - low_thr) * (y2 - low_thr) < 0) {
            alpha_low <- (low_thr - y1) / (y2 - y1)
            x_low <- x1 + alpha_low * (x2 - x1)
            xs <- c(xs, x_low)
            ys <- c(ys, low_thr)
          }

          # Intersection with upper threshold
          if ((y1 - high_thr) * (y2 - high_thr) < 0) {
            alpha_high <- (high_thr - y1) / (y2 - y1)
            x_high <- x1 + alpha_high * (x2 - x1)
            xs <- c(xs, x_high)
            ys <- c(ys, high_thr)
          }

          xs <- c(xs, x2)
          ys <- c(ys, y2)

          # Sort by time to avoid mis-ordered intersections
          ord <- order(xs)
          xs  <- xs[ord]
          ys  <- ys[ord]

          if (length(xs) >= 2L) {
            for (j in seq_len(length(xs) - 1L)) {
              x_start <- xs[j]
              x_end   <- xs[j + 1L]
              if (x_end <= x_start) next

              y_start <- ys[j]
              y_end   <- ys[j + 1L]
              y_mid   <- (y_start + y_end) / 2

              zone <- if (y_mid < low_thr) {
                "Below range"
              } else if (y_mid > high_thr) {
                "Above range"
              } else {
                "In range"
              }

              segs[[length(segs) + 1L]] <- data.frame(
                date_label = d_day$date_label[1],
                x_start    = x_start,
                x_end      = x_end,
                y_start    = y_start,
                y_end      = y_end,
                zone       = zone,
                stringsAsFactors = FALSE
              )
            }
          }
        }

        if (!length(segs)) return(NULL)
        do.call(rbind, segs)
      }

      # Generate segments for each day and combine them
      seg_list <- lapply(
        split(df_sid, df_sid$date_label),
        build_segments_one_day,
        low_thr  = low_thr,
        high_thr = high_thr
      )
      seg_df <- do.call(
        rbind,
        seg_list[!vapply(seg_list, is.null, logical(1L))]
      )

      band_fill  <- "#E6F4FF"
      strip_fill <- "#FFFFFF"

      tir_band <- data.frame(
        xmin = 0,
        xmax = 1440,
        ymin = low_thr,
        ymax = high_thr
      )

      if (is.null(seg_df) || !nrow(seg_df)) {
        seg_df <- data.frame(
          date_label = character(),
          x_start = numeric(), x_end = numeric(),
          y_start = numeric(), y_end = numeric(),
          zone = character()
        )
      }

      p_daily <- ggplot2::ggplot() +
        # Background band for target range
        ggplot2::geom_rect(
          data        = tir_band,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill        = band_fill,
          alpha       = 1
        ) +
        # Zone-colored segments
        ggplot2::geom_segment(
          data = seg_df,
          ggplot2::aes(
            x    = .data[["x_start"]],
            xend = .data[["x_end"]],
            y    = .data[["y_start"]],
            yend = .data[["y_end"]],
            color = .data[["zone"]]
          ),
          linewidth = 0.45
        ) +
        ggplot2::facet_wrap(~ date_label, ncol = 3) +
        ggplot2::scale_x_continuous(
          limits = c(0, 1440),
          breaks = seq(0, 1440, by = 180),
          labels = c(
            "0:00", "03:00", "06:00", "09:00",
            "12:00", "15:00", "18:00", "21:00", "24:00"
          ),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_y_continuous(
          breaks = seq(yaxis[1], yaxis[2], by = 50),
          expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(ylim = yaxis) +
        ggplot2::scale_color_manual(
          values = c(
            "In range"    = "#5EBAF2",
            "Below range" = "#EB4024",
            "Above range" = "#F9D610"
          )
        ) +
        ggplot2::labs(
          # title = paste0("Daily Glucose Profiles: ", sid),
          title = "Daily Ambulatory Glucose Profile (Raw data, no smoothing)",
          x     = "Time of day (24-hour clock)",
          y     = "Glucose (mg/dL)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          panel.border     = ggplot2::element_rect(
            color = "#CCCCCC", fill = NA, linewidth = 0.4
          ),
          axis.title       = ggplot2::element_text(size = 12, face = "bold"),
          axis.text        = ggplot2::element_text(size = 9),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(
            colour = "#EEEEEE", linewidth = 0.3
          ),
          panel.spacing.x  = grid::unit(0.8, "cm"),
          panel.spacing.y  = grid::unit(0.7, "cm"),
          legend.position  = "none",
          plot.title       = ggplot2::element_text(
            size = 16, face = "bold", hjust = 0.5
          ),
          strip.background = ggplot2::element_rect(
            fill = strip_fill, colour = "#CCCCCC", linewidth = 0.4
          ),
          strip.text       = ggplot2::element_text(size = 10, face = "bold")
        )

      # C. Output the subject-level PDF:
      #    1) Selected AGP  2) Multi-day overlay  3) Daily profile cards

      sid_chr  <- as.character(sid)
      sid_safe <- gsub("[^A-Za-z0-9_.-]", "_", sid_chr)
      pdf_file_sid <- file.path(outputdir, paste0("AGP_Report_", sid_safe, ".pdf"))

      grDevices::pdf(pdf_file_sid, width = page_width, height = page_height)

      if (!is.null(subjAGP_selected)) {
        grid::grid.draw(ggplot2::ggplotGrob(subjAGP_selected))
      }

      if (!is.null(p_multiday)) {
        grid::grid.draw(ggplot2::ggplotGrob(p_multiday))
      }

      grid::grid.draw(ggplot2::ggplotGrob(p_daily))

      grDevices::dev.off()
    }
  }

  # Actually call the per-subject plotting function if requested
  if (isTRUE(per_subject)) {
    plot_daily_subject_cards(
      aggregateAGPdata = aggregateAGPdata,
      tz        = tz,
      yaxis     = yaxis,
      low_thr   = low_thr,
      high_thr  = high_thr,
      ref_lines = ref_lines,
      outputdir = outputdir,
      agg_smooth = agg_smooth,
      bin_mins   = bin_mins,
      page_width  = page_width,
      page_height = page_height
    )
  }

  # ------------------ 8. Write a PDF with selected aggregate AGP -------------

  pdf_file <- file.path(outputdir, "AGP_Report_All.pdf")

  # Default (no grouping): choose the overall aggregate AGP
  aggAGP_selected <- switch(
    agg_smooth,
    "original" = aggAGP_original,
    "tukey"    = aggAGPtukey,
    "loess"    = aggAGPloess
  )

  # If grouping is enabled and grouped plots are available, override selection
  if (!is.null(group_var) &&
      !is.null(quartiles_grouped) &&
      nrow(quartiles_grouped) > 0L &&
      (!is.null(aggAGP_original_grouped) ||
       !is.null(aggAGPtukey_grouped) ||
       !is.null(aggAGPloess_grouped))) {

    aggAGP_selected <- switch(
      agg_smooth,
      "original" = aggAGP_original_grouped,
      "tukey"    = aggAGPtukey_grouped,
      "loess"    = aggAGPloess_grouped
    )
  }

  grDevices::pdf(pdf_file, width = page_width, height = page_height)

  # Page 1: aggregate AGP (overall or grouped)
  if (!is.null(aggAGP_selected)) {
    grid::grid.draw(ggplot2::ggplotGrob(aggAGP_selected))
  }

  # Page 2: daily overlay per subject (LOESS)
  if (isTRUE(subject_overlay) && !is.null(AGPloess)) {
    grid::grid.draw(ggplot2::ggplotGrob(AGPloess))
  }

  grDevices::dev.off()

  message("AGP Report saved to: ", pdf_file)

  invisible(list(
    data                       = aggregateAGPdata,
    quartiles                  = quartiles,
    quartiles_grouped          = quartiles_grouped,
    quartiles_grouped_original = quartiles_grouped_original,
    aggAGP_original            = aggAGP_original,
    aggAGPtukey                = aggAGPtukey,
    aggAGPloess                = aggAGPloess,
    aggAGP_original_grouped    = aggAGP_original_grouped,
    aggAGPtukey_grouped        = aggAGPtukey_grouped,
    aggAGPloess_grouped        = aggAGPloess_grouped,
    AGPloess                   = AGPloess
  ))
}

