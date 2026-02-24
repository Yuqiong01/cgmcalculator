#' Calculate Comprehensive CGM Metrics
#'
#' @description
#' This function processes Continuous Glucose Monitoring (CGM) CSV files that
#' have been pre-processed via \code{tidydata()}.
#' It computes a comprehensive set of CGM metrics.
#' The metrics are organized into several key areas, allowing users to analyze
#' and understand the glycemic control of the individual:
#' \itemize{
#'   \item Basic data and statistics
#'   \item Time-based metrics (TIR, TAR, TBR)
#'   \item Glycemic variability, fluctuations, and additional indices
#'   \item Autocorrelation and rate of change metrics
#'   \item GRADE score and risk of diabetes
#' }
#'
#' All input files must be saved as CSV files and must include three specific columns:
#' \itemize{
#'   \item The first column must contain the subject ID.
#'   \item The second column must contain the timestamp (date and time of each CGM reading).
#'   \item The third column must contain the sensor glucose values.
#'   \item The column names must be exactly \code{"subjectid"}, \code{"timestamp"},
#'   and \code{"sensorglucose"} (without quotes).
#' }
#' Files can be cleaned and formatted using this package's \code{tidydata()}
#' function to ensure they are ready for processing.
#'
#' @details
#' This function computes a wide range of CGM metrics, organized into
#' the following categories:
#'
#' \strong{1. Basic Data & Statistics}:
#' \itemize{
#'   \item Includes CGM device placement/removal dates, active monitoring time
#'   percentage, total valid days, and various glucose statistics:
#'   \itemize{
#'     \item Mean, median, minimum, maximum, range, standard deviation, interquartile range (IQR),
#'     and coefficient of variation (CV) of glucose readings.
#'   }
#' }
#'
#' \strong{2. Time in Range Percentage Metrics}:
#' \itemize{
#'   \item Calculates time spent within specific glucose ranges (Target Range Time, TIR),
#'   time spent above thresholds (Time Above Range, TAR),
#'   and time spent below thresholds (Time Below Range, TBR).
#'   \item Includes time and Area Under Curve (AUC) values for the following ranges:
#'   \itemize{
#'     \item TIR: 63--140 mg/dL, 70--140 mg/dL, 70--180 mg/dL;
#'     \item TAR: 140, 180, 250, 400 mg/dL;
#'     \item TBR: 40, 54, 63, 70 mg/dL.
#'   }
#'   \item Additional metrics include event counts, duration,
#'   and average glucose values for each range.
#'   \item These metrics can be further categorized by daytime and
#'   nighttime intervals (using suffixes \code{*_day} and \code{*_night}).
#'   \item Users can specify additional time-above-range scanning
#'   using \code{tar_seq} and time-below-range scanning using \code{tbr_seq}.
#' }
#'
#' \strong{3. Glycemic Variability, Fluctuations & Additional Metrics}:
#' \itemize{
#'   \item Includes a set of metrics related to glycemic variability and
#'   fluctuations, such as:
#'   \itemize{
#'     \item MAGE (Mean Absolute Glucose Excursion)
#'     \item MODD (Mean of Daily Differences)
#'     \item CONGA (Continuous Overlap of Gaps)
#'     \item LBGI (Low Blood Glucose Index)
#'     \item HBGI (High Blood Glucose Index)
#'     \item M-value (Glycemic Variability Index)
#'     \item J-index
#'     \item ADRR (Accumulated Deviation Risk Ratio)
#'     \item MAG (Mean Amplitude of Glycemic Excursions)
#'     \item GVP (Glycemic Variability Percentage)
#'     \item ......
#'   }
#' }
#'
#' \strong{4. Autocorrelation & Rate of Change}:
#' \itemize{
#'   \item Computes metrics related to autocorrelation (ACF)
#'   and rate of change (ROC), including:
#'   \itemize{
#'     \item Mean, variance, standard deviation, and percentiles of ROC values.
#'   }
#' }
#'
#' \strong{5. GRADE Score & Diabetes Risk}:
#' \itemize{
#'   \item Calculates the GRADE score and proportions of time spent in euglycemia,
#'   hypoglycemia, and hyperglycemia.
#'   \item Includes Glycemic Risk Index (GRI) and Continuous Glucose Monitoring Index (COGI).
#' }
#'
#' \strong{6. Daytime & Nighttime Subgroup Analysis}:
#' \itemize{
#'   \item Provides the same metrics as in the previous section,
#'   but segmented into daytime and nighttime periods.
#'   \item These subgroups are identified using suffixes \code{*_day} and \code{*_night}.
#' }
#'
#' @param inputdir Character string. Path to the directory containing the
#'   pre-processed CGM CSV files (tidydata output).
#' @param outputdir Character string. Path to the directory where the results
#'   will be written. Defaults to \code{tempdir()}.
#' @param outputname Character string. Base name (without file extension) for
#'   the output file.
#' @param hyper_dur_mins Numeric scalar. Minimum duration in minutes that sensor
#'   glucose must remain above the hyperglycemia threshold to be counted as a
#'   hyperglycemia episode. Default is 15.
#' @param hypo_dur_mins Numeric scalar. Minimum duration in minutes that sensor
#'   glucose must remain below the hypoglycemia threshold to be counted as a
#'   hypoglycemia episode. Default is 15.
#' @param hyper_end_mins Numeric scalar. Number of consecutive minutes with glucose
#'   back in the non-hyperglycemic range required to declare the end of a
#'   hyperglycemia episode. Short transient returns to normal shorter than this
#'   value are bridged and treated as part of the same episode. Default is 15.
#' @param hypo_end_mins Numeric scalar. Number of consecutive minutes with glucose
#'   back in the non-hypoglycemic range required to declare the end of a
#'   hypoglycemia episode. Helps avoid over-splitting episodes due to brief
#'   returns to normal. Default is 15.
#' @param ext_hyper_mins Numeric scalar. Minimum total duration in minutes of a
#'   single hyperglycemia episode for it to be additionally classified as
#'   extended hyperglycemia. Default is 120.
#' @param ext_hypo_mins Numeric scalar. Minimum total duration in minutes of a
#'   single hypoglycemia episode for it to be additionally classified as
#'   extended hypoglycemia. Default is 120.
#' @param mage_sd Numeric scalar. SD multiplier used when computing MAGE
#'   (mean amplitude of glycemic excursions) and related indices (e.g. MGE/MGN),
#'   defined relative to \code{mage_sd * SD}. Default is 1.
#' @param conga_hrs Numeric scalar. Time interval in hours used for the CONGA
#'   metric (continuous overlapping net glycemic action).
#' @param modd_lag Integer scalar. Lag in days used for MODD (mean of absolute
#'   differences in glucose measured at the same time of day on days
#'   \code{d} and \code{d + lag}). Default is 1.
#' @param roc_lag Numeric scalar. Time window in minutes used when
#'   calculating rate-of-change (ROC) metrics. Default is 15.
#' @param ac_lag Number of lags used for AC-mean and AC-variance (default 30).
#' @param gri_targets Numeric vector of length 4. Glucose cut-points (in mg/dL)
#'   used to define very low, low, high, and very high ranges for the GRI
#'   metric, in the order \code{c(vlow_cut, low_cut, high_cut, vhigh_cut)}.
#' @param gri_weights Numeric vector of length 4. Weights applied to the
#'   very-low, low, high, and very-high ranges when computing GRI, in the same
#'   order as \code{gri_targets}. Default is \code{c(3.0, 2.4, 0.8, 1.6)}.
#' @param cogi_targets Numeric vector of length 2. Lower and upper glucose
#'   targets (in mg/dL) used to define TIR/TBR for COGI calculations.
#' @param cogi_weights Numeric vector of length 3. Weights assigned to the
#'   three COGI components (time in range, time out of range, and glycemic
#'   variability), in that order. Must sum to 1. Default is \code{c(0.5, 0.35, 0.15)}.
#' @param unit Character string. Glucose unit of the input data. If
#'   \code{"mg/dL"}, values are used as is; any other value is treated as
#'   mmol/L and multiplied by 18 to convert to mg/dL.
#' @param ndays Numeric scalar. Theoretical maximum wear period in days
#'   (e.g. device prescription length). Used as the denominator for
#'   \code{active_percent}. If \code{ndays} is smaller than the actual observed
#'   number of days, the actual value is used to ensure
#'   \code{active_percent <= 100}. Default is 14.
#' @param tz Character string. Time zone used to interpret the timestamp
#'   column (e.g. \code{"UTC"}, \code{"Asia/Shanghai"}). All date-based
#'   operations (day/night classification, daily aggregation, MODD,
#'   intraday metrics) are performed using this time zone.
#'   Default is \code{"UTC"}.
#' @param daily_avg Logical. If \code{TRUE}, event counts, event minutes, and
#'   AUC-related variables are expressed as averages per observed day using
#'   \code{num_days} as the denominator. If \code{FALSE}, total values over the
#'   full observation period are returned. Default is \code{TRUE}.
#' @param wear_cut Numeric scalar between 0 and 1. Minimum daily wear fraction
#'   used to define \code{num_days_good_data}. Days with a lower wear
#'   proportion are not counted as full valid days.
#' @param tbr_seq Optional numeric vector. Sequence of lower glucose thresholds
#'   (in mg/dL) at which TBR (time below range) will be scanned and
#'   summarized. If \code{NULL}, only default thresholds are used.
#' @param tar_seq Optional numeric vector. Sequence of upper glucose thresholds
#'   (in mg/dL) at which TAR (time above range) will be scanned and
#'   summarized. If \code{NULL}, only default thresholds are used.
#' @param calc_day_night Integer vector of length 2. Start and end hours (0--23)
#'   defining the daytime period; hours outside this range are treated as
#'   nighttime when a \code{wake} variable is not available.
#' @param id_filename Logical. If \code{TRUE}, use the file name (without
#'   extension) as subject ID. If \code{FALSE}, use the ID contained in the
#'   first row/column of the data. Default is \code{TRUE}.
#' @param printname Logical. If \code{TRUE}, print the name of each file as it
#'   is processed. Default is \code{FALSE}.
#' @param format Character string. Output format. \code{"rows"} = metrics in rows,
#'   subjects in columns. \code{"long"} = subjects in rows and metrics in columns.
#'
#' @return
#' A data.frame containing all CGM metrics.
#' If \code{format = "rows"}, metrics are in rows and subjects in columns.
#' The first column is \code{"subjectid"} (metric name), and
#' each additional column corresponds to one subject.
#'
#' If \code{format = "long"}, subjects are in rows and metrics in columns.
#' The first column is \code{"subjectid"} (subject ID), and
#' additional columns correspond to the computed metrics.
#'
#' @examples
#' \dontrun{
#' # Example: compute CGM metrics from example tidydata output
#' inputdir  <- system.file("exdata", "tidydata", package = "cgmcalculator")
#' outputdir <- tempdir()
#'
#' res <- metricsstats(
#'   inputdir   = inputdir,
#'   outputdir  = outputdir,
#'   outputname = "summary_metrics_example",
#'   format     = "long"
#' )
#'
#' # View first few rows of the resulting metrics table
#' head(res)
#' }
#' @export
metricsstats <- function(inputdir,
                         outputdir = tempdir(),
                         outputname = "summary_metrics",
                         hyper_dur_mins = 15,
                         hypo_dur_mins = 15,
                         hyper_end_mins = 15,
                         hypo_end_mins  = 15,
                         ext_hypo_mins  = 120,
                         ext_hyper_mins = 120,
                         mage_sd = 1,
                         conga_hrs = 1,
                         modd_lag = 1,
                         roc_lag = 15,
                         ac_lag = 30,
                         gri_targets = c(54, 70, 180, 250),
                         gri_weights = c(3.0, 2.4, 0.8, 1.6),
                         cogi_targets = c(70, 180),
                         cogi_weights = c(.5, .35, .15),
                         unit = "mg/dL",
                         ndays = 14,
                         tz = "UTC",
                         daily_avg = TRUE,
                         wear_cut = NULL,
                         tbr_seq = NULL,
                         tar_seq = NULL,
                         calc_day_night = NULL,
                         id_filename = TRUE,
                         printname = FALSE,
                         format = "rows") {

  # If calc_day_night is specified, extract daystart and dayend
  if (!is.null(calc_day_night)) {
    if (length(calc_day_night) != 2) {
      stop("calc_day_night must be a vector with two values: daystart and dayend.")
    }
    daystart <- calc_day_night[1]
    dayend   <- calc_day_night[2]
  } else {
    daystart <- 6
    dayend   <- 0
  }

  # Get all files in the input directory
  files <- base::list.files(path = inputdir, full.names = TRUE)
  if (!length(files)) stop("No files found in inputdir.")

  # Unified decoding of GRI thresholds: gri_targets = c(vlow, low, high, vhigh)
  cuts <- as.numeric(gri_targets)
  if (length(cuts) != 4L)
    stop("gri_targets must be a numeric vector of length 4: c(vlow, low, high, vhigh).")
  if (any(!is.finite(cuts)))
    stop("gri_targets must be finite numeric values.")

  # Set the GRI thresholds
  vlow_cut_gri  <- cuts[1]
  low_cut_gri   <- cuts[2]
  high_cut_gri  <- cuts[3]
  vhigh_cut_gri <- cuts[4]

  # Process files one by one and store results in a list

  per_file_list <- vector("list", length(files))
  file_ids <- character(length(files))

  for (f in seq_along(files)) {
    tb <- .read_and_prepare(files[f], id_filename, unit, tz = tz)
    subj_id <- tb$subjectid[1]
    file_ids[f] <- subj_id

    res <- list()
    res[["subject_id"]] <- subj_id

    # Interval & wear time
    interval <- .median_interval_secs(tb$timestamp)
    wear <- .block_wear_and_window(tb, interval, ndays = ndays)
    res <- c(res, wear)

    # Strictly define num_days_good_data
    if (!is.null(wear_cut)) {
      nd_strict <- .strict_good_days(tb, interval, wear_cut)
      res[["num_days_good_data"]] <- as.numeric(nd_strict)
    } else if (is.null(res[["num_days_good_data"]])) {
      res[["num_days_good_data"]] <- as.numeric(res[["num_days"]])
    }

    # Distribution statistics
    res <- c(res, .block_distribution(tb))

    # Overall GRADE score
    res <- c(res, .block_grade(tb$sensorglucose))

    # GRI (Glucose Risk Index)
    res[["gri"]] <- .gri_from_sg(
      tb$sensorglucose, interval,
      vlow_cut  = vlow_cut_gri,
      low_cut   = low_cut_gri,
      high_cut  = high_cut_gri,
      vhigh_cut = vhigh_cut_gri,
      w_vlow  = gri_weights[1],
      w_low   = gri_weights[2],
      w_high  = gri_weights[3],
      w_vhigh = gri_weights[4]
    )

    # COGI (Continuous Glucose Monitoring Index)
    res[["cogi"]] <- .cogi_from_sg(
      sg            = tb$sensorglucose,
      interval_secs = interval,
      targets       = cogi_targets,
      weights       = cogi_weights
    )

    # Common thresholds (fixed output, including day/night)
    res <- c(res, .block_thresholds_defaults(
      tb, interval,
      above_len_mins = hyper_dur_mins,
      below_len_mins = hypo_dur_mins,
      hyper_end_mins      = hyper_end_mins,
      hypo_end_mins       = hypo_end_mins,
      ext_hyper_mins = ext_hyper_mins,
      ext_hypo_mins  = ext_hypo_mins
    ))

    # Day/night splitting (optional)
    if (!is.null(calc_day_night)) {
      res <- c(res, .block_day_night_defaults(
        tb, interval, daystart, dayend,
        hypo_dur_mins, hyper_dur_mins,
        gri_weights = gri_weights,
        gri_targets = gri_targets
      ))
    }

    # Variability & ROC (rate of change)
    res <- c(
      res,
      .block_variability(
        tb           = tb,
        interval     = interval,
        mage_sd      = mage_sd,
        conga_hrs    = conga_hrs,
        modd_lag     = modd_lag,
        ac_lag       = ac_lag,
        roc_lag  = roc_lag
      )
    )

    # Scanned thresholds (optional)
    if (!is.null(tbr_seq) || !is.null(tar_seq)) {
      scan <- .block_scanned_intervals_optional(
        tb, interval,
        tbr_seq = tbr_seq, tar_seq = tar_seq,
        above_len_mins = hyper_dur_mins,
        below_len_mins = hypo_dur_mins,
        calc_day_night = calc_day_night,
        daystart = daystart,
        dayend   = dayend,
        hypo_dur_mins = hypo_dur_mins,
        hyper_dur_mins = hyper_dur_mins
      )
      res <- c(res, scan)
    }

    # --- add binary flags ---
    res <- .add_event_binary(res)

    # Daily averaging
    denom_days <- res[["num_days_wear"]]
    if (is.null(denom_days) || !is.finite(as.numeric(denom_days)) || as.numeric(denom_days) <= 0) {
      denom_days <- res[["num_days"]]
    }

    res <- .dailyize_metrics(
      res_list   = res,
      num_days   = denom_days,
      daily_avg  = daily_avg
    )

    per_file_list[[f]] <- res
  }

  # Combine results into a wide table
  all_keys <- unique(unlist(lapply(per_file_list, names), use.names = FALSE))
  all_keys <- all_keys[nzchar(all_keys)]

  out <- data.frame(
    "subjectid" = all_keys,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (i in seq_along(per_file_list)) {
    vals <- per_file_list[[i]]
    col  <- rep(NA_character_, length(all_keys))
    hit  <- match(names(vals), all_keys)
    col[hit] <- as.character(unlist(vals, use.names = FALSE))
    out[[file_ids[i]]] <- col
  }

  # Remove the row with subject_id
  out <- out[out$subjectid != "subject_id", , drop = FALSE]

  # If format is "rows", return as wide table
  if (identical(format, "rows")) {

    final_df <- out

    # If format is "long", transpose the data
  } else if (identical(format, "long")) {

    wide_mat <- t(out[, -1, drop = FALSE])
    wide_df  <- as.data.frame(wide_mat,
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

    colnames(wide_df) <- out$subjectid

    final_df <- data.frame(
      subjectid = rownames(wide_df),
      wide_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(final_df) <- NULL

    # Raise error if format is not recognized
  } else {
    stop('format must be "rows" or "long".')
  }

  # Save the results to a CSV file
  filename <- file.path(outputdir, paste0(outputname, ".csv"))
  utils::write.csv(final_df, file = filename, row.names = FALSE, na = "")
  invisible(final_df)
}

#' Generate a sequence of low-glucose thresholds for TBR scanning
#'
#' @description
#' Convenience wrapper around \code{seq()} to generate a sequence of
#' lower glucose thresholds (in mg/dL) for time-below-range (TBR) scanning.
#'
#' @param from Numeric. Lowest threshold (mg/dL) in the sequence.
#' @param to Numeric. Highest threshold (mg/dL) in the sequence.
#' @param by Numeric. Step size between consecutive thresholds (mg/dL).
#'
#' @return A numeric vector of glucose thresholds (mg/dL).
#' @export
tbr_seq <- function(from = 40, to = 80, by = 1) {
  if (by <= 0) stop("by must be > 0")
  seq(from = from, to = to, by = by)
}

#' Generate a sequence of high-glucose thresholds for TAR scanning
#'
#' @description
#' Convenience wrapper around \code{seq()} to generate a sequence of
#' upper glucose thresholds (in mg/dL) for time-above-range (TAR) scanning.
#'
#' @param from Numeric. Lowest threshold (mg/dL) in the sequence.
#' @param to Numeric. Highest threshold (mg/dL) in the sequence.
#' @param by Numeric. Step size between consecutive thresholds (mg/dL).
#'
#' @return A numeric vector of glucose thresholds (mg/dL).
#' @export
tar_seq <- function(from = 120, to = 400, by = 10) {
  if (by <= 0) stop("by must be > 0")
  seq(from = from, to = to, by = by)
}

## ===================== 2) Input/Output & Basic Utilities ====================
##  reading, intervals, wear-time, generic helpers)

.read_and_prepare <- function(file_path, id_filename, unit, tz = "UTC") {
  tb <- utils::read.csv(file_path, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  colnames(tb) <- tolower(colnames(tb))

  if (!all(c("subjectid","timestamp","sensorglucose") %in% names(tb))) {
    stop("Input csv must have columns: subjectid, timestamp, sensorglucose")
  }

  if (identical(id_filename, FALSE)) {
    tb$subjectid <- tb$subjectid[1]
  } else {
    # tb$subjectid[1] <- sub("[.]csv$", "", basename(file_path))
    tb$subjectid <- sub("[.]csv$", "", basename(file_path))
  }

  tb <- unique(tb)
  tb$timestamp <- parsedate::parse_date(
    tb$timestamp,
    approx = FALSE
  )
  attr(tb$timestamp, "tzone") <- tz

  tb$sensorglucose[tb$sensorglucose == "Low"]  <- 40
  tb$sensorglucose[tb$sensorglucose == "High"] <- 400

  if (!identical(unit, "mg/dL")) {
    tb$sensorglucose <- suppressWarnings(as.numeric(tb$sensorglucose) * 18)
  } else {
    tb$sensorglucose <- suppressWarnings(as.numeric(tb$sensorglucose))
  }

  tb <- tb[!is.na(tb$timestamp) & !is.na(tb$sensorglucose), ]
  tb
}

.median_interval_secs <- function(timestamp_vec) {
  interval <- round(stats::median(diff(as.numeric(timestamp_vec)), na.rm = TRUE), 0)
  abs(interval)
}

.block_wear_and_window <- function(tb, interval, ndays = 14) {
  start_ts <- min(tb$timestamp, na.rm = TRUE)
  end_ts   <- max(tb$timestamp, na.rm = TRUE)

  ap <- active_percent(tb, interval = interval, ndays = ndays)

  list(
    date_cgm_placement   = as.character(start_ts),
    date_cgm_remove      = as.character(end_ts),
    active_percent       = ap$active_percent,
    num_days             = ap$num_days,
    num_days_good_data   = ap$num_days,
    total_glu_readings   = ap$total_glu_readings
  )
}

.effective_days <- function(tb, interval) {
  per_day <- split(tb, as.Date(tb$timestamp))
  theory  <- 86400 / interval
  if (!is.finite(theory) || theory <= 0) return(NA_real_)
  rates <- sapply(per_day, function(df) {
    sum(!is.na(df$sensorglucose)) / theory
  })
  sum(pmin(rates, 1), na.rm = TRUE)
}

# Count valid readings (non-NA, finite)
.n_valid <- function(x) sum(is.finite(as.numeric(x)))

# wear days based on number of valid readings (additive)
.wear_days_from_n <- function(n_valid, interval_secs) {
  if (!is.finite(interval_secs) || interval_secs <= 0) return(NA_real_)
  as.numeric(n_valid) * (interval_secs / 86400)
}

# Effective days within a period (day or night), based on daystart/dayend hours
.effective_days_period <- function(ts, sg, interval, daystart, dayend, tag = c("day","night")) {
  tag <- match.arg(tag)

  if (!length(ts) || length(ts) != length(sg)) return(NA_real_)
  if (!is.finite(interval) || interval <= 0) return(NA_real_)

  hr <- as.integer(format(ts, "%H"))

  hours_in_day <-
    if (daystart == dayend) {
      0:23
    } else if (daystart < dayend) {
      daystart:dayend
    } else {
      c(daystart:23, 0:dayend)
    }

  is_day <- hr %in% hours_in_day
  keep <- if (tag == "day") is_day else !is_day

  ts2 <- ts[keep]
  sg2 <- sg[keep]

  if (!length(ts2)) return(0)

  per_day <- split(seq_along(ts2), as.Date(ts2))

  # theoretical points per day within this period
  n_hours <- if (daystart == dayend) 24 else length(unique(hours_in_day))
  theory_pts <- (n_hours * 3600) / interval
  if (!is.finite(theory_pts) || theory_pts <= 0) return(NA_real_)

  rates <- vapply(per_day, function(idx) {
    .n_valid(sg2[idx]) / theory_pts
  }, numeric(1))

  sum(pmin(rates, 1), na.rm = TRUE)
}

.strict_good_days <- function(tb, interval, wear_cut) {
  if (is.null(wear_cut) || is.na(wear_cut) || wear_cut <= 0 || wear_cut > 1) return(NA_real_)
  per_day <- split(tb, as.Date(tb$timestamp))
  theory <- 86400 / interval
  rates <- sapply(per_day, function(df) length(which(!is.na(df$sensorglucose))) / theory)
  sum(rates >= wear_cut, na.rm = TRUE)
}

# active_percent
active_percent <- function(tb, interval = NULL, ndays = 14) {
  if (!("timestamp" %in% names(tb)) || !("sensorglucose" %in% names(tb))) {
    stop("tb must contain 'timestamp' and 'sensorglucose'.")
  }

  if (is.null(interval)) {
    interval <- .median_interval_secs(tb$timestamp)
  }

  num_days_eff <- .effective_days(tb, interval)
  total_glu_readings <- sum(!is.na(tb$sensorglucose))

  if (!is.finite(num_days_eff) || num_days_eff <= 0) {
    return(list(
      active_percent     = NA_real_,
      num_days           = NA_real_,
      total_glu_readings = total_glu_readings
    ))
  }

  ndays_use <- ndays
  if (ndays_use < num_days_eff) {
    warning(
      sprintf("ndays (%.3f) < effective days (%.3f); using effective days instead.",
              ndays_use, num_days_eff)
    )
    ndays_use <- num_days_eff
  }

  theory_pts_per_day <- 86400 / interval
  denom_pts <- ndays_use * theory_pts_per_day

  ap <- 100 * total_glu_readings / denom_pts
  ap <- max(0, min(100, ap))

  list(
    active_percent     = as.numeric(ap),
    num_days           = as.numeric(round(num_days_eff, 3)),
    total_glu_readings = as.numeric(total_glu_readings)
  )
}

# Supports expanding the c(from, to, by) triplet into a sequence
._expand_seq <- function(v) {
  v <- as.numeric(v)
  if (length(v) == 3 && is.finite(v[1]) && is.finite(v[2]) && is.finite(v[3]) &&
      v[3] > 0 && v[1] <= v[2]) {
    return(seq(from = v[1], to = v[2], by = v[3]))
  }
  unique(sort(v))
}

# .percent_time <- function(flag_vec, interval_secs) {
#   ((sum(flag_vec, na.rm = TRUE) * (interval_secs / 60)) /
#      (length(flag_vec) * (interval_secs / 60))) * 100
# }

.percent_time <- function(flag_vec, interval_secs) {
  denom <- sum(!is.na(flag_vec))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  100 * sum(flag_vec == 1, na.rm = TRUE) / denom
}

## ===================== 3) Summary / Distribution Metrics ====================

.block_distribution <- function(tb) {
  sg <- tb$sensorglucose
  qs <- stats::quantile(sg, probs = c(0.25, 0.75), na.rm = TRUE)
  q1 <- qs[1]
  q3 <- qs[2]

  list(
    mean_glu   = mean(sg, na.rm = TRUE),
    median_glu = stats::median(sg, na.rm = TRUE),
    gmi        = 3.31 + (0.02392 * mean(sg, na.rm = TRUE)),
    ehba1c     = (46.7 + mean(sg, na.rm = TRUE)) / 28.7,
    p25_glu    = q1,
    p75_glu    = q3,
    iqr_glu    = q3 - q1,
    cqv_glu = if ((q3 + q1) == 0 || is.na(q3) || is.na(q1)) NA_real_
    else ((q3 - q1) / (q3 + q1)) * 100,
    sd_glu     = stats::sd(sg, na.rm = TRUE),
    cv_glu     = (stats::sd(sg, na.rm = TRUE) / base::mean(sg, na.rm = TRUE)) * 100,
    min_glu    = min(sg, na.rm = TRUE),
    max_glu    = max(sg, na.rm = TRUE),
    range_glu  = max(sg, na.rm = TRUE) - min(sg, na.rm = TRUE)
  )
}

## ======================== 4) Risk Indices / Ratings =========================
## GRADE, GRI, LBGI/HBGI, indices, COGI, etc.

# GRADE: GRADE + %hypo/%eugly/%hyper
# https://doi.org/10.1111/j.1464-5491.2007.02119.x
.block_grade <- function(sg,
                         hypo_cut = 70,
                         hyper_cut = 140,
                         min_mg = 37,
                         max_mg = 630) {

  x <- as.numeric(sg)
  out <- list(
    grade_score      = NA_real_,
    grade_hypo_prop  = NA_real_,
    grade_eugly_prop = NA_real_,
    grade_hyper_prop = NA_real_
  )
  if (!length(x) || all(is.na(x))) return(out)

  ok <- is.finite(x)
  g  <- rep(NA_real_, length(x))

  in_range <- ok & x >= min_mg & x <= max_mg
  if (any(in_range)) {
    mmol <- x[in_range] / 18
    safe <- mmol > 0
    tmp  <- rep(NA_real_, length(mmol))
    tmp[safe] <- 425 * (log10(log10(mmol[safe])) + 0.16)^2
    tmp[tmp > 50] <- 50
    g[in_range]   <- tmp
  }

  out_range <- ok & (x < min_mg | x > max_mg)
  g[out_range] <- 50

  if (all(is.na(g))) return(out)

  grade_mean  <- mean(g, na.rm = TRUE)
  grade_total <- sum(g,  na.rm = TRUE)
  if (!is.finite(grade_total) || grade_total <= 0) {
    out$grade_score <- grade_mean
    return(out)
  }

  hypo_idx  <- ok & x < hypo_cut
  eugly_idx <- ok & x >= hypo_cut & x <= hyper_cut
  hyper_idx <- ok & x > hyper_cut

  hypo_sum   <- sum(g[hypo_idx],   na.rm = TRUE)
  eugly_sum  <- sum(g[eugly_idx],  na.rm = TRUE)
  hyper_sum  <- sum(g[hyper_idx],  na.rm = TRUE)

  out$grade_score      <- as.numeric(grade_mean)
  out$grade_hypo_prop  <- as.numeric(100 * hypo_sum  / grade_total)
  out$grade_eugly_prop <- as.numeric(100 * eugly_sum / grade_total)
  out$grade_hyper_prop <- as.numeric(100 * hyper_sum / grade_total)

  out
}

# GRI: Derived from GRADE
# doi: 10.1177/19322968221085273
.gri_from_sg <- function(sg,
                         interval_secs,
                         vlow_cut  = 54,
                         low_cut   = 70,
                         high_cut  = 180,
                         vhigh_cut = 250,
                         w_vlow  = 3.0,
                         w_low   = 2.4,
                         w_high  = 0.8,
                         w_vhigh= 1.6) {

  sg <- as.numeric(sg)
  if (!length(sg) || all(is.na(sg))) return(NA_real_)

  vlow  <- .percent_time(ifelse(sg < vlow_cut,                    1, 0), interval_secs)
  low   <- .percent_time(ifelse(sg >= vlow_cut & sg < low_cut,    1, 0), interval_secs)
  high  <- .percent_time(ifelse(sg >  high_cut & sg <= vhigh_cut, 1, 0), interval_secs)
  vhigh <- .percent_time(ifelse(sg >  vhigh_cut,                  1, 0), interval_secs)

  gri <- w_vlow * vlow + w_low * low + w_high * high + w_vhigh * vhigh
  if (!is.finite(gri)) return(NA_real_)
  ifelse(gri > 100, 100, gri)
}

.lbgi_hbgi_from_sg <- function(sg) {
  sg <- as.numeric(sg)
  if (!length(sg) || all(is.na(sg))) return(list(lbgi = NA_real_, hbgi = NA_real_))
  gluctransform <- (log(sg)^1.084) - 5.381
  rBG <- 22.77 * (gluctransform^2)
  rl  <- ifelse(gluctransform <= 0, rBG, 0)
  rh  <- ifelse(gluctransform >  0, rBG, 0)
  list(
    lbgi = mean(stats::na.omit(rl)),
    hbgi = mean(stats::na.omit(rh))
  )
}

# M-value
# doi: 10.1111/j.0954-6820.1965.tb01810.x
.m_value_from_sg <- function(sg, r = 90) {
  x <- as.numeric(sg)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(1000 * abs(log10(x / r))^3, na.rm = TRUE)
}

# MAD (Median Absolute Deviation)
# https://doi.org/10.1371/ journal.pone.0248560
.mad_from_sg <- function(sg, constant = 1.4826) {
  x <- as.numeric(sg)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::mad(x, constant = constant, na.rm = TRUE)
}

# GVP: Glucose Variability Percentage
# doi: 10.1089/dia.2017.0187
.gvp_from_sg <- function(sg, interval_secs) {
  v <- as.numeric(sg)
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(NA_real_)

  if (!is.finite(interval_secs) || interval_secs <= 0) return(NA_real_)
  dt0 <- interval_secs / 60

  diffvec <- diff(v)
  if (!length(diffvec)) return(NA_real_)

  added_length <- sqrt(dt0^2 + diffvec^2)
  base_length  <- length(diffvec) * dt0

  if (!is.finite(base_length) || base_length <= 0) return(NA_real_)

  as.numeric((sum(added_length, na.rm = TRUE) / base_length - 1) * 100)
}

# AC_Mean / AC_Var from ACF
# https://doi.org/10.1038/s43856-025-00819-5
.ac_summary_from_sg <- function(sg, ac_lag = 30) {
  x <- as.numeric(sg)
  x <- x[is.finite(x)]

  out <- list(ac_mean = NA_real_, ac_var = NA_real_)
  if (length(x) < (ac_lag + 1L)) return(out)

  ac <- stats::acf(x, lag.max = ac_lag, plot = FALSE)$acf
  ac_lags <- ac[2:(ac_lag + 1L)]

  out$ac_mean <- mean(ac_lags, na.rm = TRUE)
  out$ac_var  <- stats::var(ac_lags,  na.rm = TRUE)
  out
}

# Hyperglycemia Index
# doi: 10.1089/dia.2008.0132
.hyper_index_from_sg <- function(sg,
                                 ULTR = 140,
                                 a    = 1.1,
                                 c    = 30) {
  x <- as.numeric(sg)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)

  n_all <- length(x)
  if (!n_all) return(NA_real_)

  hi_vals <- x[x > ULTR]
  if (!length(hi_vals)) return(0)

  num <- sum((hi_vals - ULTR)^a, na.rm = TRUE)
  if (!is.finite(num)) return(NA_real_)

  as.numeric(num / (n_all * c))
}

# Hypoglycemia Index
# doi: 10.1089/dia.2008.0132
.hypo_index_from_sg <- function(sg,
                                LLTR = 80,
                                b    = 2.0,
                                d    = 30) {
  x <- as.numeric(sg)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)

  n_all <- length(x)
  if (!n_all) return(NA_real_)

  lo_vals <- x[x < LLTR]
  if (!length(lo_vals)) return(0)

  num <- sum((LLTR - lo_vals)^b, na.rm = TRUE)
  if (!is.finite(num)) return(NA_real_)

  as.numeric(num / (n_all * d))
}

# COGI (Continuous Glucose Monitoring Index)
# https://doi.org/10.1177/1932296819838525
.cogi_from_sg <- function(sg,
                          interval_secs,
                          targets = c(70, 180),
                          weights = c(.5, .35, .15)) {
  sg <- as.numeric(sg)
  sg <- sg[is.finite(sg)]
  if (!length(sg)) return(NA_real_)

  targets <- sort(as.numeric(targets))
  if (length(targets) != 2L || any(!is.finite(targets))) {
    stop("cogi_targets must be a numeric vector of length 2.")
  }

  if (is.null(interval_secs) || !is.finite(interval_secs) || interval_secs <= 0) {
    interval_secs <- 300
  }

  tir_flag <- ifelse(sg >= targets[1] & sg <= targets[2], 1, 0)
  ir <- .percent_time(tir_flag, interval_secs)

  br_flag <- ifelse(sg < targets[1], 1, 0)
  br <- .percent_time(br_flag, interval_secs)

  stddev <- stats::sd(sg, na.rm = TRUE)

  weight_features <- function(feature, scale_range, weight = 1, increasing = FALSE) {
    lo <- min(scale_range)
    hi <- max(scale_range)
    if (increasing) {
      out <- as.numeric(feature > lo) * (feature - lo) / (hi - lo)
    } else {
      out <- as.numeric(feature < hi) * (feature - hi) / (lo - hi)
    }
    out <- pmin(out, 1)
    out * weight
  }

  wf_ir  <- weight_features(ir,     c(0, 100), weight = weights[1], increasing = TRUE)
  wf_br  <- weight_features(br,     c(0, 15),  weight = weights[2], increasing = FALSE)
  wf_sd  <- weight_features(stddev, c(18,108), weight = weights[3], increasing = FALSE)

  cogi_score <- (wf_ir + wf_br + wf_sd) * 100
  as.numeric(cogi_score)
}

## ============================ 5) Episode Tools ==============================

# Calculate episode (with end condition + extended)
# https://doi.org/10.1016/ S2213-8587(22)00319-9
.episode_stats_with_end <- function(values,
                                    threshold,
                                    interval_secs,
                                    direction = c("over","under"),
                                    min_len_mins,
                                    end_len_mins,
                                    extended_mins = 120) {
  direction <- match.arg(direction)
  x <- as.numeric(values)
  x[!is.finite(x)] <- NA
  # x <- x[!is.na(x)]
  n <- length(x)
  if (!n) {
    return(list(
      n            = 0L,
      total_mins   = 0,
      mean_dur     = NA_real_,
      mean_glu     = NA_real_,
      n_ext        = 0L,
      total_mins_ext = 0,
      mean_dur_ext = NA_real_,
      mean_glu_ext = NA_real_
    ))
  }

  if (direction == "over") {
    flag <- x > threshold
  } else {
    flag <- x < threshold
  }

  # min_len_pts <- (min_len_mins  * 60) / interval_secs
  # end_len_pts <- (end_len_mins  * 60) / interval_secs
  # ext_len_pts <- (extended_mins * 60) / interval_secs

  min_len_pts <- ceiling((min_len_mins * 60) / interval_secs)
  end_len_pts <- ceiling((end_len_mins * 60) / interval_secs)
  ext_len_pts <- ceiling((extended_mins * 60) / interval_secs)

  events_idx <- vector("list", 0L)
  i <- 1L

  while (i <= n) {
    while (i <= n && !flag[i]) i <- i + 1L
    if (i > n) break

    j <- i
    while (j <= n && flag[j]) j <- j + 1L
    core_start <- i
    core_end   <- j - 1L
    core_len_pts <- core_end - core_start + 1L

    if (core_len_pts <= min_len_pts) {
      i <- j
      next
    }

    last_end <- core_end
    k <- j
    while (k <= n) {
      norm_start <- k
      while (k <= n && !flag[k]) k <- k + 1L
      norm_end <- k - 1L

      if (norm_start <= norm_end) {
        norm_len_pts <- norm_end - norm_start + 1L
        if (norm_len_pts >= end_len_pts) {
          break
        } else {
          last_end <- norm_end
        }
      }

      if (k > n) break
      true_start <- k
      while (k <= n && flag[k]) k <- k + 1L
      true_end <- k - 1L
      last_end <- true_end
    }

    ep_start <- core_start
    ep_end   <- last_end
    events_idx[[length(events_idx) + 1L]] <- ep_start:ep_end

    i <- max(ep_end + 1L, i + 1L)
  }

  if (!length(events_idx)) {
    return(list(
      n            = 0L,
      total_mins   = 0,
      mean_dur     = NA_real_,
      mean_glu     = NA_real_,
      n_ext        = 0L,
      total_mins_ext = 0,
      mean_dur_ext = NA_real_,
      mean_glu_ext = NA_real_
    ))
  }

  ep_dur <- vapply(events_idx,
                   function(idx) length(idx) * interval_secs / 60,
                   numeric(1))
  ep_glu <- vapply(events_idx,
                   function(idx) mean(x[idx], na.rm = TRUE),
                   numeric(1))

  total_mins <- sum(ep_dur)

  is_ext <- ep_dur >= extended_mins
  n_ext  <- sum(is_ext)
  total_mins_ext <- sum(ep_dur[is_ext])

  list(
    n            = length(ep_dur),
    total_mins   = total_mins,
    mean_dur     = mean(ep_dur),
    mean_glu     = mean(ep_glu),
    n_ext        = n_ext,
    total_mins_ext = total_mins_ext,
    mean_dur_ext = if (n_ext) mean(ep_dur[is_ext]) else NA_real_,
    mean_glu_ext = if (n_ext) mean(ep_glu[is_ext]) else NA_real_
  )
}

## ======================== 6) Fixed Threshold Blocks =========================
##  Full Day TIR/TBR/TAR + AUC + Episode

# AUC helpers: Excess area for upper/lower threshold (mg·h/dL)
# https://doi.org/10.2337/dc17-1600
.auc_excess <- function(values, threshold, interval_secs, type = c("over","under")) {
  type <- match.arg(type)
  y <- as.numeric(values)
  if (type == "over") {
    y <- pmax(y - threshold, 0)
  } else {
    y <- pmax(threshold - y, 0)
  }
  if (all(!is.finite(y)) || sum(y, na.rm = TRUE) == 0) return(0)

  x <- seq(0, by = interval_secs / (60*60), length.out = length(y))
  auc <- tryCatch(pracma::trapz(x, y), error = function(e) NA_real_)
  if (!is.finite(auc)) 0 else auc
}

# AUC helper: Area within (TIR) range
.auc_inrange <- function(values, lo, hi, interval_secs, mode = c("clipped","excess_from_lo")) {
  mode <- match.arg(mode)
  y <- as.numeric(values)
  if (!length(y)) return(0)

  if (mode == "clipped") {
    y <- pmin(pmax(y, lo), hi)
  } else {
    y <- pmin(pmax(y - lo, 0), hi - lo)
  }

  if (all(!is.finite(y)) || sum(y, na.rm = TRUE) == 0) return(0)

  x <- seq(0, by = interval_secs / 3600, length.out = length(y))
  auc <- tryCatch(pracma::trapz(x, y), error = function(e) NA_real_)
  if (!is.finite(auc)) 0 else auc
}

# Fixed output thresholds (Total, without scanning)
.block_thresholds_defaults <- function(tb,
                                       interval,
                                       above_len_mins,
                                       below_len_mins,
                                       hyper_end_mins      = above_len_mins,
                                       hypo_end_mins       = below_len_mins,
                                       ext_hyper_mins = 120,
                                       ext_hypo_mins  = 120) {
  sg <- tb$sensorglucose
  res <- list()

  # Hyperglycemia threshold list (TAR)
  over_list  <- c(140, 180, 250, 400)

  # Hypoglycemia threshold list (TBR)
  under_list <- c(40, 54, 63, 70)

  # Hyperglycemia: TAR + AUC + episode (Guideline version)
  for (th in over_list) {

    # Percent time & AUC (keeping original logic)
    res[[paste0("tar_prop_",  th)]] <-
      .percent_time(ifelse(sg > th, 1, 0), interval)
    res[[paste0("auc_above_", th)]] <-
      .auc_excess(sg, th, interval, type = "over")

    # episode: Unified using episode_stats_with_end()
    ep_h <- .episode_stats_with_end(
      values        = tb$sensorglucose,
      threshold     = th,
      interval_secs = interval,
      direction     = "over",
      min_len_mins  = above_len_mins,
      end_len_mins  = hyper_end_mins,
      extended_mins = ext_hyper_mins
    )

    # Use episode count defined by guidelines as "event count"
    res[[paste0("hyper_event_",  th)]]  <- ep_h$n
    res[[paste0("hyper_mins_",   th)]]  <- ep_h$total_mins
    res[[paste0("hyper_avg_dur_", th)]] <- ep_h$mean_dur
    res[[paste0("hyper_avg_glu_", th)]] <- ep_h$mean_glu

    # Extended episode metrics
    res[[paste0("hyper_ext_event_",    th)]] <- ep_h$n_ext
    res[[paste0("hyper_ext_mins_",     th)]] <- ep_h$total_mins_ext
    res[[paste0("hyper_ext_avg_dur_",  th)]] <- ep_h$mean_dur_ext
    res[[paste0("hyper_ext_avg_glu_",  th)]] <- ep_h$mean_glu_ext
  }

  # Hypoglycemia: TBR + AUC + episode
  for (th in under_list) {

    res[[paste0("tbr_prop_",  th)]] <-
      .percent_time(ifelse(sg < th, 1, 0), interval)
    res[[paste0("auc_below_", th)]] <-
      .auc_excess(sg, th, interval, type = "under")

    ep_l <- .episode_stats_with_end(
      values        = tb$sensorglucose,
      threshold     = th,
      interval_secs = interval,
      direction     = "under",
      min_len_mins  = below_len_mins,
      end_len_mins  = hypo_end_mins,
      extended_mins = ext_hypo_mins
    )

    res[[paste0("hypo_event_",   th)]]  <- ep_l$n
    res[[paste0("hypo_mins_",    th)]]  <- ep_l$total_mins
    res[[paste0("hypo_avg_dur_", th)]]  <- ep_l$mean_dur
    res[[paste0("hypo_avg_glu_", th)]]  <- ep_l$mean_glu

    res[[paste0("hypo_ext_event_",    th)]] <- ep_l$n_ext
    res[[paste0("hypo_ext_mins_",     th)]] <- ep_l$total_mins_ext
    res[[paste0("hypo_ext_avg_dur_",  th)]] <- ep_l$mean_dur_ext
    res[[paste0("hypo_ext_avg_glu_",  th)]] <- ep_l$mean_glu_ext
  }

  # TIR: 63–140, 70–140, 70–180
  for (rng in list(c(63,140), c(70,140), c(70,180))) {
    lo <- rng[1]; hi <- rng[2]
    inrange <- (sg >= lo & sg <= hi)
    res[[paste0("tir_prop_", lo, "_", hi)]] <-
      .percent_time(inrange, interval)
    res[[paste0("auc_in_",   lo, "_", hi)]] <-
      .auc_inrange(sg, lo, hi, interval, mode = "clipped")
  }

  res
}

## ======================== 7) Day/Night Split Blocks =========================
## Day/Night TIR/TBR/TAR, Summary, Variability, Indices

.block_day_night_defaults <- function(tb, interval, daystart, dayend,
                                      hypo_dur_mins, hyper_dur_mins,
                                      gri_weights = c(3.0, 2.4, 0.8, 1.6),
                                      gri_targets = c(54, 70, 180, 250)) {

  res <- list()

  if ("wake" %in% colnames(tb)) {
    daytime_idx   <- which(tb$wake == 1)
    nighttime_idx <- which(tb$wake == 0)
  } else {
    hr <- as.integer(format(tb$timestamp, "%H"))

    hours_in_day <-
      if (daystart == dayend) {
        0:23
      } else if (daystart < dayend) {
        daystart:dayend
      } else {
        c(daystart:23, 0:dayend)
      }

    daytime_idx   <- which(hr %in% hours_in_day)
    nighttime_idx <- setdiff(seq_along(hr), daytime_idx)
  }

  # ---- day/night reading counts + effective days (based on daystart/dayend) ----
  res$total_glu_readings_day   <- .n_valid(tb$sensorglucose[daytime_idx])
  res$total_glu_readings_night <- .n_valid(tb$sensorglucose[nighttime_idx])

  # ---- additive wear-days (STRICT: day + night == all) ----
  res$num_days_day   <- .wear_days_from_n(res$total_glu_readings_day,   interval)
  res$num_days_night <- .wear_days_from_n(res$total_glu_readings_night, interval)

  # additive wear-days (day + night == all), used for daily_avg denominator
  res$num_days_wear  <- res$num_days_day + res$num_days_night

  .one <- function(vals, ts, tag, hypo_dur_mins, hyper_dur_mins) {
    if (!length(vals) || !length(ts)) return(list())

    vals <- as.numeric(vals)

    ok <- is.finite(vals) & !is.na(ts)
    vals <- vals[ok]
    ts   <- ts[ok]

    if (!length(vals)) return(list())

    out <- list()

    ## TIR (63–140, 70–140, 70–180)
    for (rng in list(c(63,140), c(70,140), c(70,180))) {
      lo <- rng[1]; hi <- rng[2]
      flag <- ifelse(vals >= lo & vals <= hi, 1, 0)
      out[[paste0("tir_prop_", lo, "_", hi, "_", tag)]] <-
        .percent_time(flag, interval)
      out[[paste0("auc_in_", lo, "_", hi, "_", tag)]] <-
        .auc_inrange(vals, lo, hi, interval, mode = "clipped")
    }

    ## Hypoglycemia (TBR: Threshold Below Range)
    tbr_list <- c(40, 54, 63, 70)
    for (thr in tbr_list) {
      flag <- ifelse(vals < thr, 1, 0)
      out[[paste0("tbr_prop_", thr, "_", tag)]] <-
        .percent_time(flag, interval)
      out[[paste0("auc_below_",     thr, "_", tag)]] <-
        .auc_excess(vals, thr, interval, type = "under")

      ### Day/Night episode: Using episode_stats_with_end, end_len=0 (no bridging)
      ep_l <- .episode_stats_with_end(
        values        = vals,
        threshold     = thr,
        interval_secs = interval,
        direction     = "under",
        min_len_mins  = hypo_dur_mins,
        end_len_mins  = 0,
        extended_mins = 1e9
      )

      out[[paste0("hypo_event_",    thr, "_", tag)]] <- ep_l$n
      out[[paste0("hypo_mins_",     thr, "_", tag)]] <- ep_l$total_mins
      out[[paste0("hypo_avg_dur_",  thr, "_", tag)]] <- ep_l$mean_dur
    }

    ## Hyperglycemia (TAR: Threshold Above Range)
    tar_list <- c(140, 180, 250, 400)
    for (thr in tar_list) {
      flag <- ifelse(vals > thr, 1, 0)
      flag <- ifelse(vals > thr, 1, 0)
      out[[paste0("tar_prop_", thr, "_", tag)]] <-
        .percent_time(flag, interval)
      out[[paste0("auc_above_",     thr, "_", tag)]] <-
        .auc_excess(vals, thr, interval, type = "over")

      ep_h <- .episode_stats_with_end(
        values        = vals,
        threshold     = thr,
        interval_secs = interval,
        direction     = "over",
        min_len_mins  = hyper_dur_mins,
        end_len_mins  = 0,
        extended_mins = 1e9
      )

      out[[paste0("hyper_event_",   thr, "_", tag)]] <- ep_h$n
      out[[paste0("hyper_mins_",    thr, "_", tag)]] <- ep_h$total_mins
      out[[paste0("hyper_avg_dur_", thr, "_", tag)]] <- ep_h$mean_dur
    }

    ## Summary statistics
    out[[paste0("mean_glu_", tag)]]   <- mean(vals)
    out[[paste0("median_glu_", tag)]] <- stats::median(vals)

    out[[paste0("gmi_", tag)]]    <- 3.31 + (0.02392 * mean(vals, na.rm = TRUE))
    out[[paste0("ehba1c_", tag)]] <- (46.7 + mean(vals, na.rm = TRUE)) / 28.7

    qs <- stats::quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
    out[[paste0("p25_glu_", tag)]]    <- qs[1]
    out[[paste0("p75_glu_", tag)]]    <- qs[2]
    out[[paste0("iqr_glu_",  tag)]]   <- qs[2] - qs[1]

    cqv <- if ((qs[2] + qs[1]) == 0 || is.na(qs[2]) || is.na(qs[1])) NA_real_
    else ((qs[2] - qs[1]) / (qs[2] + qs[1])) * 100
    out[[paste0("cqv_glu_", tag)]] <- cqv

    out[[paste0("sd_glu_",   tag)]]   <- stats::sd(vals)
    mu <- mean(vals)
    out[[paste0("cv_glu_", tag)]] <-
      if (is.na(mu) || mu == 0) NA_real_ else (stats::sd(vals) / mu) * 100
    out[[paste0("min_glu_",  tag)]]   <- min(vals)
    out[[paste0("max_glu_",  tag)]]   <- max(vals)
    out[[paste0("range_glu_",tag)]]   <- max(vals) - min(vals)

    ## Interday variability (based on daily mean)
    if (!all(is.na(vals))) {
      daily <- stats::aggregate(
        vals,
        by = list(day = as.Date(ts)),
        FUN = function(x) mean(x, na.rm = TRUE)
      )
      daily_means <- daily$x
      if (length(stats::na.omit(daily_means)) >= 2) {
        inter_sd <- stats::sd(daily_means, na.rm = TRUE)
        inter_mu <- base::mean(daily_means, na.rm = TRUE)
        inter_cv <- if (is.na(inter_mu) || inter_mu == 0) NA_real_
        else inter_sd / inter_mu * 100
      } else {
        inter_sd <- NA_real_
        inter_cv <- NA_real_
      }
    } else {
      inter_sd <- NA_real_
      inter_cv <- NA_real_
    }

    out[[paste0("interdaysd_", tag)]] <- inter_sd
    out[[paste0("interdaycv_", tag)]] <- inter_cv

    ## Intraday variability (Summary of daily sd / cv)
    if (!all(is.na(vals))) {
      by_day <- split(vals, as.Date(ts))

      daily_sd <- vapply(by_day, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2) return(NA_real_)
        stats::sd(x, na.rm = TRUE)
      }, numeric(1))

      daily_cv <- vapply(by_day, function(x) {
        x <- x[is.finite(x)]
        if (length(x) < 2) return(NA_real_)
        mu2 <- base::mean(x, na.rm = TRUE)
        if (is.na(mu2) || mu2 == 0) return(NA_real_)
        stats::sd(x, na.rm = TRUE) / mu2 * 100
      }, numeric(1))

      sd_ok <- daily_sd[is.finite(daily_sd)]
      cv_ok <- daily_cv[is.finite(daily_cv)]

      if (length(sd_ok)) {
        out[[paste0("intradaysd_mean_",   tag)]] <- base::mean(sd_ok, na.rm = TRUE)
        out[[paste0("intradaysd_median_", tag)]] <- stats::median(sd_ok, na.rm = TRUE)
        out[[paste0("intradaysd_sd_",     tag)]] <- if (length(sd_ok) >= 2) stats::sd(sd_ok, na.rm = TRUE) else NA_real_
      } else {
        out[[paste0("intradaysd_mean_",   tag)]] <- NA_real_
        out[[paste0("intradaysd_median_", tag)]] <- NA_real_
        out[[paste0("intradaysd_sd_",     tag)]] <- NA_real_
      }

      if (length(cv_ok)) {
        out[[paste0("intradaycv_mean_",   tag)]] <- base::mean(cv_ok, na.rm = TRUE)
        out[[paste0("intradaycv_median_", tag)]] <- stats::median(cv_ok, na.rm = TRUE)
        out[[paste0("intradaycv_sd_",     tag)]] <- if (length(cv_ok) >= 2) stats::sd(cv_ok, na.rm = TRUE) else NA_real_
      } else {
        out[[paste0("intradaycv_mean_",   tag)]] <- NA_real_
        out[[paste0("intradaycv_median_", tag)]] <- NA_real_
        out[[paste0("intradaycv_sd_",     tag)]] <- NA_real_
      }
    } else {
      out[[paste0("intradaysd_mean_",   tag)]] <- NA_real_
      out[[paste0("intradaysd_median_", tag)]] <- NA_real_
      out[[paste0("intradaysd_sd_",     tag)]] <- NA_real_
      out[[paste0("intradaycv_mean_",   tag)]] <- NA_real_
      out[[paste0("intradaycv_median_", tag)]] <- NA_real_
      out[[paste0("intradaycv_sd_",     tag)]] <- NA_real_
    }

    ## GRI (Day/Night)
    cuts <- as.numeric(gri_targets)
    if (length(cuts) != 4L)
      stop("gri_targets must be a numeric vector of length 4: c(vlow, low, high, vhigh).")
    if (any(!is.finite(cuts)))
      stop("gri_targets must be finite numeric values.")

    vlow_cut  <- cuts[1]
    low_cut   <- cuts[2]
    high_cut  <- cuts[3]
    vhigh_cut <- cuts[4]

    out[[paste0("gri_", tag)]] <- .gri_from_sg(
      vals, interval,
      vlow_cut  = vlow_cut,
      low_cut   = low_cut,
      high_cut  = high_cut,
      vhigh_cut = vhigh_cut,
      w_vlow  = gri_weights[1],
      w_low   = gri_weights[2],
      w_high  = gri_weights[3],
      w_vhigh = gri_weights[4]
    )

    ## LBGI / HBGI (Day/Night)
    tmp <- .lbgi_hbgi_from_sg(vals)
    out[[paste0("lbgi_", tag)]] <- tmp$lbgi
    out[[paste0("hbgi_", tag)]] <- tmp$hbgi

    ## GRADE (Day/Night)
    g <- .block_grade(vals)
    out[[paste0("grade_score_",      tag)]] <- g$grade_score
    out[[paste0("grade_hypo_prop_",  tag)]] <- g$grade_hypo_prop
    out[[paste0("grade_eugly_prop_", tag)]] <- g$grade_eugly_prop
    out[[paste0("grade_hyper_prop_", tag)]] <- g$grade_hyper_prop

    ## M-value (Day/Night)
    out[[paste0("m_value_", tag)]] <- .m_value_from_sg(vals, r = 90)

    ## MAD (Day/Night)
    out[[paste0("mad_", tag)]] <- .mad_from_sg(vals, constant = 1.4826)

    ## Hyperindex / Hypoindex / IGC (Day/Night)
    hi <- .hyper_index_from_sg(vals)
    ho <- .hypo_index_from_sg(vals)

    out[[paste0("hyper_index_", tag)]] <- hi
    out[[paste0("hypo_index_",  tag)]] <- ho

    if (all(is.na(c(hi, ho)))) {
      out[[paste0("igc_", tag)]] <- NA_real_
    } else {
      out[[paste0("igc_", tag)]] <- sum(hi, ho, na.rm = TRUE)
    }

    out
  }

  day_out   <- .one(tb$sensorglucose[daytime_idx],   tb$timestamp[daytime_idx],   "day",   hypo_dur_mins, hyper_dur_mins)
  night_out <- .one(tb$sensorglucose[nighttime_idx], tb$timestamp[nighttime_idx], "night", hypo_dur_mins, hyper_dur_mins)
  c(res, day_out, night_out)
}

## ====================== 8) Variability Metrics Blocks =======================
## MAGE, MODD, CONGA, GVP, ADRR, MAG, ROC, etc.

# MODD (Mean of Daily Differences)
# https:// doi.org/10.2337/diacare.3.1.58
.compute_modd <- function(tb, lag = 1, interval = NULL) {
  ts <- tb$timestamp
  sg <- tb$sensorglucose

  if (!length(ts) || length(ts) != length(sg)) return(NA_real_)

  if (is.null(interval)) {
    interval <- suppressWarnings(
      stats::median(diff(as.numeric(ts)), na.rm = TRUE)
    )
  }
  if (!is.finite(interval) || interval <= 0) return(NA_real_)

  day_secs    <- 24 * 3600
  pts_per_day <- round(day_secs / interval)
  if (!is.finite(pts_per_day) || pts_per_day <= 0) return(NA_real_)

  tb$day <- as.Date(ts)

  tz <- attr(ts, "tzone")
  if (is.null(tz) || is.na(tz)) tz <- "UTC"

  tb$pos <- as.numeric(difftime(
    ts,
    as.POSIXct(tb$day, tz = tz),
    units = "secs"
  )) / interval + 1
  tb$pos <- round(tb$pos)

  modd_vals <- numeric(0)

  for (p in 1:pts_per_day) {
    grp <- tb[tb$pos == p, ]
    if (nrow(grp) < (lag + 1)) next

    grp <- grp[order(grp$day), ]
    diffs_p <- abs(diff(grp$sensorglucose, lag = lag))
    if (length(diffs_p)) {
      modd_vals <- c(modd_vals, diffs_p)
    }
  }

  if (length(modd_vals)) {
    mean(modd_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
}

# ROC (Rate Of Change) summary (mg/dL per min)
# https://doi.org/10.1089/dia.2008.0138
.roc_summary_from_tb <- function(tb,
                                 interval_secs,
                                 timelag = 15) {

  out <- list(
    roc_mean      = NA_real_,
    roc_sd        = NA_real_,
    roc_mean_abs  = NA_real_,
    roc_p95_abs   = NA_real_,
    roc_max       = NA_real_,
    roc_min       = NA_real_
  )

  sg <- as.numeric(tb$sensorglucose)
  sg <- sg[is.finite(sg)]
  n  <- length(sg)

  if (!n || !is.finite(interval_secs) || interval_secs <= 0) {
    return(out)
  }

  dt_min <- interval_secs / 60

  lag_pts <- round(timelag / dt_min)
  if (!is.finite(lag_pts) || lag_pts < 1L) lag_pts <- 1L
  if (lag_pts >= n) return(out)

  roc_vec <- c(
    rep(NA_real_, lag_pts),
    (sg[(lag_pts + 1):n] - sg[1:(n - lag_pts)]) / timelag
  )

  roc_valid <- roc_vec[is.finite(roc_vec)]
  if (!length(roc_valid)) return(out)

  abs_roc <- abs(roc_valid)

  out$roc_mean      <- mean(roc_valid, na.rm = TRUE)
  out$roc_sd        <- stats::sd(roc_valid, na.rm = TRUE)
  out$roc_mean_abs  <- mean(abs_roc, na.rm = TRUE)
  out$roc_p95_abs   <- as.numeric(stats::quantile(abs_roc, 0.95, na.rm = TRUE))
  out$roc_max       <- max(roc_valid, na.rm = TRUE)
  out$roc_min       <- min(roc_valid, na.rm = TRUE)

  out
}

# Glucose variability (main block)
.block_variability <- function(tb,
                               interval,
                               mage_sd,
                               conga_hrs,
                               modd_lag,
                               ac_lag = 30,
                               roc_lag = 15) {
  res <- list()

  tb$smoothed <- as.numeric(zoo::rollapply(
    zoo::zoo(tb$sensorglucose), 9,
    function(x) { c(1,2,4,8,16,8,4,2,1) %*% (x/46) },
    fill = NA))

  tb$smoothed[1:4] <- mean(stats::na.omit(tb$sensorglucose[1:4]))
  tb$smoothed[(length(tb$smoothed) - 3):length(tb$smoothed)] <-
    mean(tb$sensorglucose[(length(tb$sensorglucose) - 3):length(tb$sensorglucose)])

  mu_all <- base::mean(tb$sensorglucose, na.rm = TRUE)
  sd_all <- stats::sd(tb$sensorglucose,  na.rm = TRUE)

  # Interday SD / CV
  # doi: 10.1089/dia.2009.0015
  # https://doi.org/ 10.1016/j.amjms.2018.09.010
  if (!all(is.na(tb$sensorglucose))) {
    daily <- stats::aggregate(
      tb$sensorglucose,
      by = list(day = as.Date(tb$timestamp)),
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    daily_means <- daily$x

    if (length(stats::na.omit(daily_means)) >= 2) {
      inter_sd <- stats::sd(daily_means, na.rm = TRUE)
      inter_mu <- base::mean(daily_means, na.rm = TRUE)
      inter_cv <- if (is.na(inter_mu) || inter_mu == 0) NA_real_
      else inter_sd / inter_mu * 100
    } else {
      inter_sd <- NA_real_
      inter_cv <- NA_real_
    }
  } else {
    inter_sd <- NA_real_
    inter_cv <- NA_real_
  }

  res$interdaysd <- inter_sd
  res$interdaycv <- inter_cv

  # Intraday SD / CV Summary
  # doi: 10.1089/dia.2009.0015
  # https://doi.org/ 10.1016/j.amjms.2018.09.010
  # 73-LB: expanding the definition of Intraday glucose variability
  if (!all(is.na(tb$sensorglucose))) {
    by_day <- split(tb$sensorglucose, as.Date(tb$timestamp))

    daily_sd <- vapply(by_day, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NA_real_)
      stats::sd(x, na.rm = TRUE)
    }, numeric(1))

    daily_cv <- vapply(by_day, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NA_real_)
      mu <- base::mean(x, na.rm = TRUE)
      if (is.na(mu) || mu == 0) return(NA_real_)
      stats::sd(x, na.rm = TRUE) / mu * 100
    }, numeric(1))

    sd_ok <- daily_sd[is.finite(daily_sd)]
    cv_ok <- daily_cv[is.finite(daily_cv)]

    if (length(sd_ok)) {
      res$intradaysd_mean    <- base::mean(sd_ok, na.rm = TRUE)
      res$intradaysd_median  <- stats::median(sd_ok, na.rm = TRUE)
      res$intradaysd_sd      <- if (length(sd_ok) >= 2) stats::sd(sd_ok, na.rm = TRUE) else NA_real_
    } else {
      res$intradaysd_mean    <- NA_real_
      res$intradaysd_median  <- NA_real_
      res$intradaysd_sd      <- NA_real_
    }

    if (length(cv_ok)) {
      res$intradaycv_mean    <- base::mean(cv_ok, na.rm = TRUE)
      res$intradaycv_median  <- stats::median(cv_ok, na.rm = TRUE)
      res$intradaycv_sd      <- if (length(cv_ok) >= 2) stats::sd(cv_ok, na.rm = TRUE) else NA_real_
    } else {
      res$intradaycv_mean    <- NA_real_
      res$intradaycv_median  <- NA_real_
      res$intradaycv_sd      <- NA_real_
    }
  } else {
    res$intradaysd_mean   <- NA_real_
    res$intradaysd_median <- NA_real_
    res$intradaysd_sd     <- NA_real_
    res$intradaycv_mean   <- NA_real_
    res$intradaycv_median <- NA_real_
    res$intradaycv_sd     <- NA_real_
  }

  # MAGE
  # https:// doi.org/10.2337/diacare.3.1.5
  tp    <- pastecs::turnpoints(tb$smoothed)
  peaks <- tp$pos[tp$peaks]
  pits  <- tp$pos[tp$pits]

  if (isTRUE(tp[["firstispeak"]]) && length(peaks) != length(pits)) {
    peaks <- peaks[2:length(peaks)]
  } else if (!isTRUE(tp[["firstispeak"]]) && length(peaks) != length(pits)) {
    pits <- pits[1:(length(pits) - 1)]
  }

  diffs <- tb$sensorglucose[peaks] - tb$sensorglucose[pits]
  res$r_mage <- mean(stats::na.omit(diffs[diffs > (mage_sd * sd_all)]))

  # MGE / MGN + PTIR / PTOR
  # doi: 10.2337/db12-1396
  # 10.1109/OJEMB.2021.3105816
  # doi: 10.1177/1932296819826111
  if (is.na(mu_all) || is.na(sd_all)) {
    res$mge       <- NA_real_
    res$mgn       <- NA_real_
    res$ptir_stat <- NA_real_
    res$ptor_stat <- NA_real_
  } else {
    up   <- mu_all + mage_sd * sd_all
    dw   <- mu_all - mage_sd * sd_all
    vals <- tb$sensorglucose

    outside <- vals[vals >= up | vals <= dw]
    inside  <- vals[vals >  dw & vals <  up]

    res$mge <- if (length(outside)) base::mean(outside, na.rm = TRUE) else NA_real_
    res$mgn <- if (length(inside))  base::mean(inside,  na.rm = TRUE) else NA_real_

    inside_flag  <- vals >= dw & vals <= up
    outside_flag <- !inside_flag
    total_n <- sum(!is.na(vals))

    if (total_n > 0) {
      res$ptir_stat <- 100 * sum(inside_flag,  na.rm = TRUE) / total_n
      res$ptor_stat <- 100 * sum(outside_flag, na.rm = TRUE) / total_n
    } else {
      res$ptir_stat <- NA_real_
      res$ptor_stat <- NA_real_
    }
  }

  # J-index
  # https://doi.org/10.1055/s-2007979906
  mu  <- base::mean(tb$sensorglucose, na.rm = TRUE)
  sdv <- stats::sd(tb$sensorglucose, na.rm = TRUE)
  if (is.na(mu) || is.na(sdv)) {
    res$j_index <- NA_real_
  } else {
    res$j_index <- 0.001 * (mu + sdv)^2
  }

  # CONGA
  # https://doi.org/10.1089/dia.2005.7.253
  n <- (conga_hrs * 3600)
  conga.times <- tb$timestamp + n
  conga.times <- conga.times[!is.na(conga.times)]
  conga.times <- conga.times[order(conga.times)]
  conga.times <- conga.times[conga.times %in% tb$timestamp]
  begin.times <- conga.times - n
  congas <- suppressWarnings({
    tb$sensorglucose[tb$timestamp %in% conga.times] -
      tb$sensorglucose[tb$timestamp %in% begin.times]
  })
  res[[paste0("conga_", conga_hrs, "h")]] <- stats::sd(congas, na.rm = TRUE)

  # MODD
  # https:// doi.org/10.2337/diacare.3.1.58
  res$modd <- .compute_modd(tb, lag = modd_lag, interval = interval)

  # LBGI / HBGI (Global)
  # https://doi.org/10.2337/dc06-1085
  tmp <- .lbgi_hbgi_from_sg(tb$sensorglucose)
  res$lbgi <- tmp$lbgi
  res$hbgi <- tmp$hbgi

  # M-value
  # https://doi.org/10.1111/j.0954-6820.1965.tb01810.x
  res$m_value <- .m_value_from_sg(tb$sensorglucose, r = 90)

  # MAD (Median Absolute Deviation)
  # https://doi.org/10.1371/ journal.pone.0248560
  res$mad <- .mad_from_sg(tb$sensorglucose, constant = 1.4826)

  # ADRR
  # https://doi.org/10.2337/dc06-1085
  tb$gluctransform <- (log(tb$sensorglucose)^1.084) - 5.381
  tb$rBG <- 22.77 * (tb$gluctransform^2)
  tb$rl  <- ifelse(tb$gluctransform <= 0, tb$rBG, 0)
  tb$rh  <- ifelse(tb$gluctransform >  0, tb$rBG, 0)

  tb$Date <- as.Date(tb$timestamp)
  dates   <- unique(tb$Date)
  adrr <- sapply(dates, function(d) {
    rl_d <- tb$rl[tb$Date == d]
    rh_d <- tb$rh[tb$Date == d]
    if (all(is.na(rl_d)) && all(is.na(rh_d))) return(NA_real_)
    max(rl_d, na.rm = TRUE) + max(rh_d, na.rm = TRUE)
  })
  res$adrr <- mean(adrr, na.rm = TRUE)

  # MAG
  # https:// doi.org/10.1097/CCM.0b013e3181cc4be9
  delta_bg <- sum(abs(diff(tb$sensorglucose)), na.rm = TRUE)
  delta_t  <- (interval * sum(!is.na(tb$sensorglucose))) / 3600
  res$mag  <- delta_bg / delta_t

  # GVP: Glucose Variability Percentage
  # doi: 10.1089/dia.2017.0187
  res$gvp <- .gvp_from_sg(tb$sensorglucose, interval)

  # Hyperglycemia Index / Hypoglycemia Index / IGC
  # doi: 10.1089/dia.2008.0132
  res$hyper_index <- .hyper_index_from_sg(tb$sensorglucose)
  res$hypo_index  <- .hypo_index_from_sg(tb$sensorglucose)

  if (all(is.na(c(res$hyper_index, res$hypo_index)))) {
    res$igc <- NA_real_
  } else {
    res$igc <- sum(res$hyper_index, res$hypo_index, na.rm = TRUE)
  }

  # PGS (Personal Glycemic State)
  # doi: 10.1089/dia.2017.0080
  res$pgs <- tryCatch({
    sg <- as.numeric(tb$sensorglucose)
    sg <- sg[is.finite(sg)]
    if (!length(sg)) return(NA_real_)

    GVP <- .gvp_from_sg(tb$sensorglucose, interval)
    MG  <- mean(tb$sensorglucose, na.rm = TRUE)
    PTIR <- .percent_time(
      ifelse(tb$sensorglucose >= 70 & tb$sensorglucose <= 180, 1, 0),
      interval
    )

    ndays_eff <- .effective_days(tb, interval)

    if (!is.finite(ndays_eff) || ndays_eff <= 0) {
      N54 <- 0
      N70 <- 0
    } else {
      ep54 <- .episode_stats_with_end(
        values        = tb$sensorglucose,
        threshold     = 54,
        interval_secs = interval,
        direction     = "under",
        min_len_mins  = 20,
        end_len_mins  = 30,
        extended_mins = 1e9
      )
      N54 <- (ep54$n / ndays_eff) * 7

      ep70 <- .episode_stats_with_end(
        values        = tb$sensorglucose,
        threshold     = 70,
        interval_secs = interval,
        direction     = "under",
        min_len_mins  = 20,
        end_len_mins  = 30,
        extended_mins = 1e9
      )
      N70_only <- max(ep70$n - ep54$n, 0)
      N70 <- (N70_only / ndays_eff) * 7
    }

    f_gvp <- 1 + (9 / (1 + exp(-0.049  * (GVP - 65.47))))
    f_ptir <- 1 + (9 / (1 + exp( 0.0833 * (PTIR - 55.04))))
    f_mg <- 1 + 9 * (
      1 / (1 + exp( 0.1139  * (MG - 72.08))) +
        1 / (1 + exp(-0.09195 * (MG - 157.57)))
    )

    f_h54 <- 0.5 + 4.5 * (1 - exp(-0.81093 * N54))
    f_h70 <- ifelse(N70 <= 7.65,
                    0.5714 * N70 + 0.625,
                    5)

    as.numeric(f_gvp + f_ptir + f_mg + f_h54 + f_h70)
  }, error = function(e) NA_real_)

  # AC_Mean / AC_Var from ACF
  # https://doi.org/10.1038/s43856-025-00819-5
  ac_info <- .ac_summary_from_sg(tb$sensorglucose, ac_lag = ac_lag)
  res$ac_mean <- ac_info$ac_mean
  res$ac_var  <- ac_info$ac_var

  # ROC summary
  res <- c(res,
           .roc_summary_from_tb(tb,
                                interval_secs = interval,
                                timelag       = roc_lag))
  res
}

## ======================= 9) Scanned Threshold Blocks ========================
## tbr_seq / tar_seq Scanning (Full Day + Day/Night)
.block_scanned_intervals_optional <- function(tb, interval, tbr_seq, tar_seq,
                                              above_len_mins, below_len_mins,
                                              calc_day_night = NULL,
                                              daystart = 6, dayend = 0,
                                              hypo_dur_mins = below_len_mins,
                                              hyper_dur_mins = above_len_mins) {
  sg <- tb$sensorglucose
  res <- list()

  .scan_one <- function(vals, tag = NULL) {
    out <- list()

    ## Hypoglycemia scanning (TBR: Threshold Below Range)
    if (!is.null(tbr_seq)) {
      lo_seq <- ._expand_seq(tbr_seq)
      for (lo in lo_seq) {
        flag <- ifelse(vals < lo, 1, 0)
        key <- function(k) if (is.null(tag)) paste0(k, "_", lo) else paste0(k, "_", lo, "_", tag)

        out[[key("tbr_prop")]]   <- .percent_time(flag, interval)
        out[[key("auc_below")]]  <- .auc_excess(vals, lo, interval, type = "under")

        ep_l <- .episode_stats_with_end(
          values        = vals,
          threshold     = lo,
          interval_secs = interval,
          direction     = "under",
          min_len_mins  = hypo_dur_mins,
          end_len_mins  = 0,
          extended_mins = 0
        )

        out[[key("hypo_event")]]   <- ep_l$n
        out[[key("hypo_mins")]]    <- ep_l$total_mins
        out[[key("hypo_avg_dur")]] <- ep_l$mean_dur
      }
    }

    ## Hyperglycemia scanning (TAR: Threshold Above Range)
    if (!is.null(tar_seq)) {
      hi_seq <- ._expand_seq(tar_seq)
      for (hi in hi_seq) {
        flag <- ifelse(vals > hi, 1, 0)
        key <- function(k) if (is.null(tag)) paste0(k, "_", hi) else paste0(k, "_", hi, "_", tag)

        out[[key("tar_prop")]]   <- .percent_time(flag, interval)
        out[[key("auc_above")]]  <- .auc_excess(vals, hi, interval, type = "over")

        ep_h <- .episode_stats_with_end(
          values        = vals,
          threshold     = hi,
          interval_secs = interval,
          direction     = "over",
          min_len_mins  = hyper_dur_mins,
          end_len_mins  = 0,
          extended_mins = 0
        )

        out[[key("hyper_event")]]   <- ep_h$n
        out[[key("hyper_mins")]]    <- ep_h$total_mins
        out[[key("hyper_avg_dur")]] <- ep_h$mean_dur
      }
    }

    ## TIR (Target Range): Only when both tbr_seq and tar_seq are provided, scanning will occur
    if (!is.null(tbr_seq) && !is.null(tar_seq)) {
      lo_seq <- ._expand_seq(tbr_seq)
      hi_seq <- ._expand_seq(tar_seq)
      for (lo in lo_seq) for (hi in hi_seq) {
        if (lo > hi) next
        inrange <- (vals >= lo & vals <= hi)
        key <- function(k) {
          if (is.null(tag)) paste0(k, "_", lo, "_", hi)
          else paste0(k, "_", lo, "_", hi, "_", tag)
        }
        out[[key("tir_prop")]] <- .percent_time(inrange, interval)
      }
    }
    out
  }

  ## Full day scan
  res <- c(res, .scan_one(sg, tag = NULL))

  ## If day/night classification is provided, perform day/night scanning
  if (!is.null(calc_day_night)) {
    if ("wake" %in% names(tb)) {
      day_idx   <- which(tb$wake == 1)
      night_idx <- which(tb$wake == 0)
    } else {
      hr <- as.integer(format(tb$timestamp, "%H"))
      hours_in_day <- if (daystart == dayend) 0:23
      else if (daystart < dayend) daystart:dayend
      else c(daystart:23, 0:dayend)
      day_idx   <- which(hr %in% hours_in_day)
      night_idx <- setdiff(seq_along(hr), day_idx)
    }
    if (length(day_idx))   res <- c(res, .scan_one(sg[day_idx],   tag = "day"))
    if (length(night_idx)) res <- c(res, .scan_one(sg[night_idx], tag = "night"))
  }

  res
}

## =========================== 10) Post-processing ============================

# Create binary flags for hypo/hyper events: 1 if event count > 0 else 0
.add_event_binary <- function(res_list) {
  nn <- names(res_list)
  if (!length(nn)) return(res_list)

  for (key in nn) {
    if (!nzchar(key)) next
    val <- suppressWarnings(as.numeric(res_list[[key]]))

    if (!is.finite(val)) next

    # regular hypo / hyper
    if (grepl("^hypo_event_", key)) {
      res_list[[sub("^hypo_event_", "hypo_binary_", key)]] <- as.integer(val > 0)
    }
    if (grepl("^hyper_event_", key)) {
      res_list[[sub("^hyper_event_", "hyper_binary_", key)]] <- as.integer(val > 0)
    }

    # extended hypo / hyper
    if (grepl("^hypo_ext_event_", key)) {
      res_list[[sub("^hypo_ext_event_", "hypo_ext_binary_", key)]] <- as.integer(val > 0)
    }
    if (grepl("^hyper_ext_event_", key)) {
      res_list[[sub("^hyper_ext_event_", "hyper_ext_binary_", key)]] <- as.integer(val > 0)
    }
  }

  res_list
}

## Daily Averaging, Post-processing

# Divides event/minutes/AUC variables by num_days if needed
.dailyize_metrics <- function(res_list, num_days, daily_avg) {
  nd <- suppressWarnings(as.numeric(num_days))
  if (is.na(nd) || nd <= 0) nd <- 1
  if (!isTRUE(daily_avg)) return(res_list)

  pats <- c(
    "^hypo_event_",      "^hyper_event_",
    "^hypo_mins_",       "^hyper_mins_",
    "^hypo_ext_event_",  "^hyper_ext_event_",
    "^hypo_ext_mins_",   "^hyper_ext_mins_",
    "^auc_above_",  "^auc_below_", "^auc_in_"
  )

  nn <- names(res_list)
  for (i in seq_along(nn)) {
    key <- nn[i]
    if (!nzchar(key)) next
    if (any(grepl(paste(pats, collapse="|"), key))) {
      val <- suppressWarnings(as.numeric(res_list[[i]]))
      if (is.finite(val)) res_list[[i]] <- val / nd
    }
  }
  res_list
}
