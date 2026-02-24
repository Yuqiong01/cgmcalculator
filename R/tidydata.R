#' Tidy Continuous Glucose Monitoring Data
#'
#' Read raw CGM exports from multiple devices and write standardized,
#' analysis-ready CSV files with consistent columns and timestamps. The function
#' optionally processes missing glucose segments (short-gap interpolation and
#' long-gap handling) and can trim recordings to complete 24-hour days.
#'
#' @description
#' `tidydata()` streamlines pre-processing of CGM datasets for downstream
#' analyses. It detects device-specific export formats, standardizes column
#' names, parses timestamps, performs optional gap processing, and writes cleaned
#' CSV files. Unless manually pre-processed (e.g., following `metricsstat()`),
#' input files should remain in their original export formats. For Diasend and
#' Sibionics Excel files, CGM data must be stored in the first worksheet/tab.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#' \item Read raw CGM files from `inputdir` (recursively).
#' \item Parse timestamps and standardize to two core columns:
#' `timestamp` and `sensorglucose`.
#' \item Optionally skip the initial `skiphours` hours of data.
#' \item Optionally process gaps when `gapprocess = TRUE`:
#'   \itemize{
#'   \item Short gaps (<= `maxgapmin` minutes) are filled using `shortgapsm`
#'         (`"linear"` or `"spline"`).
#'   \item Long gaps (> `maxgapmin` minutes) are handled using `longgapsm`
#'         (`"spline"`, `"dailymean"`, or `"remove"`).
#'   }
#' \item Detect glucose units from the observed range and convert to `unit`
#' (default `"mg/dL"`). Internally, downstream analyses in this package assume
#' mg/dL.
#' \item Optionally trim to complete 24-hour periods when `trimdays = TRUE`.
#' \item Write standardized CSV files to `outputdir`.
#' }
#'
#' Files are skipped and reported if they have too few non-missing glucose values
#' (< 10), if trimming results in an empty dataset, if `trimdays = TRUE` but the
#' total duration is < 24 hours, if the sampling interval cannot be determined,
#' or if the recording is shorter than or equal to `skiphours`.
#'
#' @param inputdir Character. Directory containing raw CGM files to tidy.
#' @param outputdir Character. Directory where cleaned CSV files will be written.
#'   Defaults to `tempdir()`.
#' @param skiphours Numeric. Number of hours to skip at the beginning of each
#'   file. Defaults to 0.
#' @param gapprocess Logical. Whether to detect and handle missing glucose
#'   segments. If `TRUE`, short and long gaps are processed according to
#'   `shortgapsm` and `longgapsm`. If `FALSE`, missing values are preserved.
#' @param maxgapmin Numeric. Maximum gap length (in minutes) treated as a short
#'   gap for interpolation. Defaults to 30.
#' @param shortgapsm Character. Method for filling short gaps (only used when
#'   `gapprocess = TRUE`). Supported values are `"linear"` (default) and
#'   `"spline"`.
#' @param longgapsm Character. Method for handling long gaps (only used when
#'   `gapprocess = TRUE`). Supported values are `"spline"` (default),
#'   `"dailymean"`, and `"remove"`.
#' @param trimdays Logical. Whether to trim data to complete 24-hour periods.
#'   If `TRUE`, partial leading and trailing days are removed. Defaults to
#'   `FALSE`.
#' @param logdisplay Logical. Whether to print progress messages for each file.
#'   Defaults to `TRUE`.
#' @param unit Character. Output glucose unit. Supported values are `"mg/dL"`
#'   (default) and `"mmol/L"`. If conversion is needed, `"mg/dL"` is converted to
#'   `"mmol/L"` by dividing by 18, and `"mmol/L"` is converted to `"mg/dL"` by
#'   multiplying by 18.
#'
#' @return
#' Invisibly returns a named list of skipped/failed file names:
#' \itemize{
#' \item `skipped_few_glucose`: files with fewer than 10 non-missing glucose values.
#' \item `skipped_trim_empty`: files that became empty after trimming.
#' \item `skipped_less_than_1day`: files shorter than 24 hours when `trimdays = TRUE`.
#' \item `skipped_bad_interval`: files with invalid or undeterminable sampling intervals.
#' \item `skipped_skiphours_too_short`: files whose total duration is <= `skiphours`.
#' \item `skipped_errors`: files that failed due to unexpected errors.
#' }
#'
#' @seealso [metricsstats()], [agpanalyze()]
#'
#' @examples
#' \dontrun{
#' tidydata(
#'   inputdir = "rawdata",
#'   outputdir = "tidydata",
#'   skiphours = 0,
#'   gapprocess = TRUE,
#'   maxgapmin = 30,
#'   shortgapsm = "linear",
#'   longgapsm = "spline",
#'   trimdays = FALSE,
#'   logdisplay = TRUE,
#'   unit = "mg/dL"
#' )
#' }
#'
#' @importFrom readr guess_encoding
#' @importFrom readxl read_excel
#' @importFrom XML xmlParse xmlToList
#' @importFrom parsedate parse_date
#' @importFrom pracma Mode
#' @importFrom zoo na.approx na.spline
#'
#' @export
tidydata <- function(inputdir,
                     outputdir = tempdir(),
                     skiphours = 0,
                     gapprocess = TRUE,
                     maxgapmin = 30,
                     shortgapsm = "linear",
                     longgapsm = "spline",
                     trimdays = FALSE,
                     logdisplay = TRUE,
                     unit = "mg/dL") {

  # ----------------------- Step 1: Prepare file list --------------------------
  files <- base::list.files(path = inputdir, full.names = TRUE, recursive = TRUE)
  base::dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

  skipped_few_glucose <- character()
  skipped_trim_empty <- character()
  skipped_less_than_1day <- character()
  skipped_bad_interval <- character()
  skipped_skiphours_too_short <- character()
  error_files <- character()

  # A small helper to standardize selecting columns
  safe_grep <- function(pattern, x, ignore.case = TRUE) {
    col <- grep(pattern, x, ignore.case = ignore.case, value = TRUE)
    if (length(col) == 0) return(NA_character_)
    col[1]
  }

  # Process one file; returns a list(status=..., file=..., msg=...)
  process_one <- function(file_path) {

    if (isTRUE(logdisplay)) {
      base::cat(">> ", base::basename(file_path), "\n")
    }

    ext <- tools::file_ext(file_path)

    # ---------------- Step 2a: Read file based on extension -------------------
    if (ext == "txt") {
      enc <- base::as.character(readr::guess_encoding(file_path)[1, 1])
      table <- utils::read.table(
        file_path,
        sep = "\t",
        skipNul = TRUE,
        header = TRUE,
        stringsAsFactors = FALSE,
        na.strings = "",
        fileEncoding = enc,
        comment.char = ""
      )

    } else if (ext == "csv") {
      table <- utils::read.csv(
        file_path,
        stringsAsFactors = FALSE,
        header = TRUE,
        na.strings = ""
      )
      if (base::ncol(table) <= 2) {
        table <- utils::read.csv(
          file_path,
          sep = ";",
          stringsAsFactors = FALSE,
          header = TRUE,
          na.strings = ""
        )
      }

    } else if (ext %in% c("xls", "xlsx", "xlsm")) {
      table <- suppressMessages(readxl::read_excel(file_path, col_types = "text"))
      table <- base::as.data.frame(table)

    } else if (ext == "xml") {
      doc <- XML::xmlParse(file_path)
      l <- XML::xmlToList(doc)
      id <- l$.attrs[["Id"]]
      l <- l[["GlucoseReadings"]]
      times <- base::lapply(l, function(x) x[["DisplayTime"]])
      times <- base::do.call(rbind, times)
      sensor <- base::lapply(l, function(x) x[["Value"]])
      sensor <- base::do.call(rbind, sensor)
      table <- base::cbind(times, sensor)
      table <- base::as.data.frame(table)
      base::colnames(table) <- c("timestamp", "sensorglucose")
      table$subjectid <- NA
      table$subjectid[1] <- id
      table <- table[, c("subjectid", "timestamp", "sensorglucose")]

    } else if (ext == "ASC") {
      table <- utils::read.delim(file_path)

    } else {
      stop("Unsupported file extension: ", ext)
    }

    # ---------------- Step 2b: Detect CGM type --------------------------------
    cgmtype <- NULL
    ncol_table <- base::ncol(table)
    colnames_table <- base::colnames(table)

    if (ncol_table == 12 && colnames_table[1] == "Index") {
      cgmtype <- "dexcomg6"

    } else if (ncol_table == 3 && colnames_table[1] != "subjectid" &&
               (colnames_table[1] != "Name" || ncol_table == 2)) {
      cgmtype <- "diasend"

    } else if (ncol_table %in% c(18, 19) || ncol_table == 4) {
      if (ncol_table %in% c(18, 19) && table[2, 1] == "Device") {
        cgmtype <- "librepro"
      } else {
        cgmtype <- "libre"
      }

    } else if (ncol_table %in% c(13, 14)) {
      cgmtype <- "dexcomg5"

    } else if (ncol_table >= 47) {
      cgmtype <- "carelink"

    } else if (ncol_table == 3 && colnames_table[2] == "timestamp" &&
               colnames_table[3] == "sensorglucose") {
      cgmtype <- "manual"

    } else if (ncol_table %in% c(17, 22, 34)) {
      cgmtype <- "ipro"

    } else if (ncol_table == 6 && ext == "ASC") {
      cgmtype <- "asc"

    } else if (ext %in% c("xlsx", "xls", "csv") &&
               any(grepl("\u8461\u8404\u7CD6", colnames_table, ignore.case = TRUE)) ||
               any(grepl("\u8840\u7CD6", colnames_table, ignore.case = TRUE)) ||
               any(grepl("\u91C7\u96C6\u65F6\u95F4", colnames_table, ignore.case = TRUE)) ||
               any(grepl("\u65F6\u95F4", colnames_table, ignore.case = TRUE))) {
      cgmtype <- "sibionics"

    } else if (ncol_table == 6) {
      cgmtype <- "tslimg4"

    } else if (ncol_table %in% 39:41) {
      cgmtype <- "tandem"

    } else if (!is.null(table) && base::nrow(table) >= 7 && table[7, 1] == "t:slim X2") {
      cgmtype <- "tslimx2"

    } else if (ncol_table == 3 && base::nrow(table) >= 3 && table[3, 3] == "Serial Number") {
      cgmtype <- "omnipod5"

    } else {
      stop("File formatted incorrectly and cannot be read: ", base::basename(file_path))
    }

    # ---------------- Step 3: Standardize columns -----------------------------
    id <- sub("\\..*", "", base::basename(file_path))

    if (cgmtype == "diasend") {
      table <- table[-c(1:which(table[, 1] == "Time")), ]
      base::colnames(table) <- c("timestamp", "sensorglucose")

    } else if (cgmtype == "carelink") {
      sensor_row <- which(table[, 3] == "Sensor")[1] + 1
      base::colnames(table) <- table[sensor_row, ]
      table <- table[-c(1:sensor_row), ]
      table$timestamp <- paste(table$Date, table$Time)
      table$timestamp <- gsub(".{3}$", "", table$timestamp)
      table <- table[, grep("timestamp|Sensor Glucose", base::colnames(table))]
      table <- table[, sort(base::colnames(table), decreasing = TRUE)]
      base::colnames(table) <- c("timestamp", "sensorglucose")

    } else if (cgmtype %in% c("dexcomg5", "dexcomg6")) {
      ts_col <- safe_grep("timestamp", base::colnames(table))
      gl_col <- safe_grep("glucose", base::colnames(table))
      table$timestamp <- table[, ts_col]
      table$sensorglucose <- table[, gl_col]
      table <- table[, c("timestamp", "sensorglucose")]
      table$timestamp <- sub("T", " ", table$timestamp)

    } else if (cgmtype %in% c("libre", "librepro")) {
      base::colnames(table) <- table[2, ]
      table <- table[-c(1:2), ]
      ts_col <- safe_grep("time|timestamp", base::colnames(table))
      gl_col <- safe_grep("glucose|sensor", base::colnames(table))
      table <- table[, c(ts_col, gl_col)]
      base::colnames(table) <- c("timestamp", "sensorglucose")

    } else if (cgmtype == "manual") {
      table$sensorglucose <- suppressWarnings(as.numeric(as.character(table$sensorglucose)))
      table <- table[1:max(which(!is.na(table$sensorglucose))), ]
      table <- table[, -1]

    } else if (cgmtype == "ipro") {
      base::colnames(table) <- table[11, ]
      table <- table[-c(1:11), ]
      if (!grepl("- | /", table$Timestamp[1])) {
        table$Timestamp <- as.POSIXct(as.numeric(table$Timestamp) * (60 * 60 * 24),
                                      origin = "1899-12-30", tz = "UTC")
      }
      table$Timestamp <- sub("[.]00", "", table$Timestamp)
      table <- table[, c("Timestamp", "Sensor Glucose (mg/dL)")]
      base::colnames(table) <- c("timestamp", "sensorglucose")

    } else if (cgmtype == "asc") {
      table$timestamp <- paste(table$Date, table$Time)
      table$sensorglucose <- table$Value
      table <- table[, c("timestamp", "sensorglucose")]

    } else if (cgmtype == "sibionics") {
      ts_col <- safe_grep("time|timestamp|\u65F6\u95F4|\u91C7\u96C6\u65F6\u95F4", base::colnames(table))
      gl_col <- safe_grep("glucose|sensor|\u8461\u8404\u7CD6|\u8840\u7CD6", base::colnames(table))
      table <- table[, c(ts_col, gl_col)]
      base::colnames(table) <- c("timestamp", "sensorglucose")

    } else if (cgmtype == "tslimg4") {
      row <- which(table[, 1] == "DeviceType")
      base::colnames(table) <- table[row, ]
      table <- table[-c(1:row), ]
      table$timestamp <- table$EventDateTime
      table$sensorglucose <- as.numeric(table$`Readings (CGM / BGM)`)
      table <- table[, c("timestamp", "sensorglucose")]

    } else if (cgmtype %in% c("tandem", "tslimx2", "omnipod5")) {

      if (cgmtype == "tandem") {
        base::colnames(table) <- table[6, ]
        table <- table[-c(1:6), ]
        table <- table[1:max(which(table[, 3] == "EGV")), 4:5]
        base::colnames(table) <- c("timestamp", "sensorglucose")
        table$timestamp <- sub("T", " ", table$timestamp)

      } else if (cgmtype == "tslimx2") {
        row <- which(table[, 1] == "DeviceType")[1]
        base::colnames(table) <- table[row, ]
        table <- table[-c(1:row), ]
        table <- table[, c("EventDateTime", "Readings (mg/dL)")]
        base::colnames(table) <- c("timestamp", "sensorglucose")
        table$timestamp <- sub("T", " ", table$timestamp)

      } else if (cgmtype == "omnipod5") {
        ts_row <- which(table[, 1] == "Timestamp")[1]
        base::colnames(table) <- table[ts_row, ]
        table <- table[-c(1:ts_row), ]
        table <- table[, 1:2]
        base::colnames(table) <- c("timestamp", "sensorglucose")
        table$timestamp <- sub("T", " ", table$timestamp)
      }
    }

    # ------------------------------ Final cleanup -----------------------------
    table$timestamp <- suppressWarnings(parsedate::parse_date(table$timestamp, approx = FALSE))
    table$timestamp <- as.POSIXct(table$timestamp, origin = "1970-01-01", tz = "UTC")

    table$sensorglucose <- suppressWarnings(as.numeric(sub(",", ".", table$sensorglucose)))

    table <- table[!is.na(table$timestamp), ]
    table <- table[!(format(table$timestamp, "%H:%M:%S") == "00:00:00" &
                       is.na(table$sensorglucose)), ]
    table <- table[order(table$timestamp), ]
    table <- table[!duplicated(table), ]

    # -------------------------- Step 4: Unit conversion -----------------------
    gl_values <- table$sensorglucose[!is.na(table$sensorglucose)]
    if (length(gl_values) < 10) {
      msg <- paste0(
        "Skipped file '", base::basename(file_path),
        "': too few glucose values (", length(gl_values),
        ") for reliable unit detection."
      )
      warning(msg)
      return(list(status = "skip_few_glucose", file = base::basename(file_path), msg = msg))
    }

    glucose_sample <- gl_values[1:min(100, length(gl_values))]
    gl_min <- min(glucose_sample, na.rm = TRUE)
    gl_max <- max(glucose_sample, na.rm = TRUE)

    if (gl_max <= 70 && gl_min >= 1) {
      detected_unit <- "mmol/L"
    } else if (gl_max > 70 && gl_max <= 2000) {
      detected_unit <- "mg/dL"
    } else {
      detected_unit <- unit
      warning(paste(
        "File", base::basename(file_path),
        "unit could not be automatically determined. Using provided 'unit' parameter."
      ))
    }

    if (!identical(detected_unit, unit)) {
      if (detected_unit == "mmol/L" && unit == "mg/dL") {
        table$sensorglucose <- table$sensorglucose * 18
        if (isTRUE(logdisplay)) message("Converted mmol/L -> mg/dL for file: ", base::basename(file_path))
      } else if (detected_unit == "mg/dL" && unit == "mmol/L") {
        table$sensorglucose <- table$sensorglucose / 18
        if (isTRUE(logdisplay)) message("Converted mg/dL -> mmol/L for file: ", base::basename(file_path))
      }
    } else {
      if (isTRUE(logdisplay)) message("No conversion needed for file: ", base::basename(file_path), " (unit = ", detected_unit, ")")
    }

    # ---------------- Sampling interval (defensive) ---------------------------
    dt <- base::diff(base::as.numeric(table$timestamp))
    dt <- dt[is.finite(dt) & dt > 0]
    if (length(dt) == 0) {
      msg <- paste0(
        "Skipped file '", base::basename(file_path),
        "': cannot determine sampling interval (too few valid timestamps)."
      )
      warning(msg)
      return(list(status = "skip_bad_interval", file = base::basename(file_path), msg = msg))
    }
    interval <- pracma::Mode(dt)
    if (!is.finite(interval) || interval <= 0) {
      msg <- paste0(
        "Skipped file '", base::basename(file_path),
        "': invalid sampling interval detected (", interval, ")."
      )
      warning(msg)
      return(list(status = "skip_bad_interval", file = base::basename(file_path), msg = msg))
    }

    # ---------------- Step 5: Skip initial hours (optional) -------------------
    if (skiphours > 0) {

      # If file shorter than skiphours, skip the whole file
      t0 <- min(table$timestamp, na.rm = TRUE)
      t1 <- max(table$timestamp, na.rm = TRUE)
      duration_hours_total <- as.numeric(difftime(t1, t0, units = "hours"))

      if (!is.finite(duration_hours_total) || duration_hours_total <= skiphours) {
        msg <- paste0(
          "Skipped file '", base::basename(file_path),
          "': total duration (", round(duration_hours_total, 2),
          " h) <= skiphours (", skiphours, " h)."
        )
        warning(msg)
        return(list(status = "skip_skiphours_too_short", file = base::basename(file_path), msg = msg))
      }

      # Otherwise skip the first skiphours
      skip_seconds <- skiphours * 3600
      skip_time <- as.numeric(t0) + skip_seconds

      ts_num <- as.numeric(table$timestamp)
      skip_row <- which.min(abs(ts_num - skip_time))

      if (!is.na(skip_row) && skip_row > 0) {
        table <- table[-c(1:skip_row), , drop = FALSE]
      }

      # Defensive: if skipping makes table empty, skip
      if (base::nrow(table) == 0) {
        msg <- paste0(
          "Skipped file '", base::basename(file_path),
          "': empty after applying skiphours (", skiphours, " h)."
        )
        warning(msg)
        return(list(status = "skip_skiphours_too_short", file = base::basename(file_path), msg = msg))
      }
    }

    # --------------------- Step 6: Handle gaps robustly ----------------------
    if (isTRUE(gapprocess)) {

      maxgap_points <- base::max(1, base::floor(maxgapmin * 60 / interval))

      first_valid <- base::which(!base::is.na(table$sensorglucose))[1]
      if (!base::is.na(first_valid) && first_valid > 1) {
        table <- table[first_valid:base::nrow(table), , drop = FALSE]
      }

      if (shortgapsm == "linear") {
        table$sensorglucose <- zoo::na.approx(
          table$sensorglucose,
          na.rm = FALSE,
          maxgap = maxgap_points
        )
      } else if (shortgapsm == "spline") {
        table$sensorglucose <- zoo::na.spline(
          table$sensorglucose,
          na.rm = FALSE,
          maxgap = maxgap_points
        )
      } else {
        stop("Invalid shortgapsm. Choose 'linear' or 'spline'.")
      }

      glucose <- table$sensorglucose
      r <- base::rle(base::is.na(glucose))
      ends <- base::cumsum(r$lengths)
      starts <- ends - r$lengths + 1
      gap_runs <- base::which(r$values)
      gap_starts <- starts[gap_runs]
      gap_ends <- ends[gap_runs]

      remove_idx_accum <- base::integer(0)

      for (i in base::seq_along(gap_starts)) {
        gap_idx <- base::seq(gap_starts[i], gap_ends[i])
        gap_len_points <- base::length(gap_idx)
        if (gap_len_points <= maxgap_points) next

        if (longgapsm == "dailymean") {
          times_hms <- base::format(table$timestamp, "%H:%M:%S")
          non_na_idx <- base::which(!base::is.na(table$sensorglucose))
          if (base::length(non_na_idx) > 0) {
            daily_means <- base::tapply(
              table$sensorglucose[non_na_idx],
              times_hms[non_na_idx],
              base::mean,
              na.rm = TRUE
            )
            gap_times <- times_hms[gap_idx]
            to_replace <- base::which(gap_times %in% base::names(daily_means))
            if (base::length(to_replace) > 0) {
              table$sensorglucose[gap_idx[to_replace]] <- daily_means[gap_times[to_replace]]
            }
          }

        } else if (longgapsm == "spline") {
          prev_idx <- gap_idx[1] - 1
          next_idx <- gap_idx[base::length(gap_idx)] + 1

          can_interp <- (
            prev_idx >= 1 &&
              next_idx <= base::nrow(table) &&
              !base::is.na(table$sensorglucose[prev_idx]) &&
              !base::is.na(table$sensorglucose[next_idx])
          )

          if (can_interp) {
            x_anchor <- base::as.numeric(table$timestamp[c(prev_idx, next_idx)])
            y_anchor <- table$sensorglucose[c(prev_idx, next_idx)]
            xout <- base::as.numeric(table$timestamp[gap_idx])
            interp_vals <- stats::spline(x = x_anchor, y = y_anchor, xout = xout)$y
            table$sensorglucose[gap_idx] <- interp_vals
          } else {
            times_hms <- base::format(table$timestamp, "%H:%M:%S")
            non_na_idx <- base::which(!base::is.na(table$sensorglucose))
            if (base::length(non_na_idx) > 0) {
              daily_means <- base::tapply(
                table$sensorglucose[non_na_idx],
                times_hms[non_na_idx],
                base::mean,
                na.rm = TRUE
              )
              gap_times <- times_hms[gap_idx]
              to_replace <- base::which(gap_times %in% base::names(daily_means))
              if (base::length(to_replace) > 0) {
                table$sensorglucose[gap_idx[to_replace]] <- daily_means[gap_times[to_replace]]
              }
            }
          }

        } else if (longgapsm == "remove") {
          remove_idx_accum <- base::c(remove_idx_accum, gap_idx)

        } else {
          stop("Invalid longgapsm. Choose 'dailymean', 'spline', or 'remove'.")
        }
      }

      if (base::length(remove_idx_accum) > 0) {
        remove_idx_accum <- base::sort(base::unique(remove_idx_accum))
        table <- table[-remove_idx_accum, , drop = FALSE]
        if (base::nrow(table) == 0) stop("All data removed after long gap removal.")
      }
    }

    # -------------------- Step 7: Trim to full 24-hour days -------------------
    if (isTRUE(trimdays)) {

      start_time <- min(table$timestamp, na.rm = TRUE)
      end_time <- max(table$timestamp, na.rm = TRUE)

      duration_hours <- as.numeric(difftime(end_time, start_time, units = "hours"))
      if (!is.finite(duration_hours) || duration_hours < 24) {
        msg <- paste0(
          "Skipped file '", base::basename(file_path),
          "': less than one full 24-hour day of data (trimdays=TRUE)."
        )
        warning(msg)
        return(list(status = "skip_less_than_1day", file = base::basename(file_path), msg = msg))
      }

      start_midnight <- as.POSIXct(format(start_time, "%Y-%m-%d 00:00:00"), tz = "UTC")
      end_midnight <- as.POSIXct(format(end_time, "%Y-%m-%d 00:00:00"), tz = "UTC")

      keep_start <- start_midnight + 86400  # drop partial leading day
      keep_end <- end_midnight              # up to last midnight

      table <- table[table$timestamp >= keep_start & table$timestamp < keep_end, ]

      if (base::nrow(table) == 0) {
        msg <- paste0(
          "Skipped file '", base::basename(file_path),
          "': trimming produced empty dataset."
        )
        warning(msg)
        return(list(status = "skip_trim_empty", file = base::basename(file_path), msg = msg))
      }
    }

    # -------------------- Step 8: Finalize and write CSV ----------------------
    table$subjectid <- id
    table$subjectid <- as.character(table$subjectid)
    table <- table[, c("subjectid", "timestamp", "sensorglucose")]

    table$timestamp <- as.POSIXct(table$timestamp, tz = "UTC")
    options(digits.secs = 0)

    filename <- base::paste0(
      outputdir, "/", tools::file_path_sans_ext(base::basename(file_path)), ".csv"
    )
    utils::write.csv(table, file = filename, row.names = FALSE, quote = FALSE)

    return(list(status = "ok", file = base::basename(file_path), msg = ""))
  }

  # ----------------------- Process each file ---------------------------------
  for (f in seq_along(files)) {
    file_path <- files[f]

    res <- base::tryCatch(
      process_one(file_path),
      error = function(e) {
        msg <- paste0(
          "Skipping file '", base::basename(file_path),
          "' due to error: ", e$message
        )
        warning(msg)
        list(status = "error", file = base::basename(file_path), msg = msg)
      }
    )

    if (!is.list(res) || is.null(res$status)) next

    if (identical(res$status, "skip_few_glucose")) {
      skipped_few_glucose <- c(skipped_few_glucose, res$file)
      next
    }
    if (identical(res$status, "skip_trim_empty")) {
      skipped_trim_empty <- c(skipped_trim_empty, res$file)
      next
    }
    if (identical(res$status, "skip_less_than_1day")) {
      skipped_less_than_1day <- c(skipped_less_than_1day, res$file)
      next
    }
    if (identical(res$status, "skip_bad_interval")) {
      skipped_bad_interval <- c(skipped_bad_interval, res$file)
      next
    }
    if (identical(res$status, "skip_skiphours_too_short")) {
      skipped_skiphours_too_short <- c(skipped_skiphours_too_short, res$file)
      next
    }
    if (identical(res$status, "error")) {
      error_files <- c(error_files, res$file)
      next
    }
  }

  # ----------------------- Summary: print skipped files -----------------------
  skipped_few_glucose <- unique(skipped_few_glucose)
  skipped_trim_empty <- unique(skipped_trim_empty)
  skipped_less_than_1day <- unique(skipped_less_than_1day)
  skipped_bad_interval <- unique(skipped_bad_interval)
  skipped_skiphours_too_short <- unique(skipped_skiphours_too_short)
  error_files <- unique(error_files)

  if (length(skipped_few_glucose) > 0) {
    message("\nFiles skipped due to too few glucose values (<10):")
    message(paste0(" - ", skipped_few_glucose, collapse = "\n"))
  }
  if (length(skipped_trim_empty) > 0) {
    message("\nFiles skipped because trimming produced empty dataset:")
    message(paste0(" - ", skipped_trim_empty, collapse = "\n"))
  }
  if (length(skipped_less_than_1day) > 0) {
    message("\nFiles skipped because they cover less than one full 24-hour day (trimdays = TRUE):")
    message(paste0(" - ", skipped_less_than_1day, collapse = "\n"))
  }
  if (length(skipped_bad_interval) > 0) {
    message("\nFiles skipped due to invalid/undeterminable sampling interval:")
    message(paste0(" - ", skipped_bad_interval, collapse = "\n"))
  }
  if (length(skipped_skiphours_too_short) > 0) {
    message("\nFiles skipped because total duration <= skiphours:")
    message(paste0(" - ", skipped_skiphours_too_short, collapse = "\n"))
  }
  if (length(error_files) > 0) {
    message("\nFiles skipped due to errors:")
    message(paste0(" - ", error_files, collapse = "\n"))
  }

  invisible(list(
    skipped_few_glucose = skipped_few_glucose,
    skipped_trim_empty = skipped_trim_empty,
    skipped_less_than_1day = skipped_less_than_1day,
    skipped_bad_interval = skipped_bad_interval,
    skipped_skiphours_too_short = skipped_skiphours_too_short,
    skipped_errors = error_files
  ))
}
