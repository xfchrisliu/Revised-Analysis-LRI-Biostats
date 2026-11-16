
# ============================================================
# Diet x Microbiome: Normality + Zero Check
# UPDATED FOR: NEW STAT DATA TABLE.xlsx
#
# Sheet types:
#   1. Germ Status + Diet + One Way Group + outcomes
#      -> Two-way model residual Shapiro-Wilk test
#      -> Two-way zero-pattern rules
#
#   2. Diet + outcomes
#      -> One-way model residual Shapiro-Wilk test
#      -> One-way zero-pattern rules
#
#   3. Percentage outcomes (e.g., Steatosis %)
#      -> Identified and flagged
#      -> Recommended transformation = Raw
#
# Zero rules:
#   Yes:
#     - All observations are zero
#     - Only 1-2 nonzero observations overall
#     - Entire Germ Status level is zero (two-way only)
#     - Entire Diet level is zero
#     - Only 1 nonzero observation in a Germ Status x Diet
#       cell (two-way)
#     - Only 1 nonzero observation within a Diet group
#       (one-way; analogous screening rule)
#
#   No:
#     - No zero-rule violation
#     - A two-way Germ Status x Diet cell is all zero
#       (retained from the original rule)
#
# Shapiro-Wilk:
#   - Attempted regardless of zero-rule result
#   - Uses residuals from the raw-scale model
#
# Outputs:
#   Diet_Microbiome_Normality_Zero_Check_All_Figures.csv
#   Diet_Microbiome_Normality_Zero_Check_By_Sheet.xlsx
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(gridExtra)

options(
  contrasts = c("contr.sum", "contr.poly")
)


# ============================================================
# 2. Input file
# ============================================================

input_file <- "NEW STAT update.xlsx"

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

sheet_names <- excel_sheets(input_file)

cat("Sheets found:\n")
print(sheet_names)


# ============================================================
# 3. Helper functions
# ============================================================

normalize_name <- function(x) {
  x <- tolower(trimws(x))
  gsub("[^a-z0-9]", "", x)
}


find_column <- function(data, candidates) {

  nm <- normalize_name(names(data))
  candidates <- normalize_name(candidates)

  idx <- match(candidates, nm, nomatch = 0)
  idx <- idx[idx > 0]

  if (length(idx) == 0) {
    return(NA_character_)
  }

  names(data)[idx[1]]
}


is_percentage_outcome <- function(outcome) {

  grepl(
    "%|percent|percentage",
    tolower(outcome)
  )
}


# Safely convert a response column to numeric
to_numeric <- function(x) {

  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  suppressWarnings(
    as.numeric(trimws(as.character(x)))
  )
}


# ============================================================
# 4. Zero-pattern check: TWO-WAY data
# ============================================================

check_zero_pattern_two_way <- function(d) {

  n <- nrow(d)
  z <- sum(d$VALUE == 0)
  nz <- sum(d$VALUE != 0)

  zp <- if (n > 0) 100 * z / n else NA_real_

  base <- list(
    n = n,
    z = z,
    nz = nz,
    zp = zp,
    too_many_zero = "Yes",
    reason = "NO_VALID_DATA",
    zero_cells = NA_integer_,
    one_nonzero_cells = NA_integer_,
    germ_all_zero = NA_character_,
    diet_all_zero = NA_character_
  )

  if (n == 0) {
    return(base)
  }

  # All observations zero
  if (nz == 0) {
    base$reason <- "ALL_ZERO"
    base$germ_all_zero <- "ALL"
    base$diet_all_zero <- "ALL"
    return(base)
  }

  # Only 1-2 nonzero observations overall
  if (nz <= 2) {
    base$reason <- "ONLY_1_2_NONZERO_OVERALL"
    return(base)
  }

  # Entire Germ Status level is zero
  gs <- d %>%
    group_by(GERM_STATUS) %>%
    summarise(
      nonzero_n = sum(VALUE != 0),
      .groups = "drop"
    )

  ga <- as.character(
    gs$GERM_STATUS[gs$nonzero_n == 0]
  )

  if (length(ga) > 0) {
    base$reason <- paste0(
      "GERM_STATUS_LEVEL_ALL_ZERO:",
      paste(ga, collapse = ";")
    )
    base$germ_all_zero <- paste(
      ga, collapse = ";"
    )
    return(base)
  }

  # Entire Diet level is zero
  ds <- d %>%
    group_by(DIET) %>%
    summarise(
      nonzero_n = sum(VALUE != 0),
      .groups = "drop"
    )

  da <- as.character(
    ds$DIET[ds$nonzero_n == 0]
  )

  if (length(da) > 0) {
    base$reason <- paste0(
      "DIET_LEVEL_ALL_ZERO:",
      paste(da, collapse = ";")
    )
    base$diet_all_zero <- paste(
      da, collapse = ";"
    )
    return(base)
  }

  # Check Germ Status x Diet cells
  cs <- d %>%
    group_by(GERM_STATUS, DIET) %>%
    summarise(
      nonzero_n = sum(VALUE != 0),
      .groups = "drop"
    )

  zero_cells <- sum(cs$nonzero_n == 0)

  one <- cs %>%
    filter(nonzero_n == 1)

  base$zero_cells <- zero_cells
  base$one_nonzero_cells <- nrow(one)

  # Exactly one nonzero observation in a cell
  if (nrow(one) > 0) {

    labs <- paste(
      one$GERM_STATUS,
      one$DIET,
      sep = " + "
    )

    base$reason <- paste0(
      "ONE_NONZERO_IN_CELL:",
      paste(labs, collapse = ";")
    )

    return(base)
  }

  # A completely zero cell is allowed by the original rule
  base$too_many_zero <- "No"

  base$reason <- if (zero_cells > 0) {
    "OK_WITH_ALL_ZERO_CELL(S)"
  } else {
    "OK"
  }

  return(base)
}


# ============================================================
# 5. Zero-pattern check: ONE-WAY data
# ============================================================

check_zero_pattern_one_way <- function(d) {

  n <- nrow(d)
  z <- sum(d$VALUE == 0)
  nz <- sum(d$VALUE != 0)

  zp <- if (n > 0) 100 * z / n else NA_real_

  base <- list(
    n = n,
    z = z,
    nz = nz,
    zp = zp,
    too_many_zero = "Yes",
    reason = "NO_VALID_DATA",
    zero_cells = NA_integer_,
    one_nonzero_cells = NA_integer_,
    germ_all_zero = NA_character_,
    diet_all_zero = NA_character_
  )

  if (n == 0) {
    return(base)
  }

  if (nz == 0) {
    base$reason <- "ALL_ZERO"
    base$diet_all_zero <- "ALL"
    return(base)
  }

  if (nz <= 2) {
    base$reason <- "ONLY_1_2_NONZERO_OVERALL"
    return(base)
  }

  # Entire Diet group is zero
  ds <- d %>%
    group_by(DIET) %>%
    summarise(
      nonzero_n = sum(VALUE != 0),
      .groups = "drop"
    )

  da <- as.character(
    ds$DIET[ds$nonzero_n == 0]
  )

  if (length(da) > 0) {
    base$reason <- paste0(
      "DIET_LEVEL_ALL_ZERO:",
      paste(da, collapse = ";")
    )
    base$diet_all_zero <- paste(
      da, collapse = ";"
    )
    return(base)
  }

  # Analogous screening rule:
  # only one nonzero observation in a Diet group
  one <- ds %>%
    filter(nonzero_n == 1)

  base$one_nonzero_cells <- nrow(one)

  if (nrow(one) > 0) {
    base$reason <- paste0(
      "ONE_NONZERO_IN_DIET:",
      paste(as.character(one$DIET), collapse = ";")
    )
    return(base)
  }

  base$too_many_zero <- "No"
  base$reason <- "OK"

  return(base)
}


# ============================================================
# 6. Check one outcome
# ============================================================

check_outcome <- function(
    dat,
    outcome,
    sheet_name,
    analysis_type,
    germ_col,
    diet_col
) {

  percentage <- is_percentage_outcome(outcome)

  # Prepare outcome and metadata
  d <- data.frame(
    VALUE = to_numeric(dat[[outcome]])
  )

  d$DIET <- as.character(dat[[diet_col]])

  if (analysis_type == "Two-way ANOVA") {
    d$GERM_STATUS <- as.character(dat[[germ_col]])
  }

  # Remove incomplete observations
  if (analysis_type == "Two-way ANOVA") {

    d <- d %>%
      filter(
        !is.na(VALUE),
        !is.na(GERM_STATUS),
        trimws(GERM_STATUS) != "",
        !is.na(DIET),
        trimws(DIET) != ""
      )

    d$GERM_STATUS <- factor(d$GERM_STATUS)
    d$DIET <- factor(d$DIET)

    d <- droplevels(d)

    q <- check_zero_pattern_two_way(d)

  } else {

    d <- d %>%
      filter(
        !is.na(VALUE),
        !is.na(DIET),
        trimws(DIET) != ""
      )

    d$DIET <- factor(d$DIET)

    d <- droplevels(d)

    q <- check_zero_pattern_one_way(d)
  }


  # Common output fields
  result <- data.frame(
    Sheet = sheet_name,
    Outcome = outcome,
    Outcome_Type = if (percentage) {
      "Percentage"
    } else {
      "Numeric"
    },
    Analysis_Type = analysis_type,

    N = q$n,
    Zero_N = q$z,
    Nonzero_N = q$nz,
    Zero_Percent = q$zp,

    `too many zero?` = q$too_many_zero,

    Zero_Cell_Count = q$zero_cells,
    One_Nonzero_Cell_Count = q$one_nonzero_cells,

    GermStatus_All_Zero = q$germ_all_zero,
    Diet_All_Zero = q$diet_all_zero,

    Analysis_Status = "SKIP",
    Analysis_Reason = q$reason,

    Raw_Shapiro_W = NA_real_,
    Raw_Shapiro_p = NA_real_,

    Recommended_Transformation = NA_character_,

    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  # ----------------------------------------------------------
  # Factor-level checks
  # ----------------------------------------------------------

  if (analysis_type == "Two-way ANOVA") {

    if (
      nlevels(d$GERM_STATUS) < 2 ||
      nlevels(d$DIET) < 2
    ) {
      result$Analysis_Reason <- paste0(
        q$reason,
        "; INSUFFICIENT_FACTOR_LEVELS"
      )
      return(result)
    }

    formula <- VALUE ~ GERM_STATUS * DIET

  } else {

    if (nlevels(d$DIET) < 2) {
      result$Analysis_Reason <- paste0(
        q$reason,
        "; INSUFFICIENT_DIET_LEVELS"
      )
      return(result)
    }

    formula <- VALUE ~ DIET
  }


  # ----------------------------------------------------------
  # Fit raw-scale model regardless of zero-rule result
  # ----------------------------------------------------------

  fit <- tryCatch(
    lm(formula, data = d),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    result$Analysis_Reason <- paste0(
      q$reason,
      "; RAW_MODEL_FAILED"
    )
    return(result)
  }


  # ----------------------------------------------------------
  # Shapiro-Wilk test on model residuals
  #
  # Attempt regardless of "too many zero?".
  # Shapiro-Wilk requires 3-5000 observations.
  # ----------------------------------------------------------

  r <- residuals(fit)

# ----------------------------------------------------------
# Residual Q-Q plot
# ----------------------------------------------------------

if (
  length(r) >= 3 &&
  length(unique(r)) > 2 &&
  all(is.finite(r))
) {

  qq_dat <- data.frame(
    Residuals = r
  )

  p <- ggplot(
    qq_dat,
    aes(sample = Residuals)
  ) +

    stat_qq(
      size = 1.2
    ) +

    stat_qq_line(
      linewidth = 0.6
    ) +

    labs(
      title = paste(sheet_name, outcome, sep = ": "),
      subtitle = paste("Residual n =", length(r)),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +

    theme_bw(
      base_size = 10
    ) +

    theme(

      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 10
      ),

      plot.subtitle = element_text(
        hjust = 0.5,
        size = 8
      ),

      axis.title = element_text(
        face = "bold",
        size = 8
      ),

      axis.text = element_text(
        size = 7
      ),

      plot.margin = margin(
        5, 5, 5, 5
      )
    )

  # Save plot
  qq_plot_list[[length(qq_plot_list) + 1]] <<- p

}


  sw <- NULL

  if (
    length(r) >= 3 &&
    length(r) <= 5000 &&
    length(unique(r)) > 2 &&
    all(is.finite(r))
  ) {

    sw <- tryCatch(
      shapiro.test(r),
      error = function(e) NULL
    )
  }


  # ----------------------------------------------------------
  # Recommend transformation
  # ----------------------------------------------------------

  if (percentage) {

    # Preserve percentages on their original scale.
    result$Recommended_Transformation <- "Raw"

  } else if (is.null(sw)) {

    result$Recommended_Transformation <- NA_character_

  } else if (sw$p.value >= 0.05) {

    result$Recommended_Transformation <- "Raw"

  } else if (all(d$VALUE >= 0)) {

    result$Recommended_Transformation <- "Log10_plus1"

  } else {

    result$Recommended_Transformation <- "Raw"

    result$Analysis_Reason <- paste0(
      result$Analysis_Reason,
      "; LOG_TRANSFORMATION_NOT_APPLICABLE_NEGATIVE_VALUES"
    )
  }


  # ----------------------------------------------------------
  # Record Shapiro results
  # ----------------------------------------------------------

  result$Analysis_Status <- "OK"

  if (!is.null(sw)) {

    result$Raw_Shapiro_W <- unname(sw$statistic)
    result$Raw_Shapiro_p <- sw$p.value

  } else {

    result$Analysis_Reason <- paste0(
      result$Analysis_Reason,
      "; SHAPIRO_NOT_AVAILABLE"
    )
  }

  return(result)
}


# ============================================================
# 7. Process every sheet
# ============================================================

# Store residual Q-Q plots
qq_plot_list <- list()
all_results <- list()
idx <- 1


for (sheet in sheet_names) {

  cat("\nProcessing sheet:", sheet, "\n")

  dat <- read_excel(
    input_file,
    sheet = sheet
  )

  # Remove completely empty rows
  dat <- dat[
    rowSums(!is.na(dat)) > 0,
    ,
    drop = FALSE
  ]

  if (ncol(dat) < 2) {
    warning("Too few columns in sheet: ", sheet)
    next
  }


  # Identify metadata columns by name
  germ_col <- find_column(
    dat,
    c("Germ Status", "Germ_Status", "GermStatus")
  )

  diet_col <- find_column(
    dat,
    c("Diet")
  )

  group_col <- find_column(
    dat,
    c("One way group", "One_way_group", "OneWayGroup")
  )


  if (is.na(diet_col)) {
    warning("Diet column not found in sheet: ", sheet)
    next
  }


  # Determine model type
  if (!is.na(germ_col)) {

    analysis_type <- "Two-way ANOVA"

    metadata_cols <- c(
      germ_col,
      diet_col,
      group_col
    )

    metadata_cols <- metadata_cols[
      !is.na(metadata_cols)
    ]

  } else {

    analysis_type <- "One-way ANOVA"

    metadata_cols <- diet_col
  }


  # Remaining columns are candidate outcomes
  outcome_cols <- setdiff(
    names(dat),
    metadata_cols
  )


  # Skip columns that contain no numeric values
  for (outcome in outcome_cols) {

    numeric_values <- to_numeric(dat[[outcome]])

    if (all(is.na(numeric_values))) {
      cat("  Skipping nonnumeric column:", outcome, "\n")
      next
    }

    res <- tryCatch(

      check_outcome(
        dat = dat,
        outcome = outcome,
        sheet_name = sheet,
        analysis_type = analysis_type,
        germ_col = germ_col,
        diet_col = diet_col
      ),

      error = function(e) {

        data.frame(
          Sheet = sheet,
          Outcome = outcome,
          Outcome_Type = if (
            is_percentage_outcome(outcome)
          ) "Percentage" else "Numeric",
          Analysis_Type = analysis_type,

          N = NA_integer_,
          Zero_N = NA_integer_,
          Nonzero_N = NA_integer_,
          Zero_Percent = NA_real_,

          `too many zero?` = NA_character_,

          Zero_Cell_Count = NA_integer_,
          One_Nonzero_Cell_Count = NA_integer_,

          GermStatus_All_Zero = NA_character_,
          Diet_All_Zero = NA_character_,

          Analysis_Status = "SKIP",
          Analysis_Reason = paste0("ERROR: ", e$message),

          Raw_Shapiro_W = NA_real_,
          Raw_Shapiro_p = NA_real_,

          Recommended_Transformation = NA_character_,

          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    )

    all_results[[idx]] <- res
    idx <- idx + 1
  }
}


# ============================================================
# 8. Combine results
# ============================================================

if (length(all_results) == 0) {
  stop("No numeric outcomes were found.")
}

final_results <- bind_rows(all_results)

# ============================================================
# 8.5 Generate Residual Q-Q Plot PDF
#
# Layout:
#   4 plots per row
#   4 rows per page
#   16 plots per page
# ============================================================

if (length(qq_plot_list) > 0) {

  plots_per_row <- 4
  rows_per_page <- 4

  plots_per_page <- plots_per_row * rows_per_page

  n_plots <- length(qq_plot_list)

  n_pages <- ceiling(
    n_plots / plots_per_page
  )

  pdf_path <- paste0(
    "Diet_Microbiome_Residual_QQPlots.pdf"
  )

  pdf(
    file = pdf_path,
    width = 14,
    height = 11
  )

  for (page in seq_len(n_pages)) {

    start_index <- (
      (page - 1) * plots_per_page
    ) + 1

    end_index <- min(
      page * plots_per_page,
      n_plots
    )

    page_plots <- qq_plot_list[
      start_index:end_index
    ]

    grid.arrange(

      grobs = page_plots,

      ncol = plots_per_row,

      nrow = rows_per_page,

      top = paste(
        "Diet x Microbiome",
        "- Model Residual Q-Q Plots"
      )
    )
  }

  dev.off()

  cat("\nResidual Q-Q plots saved:\n")
  cat(pdf_path, "\n")
  cat("Number of plots:", n_plots, "\n")
  cat("Number of pages:", n_pages, "\n")

} else {

  cat("\nNo valid residual Q-Q plots were generated.\n")

}

# ============================================================
# 9. Save combined CSV
# ============================================================

write.csv(
  final_results,
  "Diet_Microbiome_Normality_Zero_Check_All_Figures.csv",
  row.names = FALSE,
  na = ""
)


# ============================================================
# 10. Save Excel workbook
# ============================================================

out_wb <- createWorkbook()

addWorksheet(out_wb, "All Results")
writeData(out_wb, "All Results", final_results)

# Separate summary for percentage outcomes
addWorksheet(out_wb, "Percentage Outcomes")

percentage_results <- final_results %>%
  filter(Outcome_Type == "Percentage")

writeData(
  out_wb,
  "Percentage Outcomes",
  percentage_results
)

# One result sheet per original input sheet
for (sheet in unique(final_results$Sheet)) {

  output_sheet <- substr(
    gsub("[\\[\\]\\*\\?/\\\\:]", "_", sheet),
    1,
    31
  )

  base_name <- output_sheet
  counter <- 1

  while (output_sheet %in% names(out_wb)) {

    suffix <- paste0("_", counter)

    output_sheet <- paste0(
      substr(base_name, 1, 31 - nchar(suffix)),
      suffix
    )

    counter <- counter + 1
  }

  addWorksheet(out_wb, output_sheet)

  writeData(
    out_wb,
    output_sheet,
    filter(final_results, Sheet == sheet)
  )
}


saveWorkbook(
  out_wb,
  "Diet_Microbiome_Normality_Zero_Check_By_Sheet.xlsx",
  overwrite = TRUE
)


# ============================================================
# 11. Print summary
# ============================================================

cat("\n============================================\n")
cat("Normality + zero check complete\n")
cat("============================================\n")

cat("Input sheets: ", length(sheet_names), "\n")
cat("Total outcomes: ", nrow(final_results), "\n")

cat(
  "Too many zero = YES: ",
  sum(final_results$`too many zero?` == "Yes", na.rm = TRUE),
  "\n"
)

cat(
  "Too many zero = NO: ",
  sum(final_results$`too many zero?` == "No", na.rm = TRUE),
  "\n"
)

cat(
  "Shapiro p-values available: ",
  sum(!is.na(final_results$Raw_Shapiro_p)),
  "\n"
)

cat(
  "Percentage outcomes: ",
  sum(final_results$Outcome_Type == "Percentage", na.rm = TRUE),
  "\n"
)

cat("\nFiles saved:\n")
cat("1. Diet_Microbiome_Normality_Zero_Check_All_Figures.csv\n")
cat("2. Diet_Microbiome_Normality_Zero_Check_By_Sheet.xlsx\n")

print(final_results)
