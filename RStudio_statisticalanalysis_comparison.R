# Statistical analysis of binding affinities
# Three separate pairwise comparison tables
#
# Comparisons:
# 1. Reference JAK inhibitors vs Flavonoids
# 2. Other phytochemicals vs Reference JAK inhibitors
# 3. Flavonoids vs Other phytochemicals
#
# Binding-affinity values are reproduced exactly from the previously
# provided R analysis file.

required_packages <- c("dplyr", "purrr", "readr")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

# -------------------- 1. Input data --------------------

data <- data.frame(
  Ligand = c(
    "Baricitinib", "Caffeic_acid", "Catechin", "Chlorogenic_acid",
    "Coumaric_acid", "Gallic_acid", "Gingerol", "Isoliquiritigenin",
    "Kaempferol", "Naringenin", "Peficitinib", "Quercetin", "Tofacitinib"
  ),
  Group = c(
    "Reference_JAKi", "Other_phytochemical", "Flavonoid",
    "Other_phytochemical", "Other_phytochemical", "Other_phytochemical",
    "Other_phytochemical", "Other_phytochemical", "Flavonoid",
    "Flavonoid", "Reference_JAKi", "Flavonoid", "Reference_JAKi"
  ),
  JAK1 = c(
    -8.1, -6.2, -8.3, -8.1, -5.7, -5.6, -6.2,
    -7.3, -7.7, -7.8, -9.0, -7.9, -7.8
  ),
  JAK2 = c(
    -8.2, -6.5, -8.8, -8.1, -6.4, -6.2, -5.7,
    -8.0, -8.9, -8.6, -7.6, -9.3, -7.4
  ),
  JAK3 = c(
    -7.4, -6.1, -8.8, -7.5, -6.0, -5.4, -5.7,
    -8.1, -8.5, -9.0, -9.0, -8.5, -7.6
  ),
  TYK2 = c(
    -6.6, -5.4, -7.4, -7.2, -5.1, -5.2, -5.5,
    -6.4, -7.2, -7.0, -6.9, -7.2, -6.2
  ),
  stringsAsFactors = FALSE
)

data$Group <- factor(
  data$Group,
  levels = c("Reference_JAKi", "Flavonoid", "Other_phytochemical")
)

isoforms <- c("JAK1", "JAK2", "JAK3", "TYK2")

stopifnot(
  nrow(data) == 13,
  all(isoforms %in% names(data)),
  all(!is.na(data[, c("Ligand", "Group", isoforms)]))
)

# -------------------- 2. Rank-biserial effect size --------------------

rank_biserial <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  rank_sum_x <- sum(ranks[seq_len(n1)])
  U1 <- rank_sum_x - n1 * (n1 + 1) / 2
  2 * U1 / (n1 * n2) - 1
}

# -------------------- 3. Function for one comparison --------------------

run_comparison <- function(comparison_name, group1, group2) {

  map_dfr(isoforms, function(iso) {

    x <- data[data$Group == group1, iso, drop = TRUE]
    y <- data[data$Group == group2, iso, drop = TRUE]

    test <- suppressWarnings(
      wilcox.test(
        x, y,
        alternative = "two.sided",
        exact = TRUE,
        conf.int = TRUE,
        conf.level = 0.95
      )
    )

    hl <- if (!is.null(test$estimate)) unname(test$estimate) else NA_real_
    ci_low <- if (length(test$conf.int) >= 2) test$conf.int[1] else NA_real_
    ci_high <- if (length(test$conf.int) >= 2) test$conf.int[2] else NA_real_

    med1 <- median(x)
    med2 <- median(y)

    data.frame(
      Comparison = comparison_name,
      Group_1 = group1,
      Group_2 = group2,
      Isoform = iso,
      N_Group_1 = length(x),
      N_Group_2 = length(y),
      Median_Group_1 = med1,
      IQR_Group_1 = IQR(x),
      Median_Group_2 = med2,
      IQR_Group_2 = IQR(y),
      Median_Difference_Group1_minus_Group2 = med1 - med2,
      W = unname(test$statistic),
      P_value = test$p.value,
      Hodges_Lehmann_Difference = hl,
      HL_CI_95_Lower = ci_low,
      HL_CI_95_Upper = ci_high,
      Rank_Biserial_r = rank_biserial(x, y),
      stringsAsFactors = FALSE
    )
  })
}

# -------------------- 4. Three separate analyses --------------------
# These are intentionally kept as three separate result tables.
# Holm correction is applied across all 12 prespecified tests.

table1 <- run_comparison(
  "Reference JAK inhibitors vs Flavonoids",
  "Reference_JAKi", "Flavonoid"
)

table2 <- run_comparison(
  "Other phytochemicals vs Reference JAK inhibitors",
  "Other_phytochemical", "Reference_JAKi"
)

table3 <- run_comparison(
  "Flavonoids vs Other phytochemicals",
  "Flavonoid", "Other_phytochemical"
)

# Holm correction across 3 comparisons x 4 isoforms = 12 tests
all_p <- c(table1$P_value, table2$P_value, table3$P_value)
all_p_holm <- p.adjust(all_p, method = "holm")

table1$P_adjusted_Holm <- all_p_holm[1:4]
table2$P_adjusted_Holm <- all_p_holm[5:8]
table3$P_adjusted_Holm <- all_p_holm[9:12]

add_labels <- function(x) {
  x %>%
    mutate(
      Median_Direction = case_when(
        Median_Difference_Group1_minus_Group2 < 0 ~
          "Group 1 more negative",
        Median_Difference_Group1_minus_Group2 > 0 ~
          "Group 2 more negative",
        TRUE ~ "Equal medians"
      ),
      Significant_Holm_0.05 = ifelse(
        P_adjusted_Holm < 0.05, "Yes", "No"
      )
    )
}

table1 <- add_labels(table1)
table2 <- add_labels(table2)
table3 <- add_labels(table3)

# -------------------- 5. Display three separate tables --------------------

cat("\n====================================================\n")
cat("TABLE 1: Reference JAK inhibitors vs Flavonoids\n")
cat("====================================================\n")
print(table1)

cat("\n====================================================\n")
cat("TABLE 2: Other phytochemicals vs Reference JAK inhibitors\n")
cat("====================================================\n")
print(table2)

cat("\n====================================================\n")
cat("TABLE 3: Flavonoids vs Other phytochemicals\n")
cat("====================================================\n")
print(table3)

# -------------------- 6. Save three separate CSV files --------------------

dir.create("results", showWarnings = FALSE)
dir.create("results/Table_1_Reference_JAKi_vs_Flavonoids",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/Table_2_Other_phytochemicals_vs_Reference_JAKi",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/Table_3_Flavonoids_vs_Other_phytochemicals",
           recursive = TRUE, showWarnings = FALSE)

write_csv(
  table1,
  "results/Table_1_Reference_JAKi_vs_Flavonoids/Table_1_Reference_JAKi_vs_Flavonoids.csv"
)

write_csv(
  table2,
  "results/Table_2_Other_phytochemicals_vs_Reference_JAKi/Table_2_Other_phytochemicals_vs_Reference_JAKi.csv"
)

write_csv(
  table3,
  "results/Table_3_Flavonoids_vs_Other_phytochemicals/Table_3_Flavonoids_vs_Other_phytochemicals.csv"
)

# Save the exact input values separately for reproducibility.
write_csv(data, "results/input_binding_affinity_data.csv")

cat("\nAnalysis completed. Three separate tables were generated.\n")

getwd()