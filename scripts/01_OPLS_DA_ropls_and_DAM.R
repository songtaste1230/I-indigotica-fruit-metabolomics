############################################################
# Script: 01_OPLS_DA_ropls_and_DAM.R
#
# Purpose:
#   Pairwise OPLS-DA analysis, VIP extraction, and
#   differentially accumulated metabolite (DAM) statistics
#   for widely targeted metabolomics data from Isatis indigotica fruit.
#
# Manuscript:
#   Stage-resolved multi-omics links flavonoid accumulation
#   to antioxidant capacity in Isatis indigotica fruit
#
# Input:
#   The input file corresponds to Supplementary Table S2
#   in the revised manuscript.
#
# Expected input file:
#   input/Table_S2_metabolites.csv
#
# Expected columns:
#   Index, Metabolite, S1_1, S1_2, S1_3, ..., S5_3,
#   Class I, Class II, Level, Q1 (Da), Q3 (Da),
#   Ionization model, Molecular Weight (Da), Formula, CAS
#
# Main package:
#   ropls
#
# R version used:
#   R 4.2.2
############################################################


############################
# 1. Load packages
############################

required_packages <- c(
  "ropls",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "tibble",
  "purrr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package not installed:", pkg))
  }
}

library(ropls)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)


############################
# 2. Create output folders
############################

dir.create("output", showWarnings = FALSE)
dir.create("output/OPLSDA_score_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("output/OPLSDA_VIP_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/OPLSDA_model_summaries", recursive = TRUE, showWarnings = FALSE)
dir.create("output/DAM_tables", recursive = TRUE, showWarnings = FALSE)


############################
# 3. Import metabolite table
############################

input_file <- "input/Table_S2_metabolites.csv"

met_raw <- read_csv(input_file, show_col_types = FALSE)

if (!"Index" %in% colnames(met_raw)) {
  stop("The input file must contain a column named 'Index'.")
}

sample_cols <- grep("^S[1-5]_[1-3]$", colnames(met_raw), value = TRUE)

if (length(sample_cols) != 15) {
  warning("Expected 15 sample columns named S1_1 to S5_3, but found ",
          length(sample_cols), ".")
  print(sample_cols)
}

annotation_cols <- setdiff(colnames(met_raw), sample_cols)

met_anno <- met_raw %>%
  select(all_of(annotation_cols)) %>%
  rename(Feature_ID = Index)

met_abundance <- met_raw %>%
  select(Index, all_of(sample_cols)) %>%
  rename(Feature_ID = Index)

met_abundance <- met_abundance %>%
  mutate(across(all_of(sample_cols), as.numeric))

feature_ids <- met_abundance$Feature_ID

feature_sample_mat <- met_abundance %>%
  select(all_of(sample_cols)) %>%
  as.data.frame()

rownames(feature_sample_mat) <- feature_ids

sample_feature_mat <- t(feature_sample_mat) %>%
  as.data.frame()

sample_feature_mat <- sample_feature_mat %>%
  mutate(across(everything(), as.numeric))

sample_info <- data.frame(
  Sample = rownames(sample_feature_mat),
  Stage = str_extract(rownames(sample_feature_mat), "^S[1-5]")
)

sample_info$Stage <- factor(
  sample_info$Stage,
  levels = c("S1", "S2", "S3", "S4", "S5")
)

if (any(is.na(sample_info$Stage))) {
  stop("Some sample names do not contain valid stage labels S1-S5.")
}

write_csv(sample_info, "output/sample_info.csv")


############################
# 4. Data preprocessing
############################

replace_na_half_min <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  min_positive <- suppressWarnings(min(x[x > 0], na.rm = TRUE))
  if (is.infinite(min_positive)) {
    min_positive <- 1
  }
  x[is.na(x)] <- min_positive / 2
  return(x)
}

sample_feature_imputed <- sample_feature_mat %>%
  mutate(across(everything(), replace_na_half_min))

sample_feature_log2 <- log2(sample_feature_imputed + 1)

zero_var_features <- apply(
  sample_feature_log2,
  2,
  function(x) sd(x, na.rm = TRUE) == 0
)

if (any(zero_var_features)) {
  message(sum(zero_var_features),
          " zero-variance features removed before OPLS-DA.")
  sample_feature_log2 <- sample_feature_log2[, !zero_var_features, drop = FALSE]
}


############################
# 5. Define pairwise comparisons
############################

stage_levels <- levels(sample_info$Stage)

pairwise_comparisons <- combn(stage_levels, 2, simplify = FALSE)


############################
# 6. Univariate statistics
############################

calculate_univariate_stats <- function(stage_a,
                                       stage_b,
                                       abundance_original,
                                       abundance_log2,
                                       sample_info,
                                       annotation_table) {
  
  samples_a <- sample_info %>%
    filter(Stage == stage_a) %>%
    pull(Sample)
  
  samples_b <- sample_info %>%
    filter(Stage == stage_b) %>%
    pull(Sample)
  
  common_features <- colnames(abundance_log2)
  
  mean_a <- colMeans(
    abundance_original[samples_a, common_features, drop = FALSE],
    na.rm = TRUE
  )
  
  mean_b <- colMeans(
    abundance_original[samples_b, common_features, drop = FALSE],
    na.rm = TRUE
  )
  
  pseudo <- 1
  
  FC <- (mean_b + pseudo) / (mean_a + pseudo)
  log2FC <- log2(FC)
  
  p_values <- sapply(common_features, function(fid) {
    x <- abundance_log2[samples_a, fid]
    y <- abundance_log2[samples_b, fid]
    
    if (sd(x, na.rm = TRUE) == 0 &&
        sd(y, na.rm = TRUE) == 0 &&
        mean(x, na.rm = TRUE) == mean(y, na.rm = TRUE)) {
      return(1)
    }
    
    p <- tryCatch(
      t.test(x, y, alternative = "two.sided")$p.value,
      error = function(e) NA_real_
    )
    
    return(p)
  })
  
  FDR <- p.adjust(p_values, method = "BH")
  
  stat_table <- data.frame(
    Feature_ID = common_features,
    Mean_stage_a = as.numeric(mean_a[common_features]),
    Mean_stage_b = as.numeric(mean_b[common_features]),
    FC = as.numeric(FC[common_features]),
    log2FC = as.numeric(log2FC[common_features]),
    P_value = as.numeric(p_values[common_features]),
    FDR = as.numeric(FDR[common_features]),
    Stage_a = stage_a,
    Stage_b = stage_b,
    Comparison = paste0(stage_a, "_vs_", stage_b)
  )
  
  stat_table <- stat_table %>%
    left_join(annotation_table, by = "Feature_ID")
  
  return(stat_table)
}


############################
# 7. Pairwise OPLS-DA function
############################

run_pairwise_oplsda <- function(stage_a,
                                stage_b,
                                data_log2,
                                data_original,
                                sample_info,
                                annotation_table) {
  
  comparison_name <- paste0(stage_a, "_vs_", stage_b)
  message("Running OPLS-DA: ", comparison_name)
  
  selected_samples <- sample_info %>%
    filter(Stage %in% c(stage_a, stage_b))
  
  selected_data <- data_log2[selected_samples$Sample, , drop = FALSE]
  group <- droplevels(selected_samples$Stage)
  
  opls_model <- opls(
    x = selected_data,
    y = group,
    predI = 1,
    orthoI = NA,
    permI = 200,
    scaleC = "pareto",
    fig.pdfC = "none"
  )
  
  model_summary <- as.data.frame(getSummaryDF(opls_model))
  model_summary$Comparison <- comparison_name
  
  write_csv(
    model_summary,
    paste0("output/OPLSDA_model_summaries/",
           comparison_name,
           "_model_summary.csv")
  )
  
  vip_values <- getVipVn(opls_model)
  
  vip_table <- data.frame(
    Feature_ID = names(vip_values),
    VIP = as.numeric(vip_values),
    Comparison = comparison_name
  ) %>%
    arrange(desc(VIP)) %>%
    left_join(annotation_table, by = "Feature_ID")
  
  write_csv(
    vip_table,
    paste0("output/OPLSDA_VIP_tables/",
           comparison_name,
           "_VIP.csv")
  )
  
  scoreMN <- getScoreMN(opls_model)
  score_df <- as.data.frame(scoreMN)
  score_df$Sample <- rownames(score_df)
  score_df$Stage <- selected_samples$Stage[
    match(score_df$Sample, selected_samples$Sample)
  ]
  
  component_names <- colnames(scoreMN)
  
  if (length(component_names) >= 2) {
    x_comp <- component_names[1]
    y_comp <- component_names[2]
  } else {
    x_comp <- component_names[1]
    y_comp <- component_names[1]
  }
  
  p <- ggplot(score_df, aes_string(x = x_comp, y = y_comp, color = "Stage")) +
    geom_point(size = 4, alpha = 0.9) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("OPLS-DA: ", stage_a, " vs ", stage_b),
      x = x_comp,
      y = y_comp,
      color = "Stage"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  ggsave(
    paste0("output/OPLSDA_score_plots/",
           comparison_name,
           "_score_plot.pdf"),
    plot = p,
    width = 6,
    height = 5
  )
  
  stat_table <- calculate_univariate_stats(
    stage_a = stage_a,
    stage_b = stage_b,
    abundance_original = data_original,
    abundance_log2 = data_log2,
    sample_info = sample_info,
    annotation_table = annotation_table
  )
  
  dam_table <- stat_table %>%
    left_join(
      vip_table %>% select(Feature_ID, VIP),
      by = "Feature_ID"
    ) %>%
    mutate(
      Regulation = case_when(
        log2FC >= 1 & FDR < 0.05 & VIP > 1 ~ "Up",
        log2FC <= -1 & FDR < 0.05 & VIP > 1 ~ "Down",
        TRUE ~ "Not_DAM"
      ),
      DAM = ifelse(Regulation %in% c("Up", "Down"), "Yes", "No")
    )
  
  write_csv(
    dam_table,
    paste0("output/DAM_tables/",
           comparison_name,
           "_DAM_statistics.csv")
  )
  
  return(dam_table)
}


############################
# 8. Run analyses
############################

all_dam_tables <- list()

for (cmp in pairwise_comparisons) {
  
  stage_a <- cmp[1]
  stage_b <- cmp[2]
  
  comparison_name <- paste0(stage_a, "_vs_", stage_b)
  
  dam_table <- run_pairwise_oplsda(
    stage_a = stage_a,
    stage_b = stage_b,
    data_log2 = sample_feature_log2,
    data_original = sample_feature_imputed,
    sample_info = sample_info,
    annotation_table = met_anno
  )
  
  all_dam_tables[[comparison_name]] <- dam_table
}


############################
# 9. Combine DAM results
############################

dam_all <- bind_rows(all_dam_tables)

write_csv(
  dam_all,
  "output/DAM_tables/all_pairwise_DAM_statistics.csv"
)

dam_summary <- dam_all %>%
  filter(DAM == "Yes") %>%
  group_by(Comparison, Regulation) %>%
  summarise(Number = n(), .groups = "drop") %>%
  tidyr::complete(
    Comparison,
    Regulation = c("Up", "Down"),
    fill = list(Number = 0)
  )

write_csv(
  dam_summary,
  "output/DAM_tables/all_pairwise_DAM_number_summary.csv"
)


############################
# 10. Save session information
############################

sink("output/sessionInfo.txt")
sessionInfo()
sink()

message("OPLS-DA and DAM analyses completed successfully.")
