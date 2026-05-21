############################################################
# Script:
#   02_Fig5_heatmap_original_style_revised_boxed.R
#
# Purpose:
#   Generate a revised Fig. 5 heatmap using the original-style
#   group colors and group order, with enlarged fonts and 
#   legends positioned horizontally at the top, enclosed in 
#   dashed bounding boxes.
############################################################

############################
# 1. Load packages
############################

required_packages <- c(
  "readr",
  "dplyr",
  "ComplexHeatmap",
  "circlize",
  "grid"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package not installed: ", pkg))
  }
}

library(readr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)


############################
# 2. Input and output
############################

input_file <- "input/Fig5_flavonoid_genes_FPKM.tsv"

output_dir <- "output/Fig5_heatmap_original_style"

dir.create("output", showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################
# 3. Read input data
############################

dat <- read_tsv(input_file, show_col_types = FALSE)

required_cols <- c(
  "#ID",
  "S1-1_FPKM", "S1-2_FPKM", "S1-3_FPKM",
  "S2-1_FPKM", "S2-2_FPKM", "S2-3_FPKM",
  "S3-1_FPKM", "S3-2_FPKM", "S3-3_FPKM",
  "S4-1_FPKM", "S4-2_FPKM", "S4-3_FPKM",
  "S5-1_FPKM", "S5-2_FPKM", "S5-3_FPKM",
  "Gene_group"
)

missing_cols <- setdiff(required_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The following required columns are missing:\n",
      paste(missing_cols, collapse = ", ")
    )
  )
}

dat <- dat %>%
  rename(Gene_ID = `#ID`)

fpkm_cols <- c(
  "S1-1_FPKM", "S1-2_FPKM", "S1-3_FPKM",
  "S2-1_FPKM", "S2-2_FPKM", "S2-3_FPKM",
  "S3-1_FPKM", "S3-2_FPKM", "S3-3_FPKM",
  "S4-1_FPKM", "S4-2_FPKM", "S4-3_FPKM",
  "S5-1_FPKM", "S5-2_FPKM", "S5-3_FPKM"
)

dat <- dat %>%
  mutate(across(all_of(fpkm_cols), as.numeric)) %>%
  mutate(
    Gene_group = ifelse(
      is.na(Gene_group) | Gene_group == "" | Gene_group == "--",
      "Unclassified",
      Gene_group
    )
  )


############################
# 4. Original-style group order
############################

group_order <- c(
  "4CL", "ANR", "ANS", "bglB", "bglX", "C4H", "CHI", "CHS",
  "CYP81E", "DFR", "F3H", "FLS", "IF7MAT", "PAL", "PGT1",
  "PRX", "UGT73C6", "UGT75C1", "UGT78D2", "UGT79B1", "VR"
)

present_groups <- unique(dat$Gene_group)

group_order <- c(
  group_order[group_order %in% present_groups],
  setdiff(sort(present_groups), group_order)
)


############################
# 5. Original-style colors
############################

group_colors_base <- c(
  "4CL"     = "#294B52",  "ANR"     = "#E9923A",  "ANS"     = "#00A6D6",
  "bglB"    = "#A43D3F",  "bglX"    = "#79B89B",  "C4H"     = "#66629A",
  "CHI"     = "#7A7664",  "CHS"     = "#8A6A45",  "CYP81E"  = "#B7A99D",
  "DFR"     = "#9B7B91",  "F3H"     = "#005A9C",  "FLS"     = "#FF0000",
  "IF7MAT"  = "#35B64A",  "PAL"     = "#0EA5A6",  "PGT1"    = "#9B63A9",
  "PRX"     = "#FFA07A",  "UGT73C6" = "#B0002A",  "UGT75C1" = "#A8B4B5",
  "UGT78D2" = "#A23854",  "UGT79B1" = "#B58510",  "VR"      = "#294B52"
)

extra_groups <- setdiff(group_order, names(group_colors_base))

if (length(extra_groups) > 0) {
  extra_cols <- grDevices::rainbow(length(extra_groups))
  names(extra_cols) <- extra_groups
  group_colors_base <- c(group_colors_base, extra_cols)
}

group_colors_base <- group_colors_base[group_order]


############################
# 6. Calculate mean FPKM per stage
############################

stage_mean <- dat %>%
  mutate(
    S1 = rowMeans(select(., `S1-1_FPKM`, `S1-2_FPKM`, `S1-3_FPKM`), na.rm = TRUE),
    S2 = rowMeans(select(., `S2-1_FPKM`, `S2-2_FPKM`, `S2-3_FPKM`), na.rm = TRUE),
    S3 = rowMeans(select(., `S3-1_FPKM`, `S3-2_FPKM`, `S3-3_FPKM`), na.rm = TRUE),
    S4 = rowMeans(select(., `S4-1_FPKM`, `S4-2_FPKM`, `S4-3_FPKM`), na.rm = TRUE),
    S5 = rowMeans(select(., `S5-1_FPKM`, `S5-2_FPKM`, `S5-3_FPKM`), na.rm = TRUE)
  ) %>%
  select(Gene_ID, Gene_group, S1, S2, S3, S4, S5)

stage_mean <- stage_mean %>%
  mutate(
    Gene_group = factor(Gene_group, levels = group_order)
  ) %>%
  arrange(Gene_group, Gene_ID)

write_csv(
  stage_mean,
  file.path(output_dir, "Fig5_gene_stage_mean_FPKM.csv")
)


############################
# 7. Gene group counts for legend labels
############################

group_counts <- stage_mean %>%
  count(Gene_group, name = "n") %>%
  mutate(
    Gene_group = as.character(Gene_group),
    Group_label = paste0(Gene_group, " (n=", n, ")")
  ) %>%
  arrange(factor(Gene_group, levels = group_order))

write_csv(
  group_counts,
  file.path(output_dir, "Fig5_gene_group_counts.csv")
)

group_label_map <- setNames(group_counts$Group_label, group_counts$Gene_group)

stage_mean <- stage_mean %>%
  mutate(Gene_group_label = group_label_map[as.character(Gene_group)])

group_colors_with_counts <- group_colors_base[group_counts$Gene_group]
names(group_colors_with_counts) <- group_counts$Group_label


############################
# 8. Prepare expression matrix
############################

expr_mat <- stage_mean %>%
  select(S1, S2, S3, S4, S5) %>%
  as.data.frame()

rownames(expr_mat) <- stage_mean$Gene_ID
expr_log <- log2(expr_mat + 1)

zscore_by_gene <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

expr_z <- t(apply(expr_log, 1, zscore_by_gene))
expr_z <- as.matrix(expr_z)
heatmap_mat <- t(expr_z)
heatmap_mat <- heatmap_mat[, stage_mean$Gene_ID, drop = FALSE]


############################
# 9. Top annotation bar
############################

gene_group_label <- as.character(stage_mean$Gene_group_label)
names(gene_group_label) <- stage_mean$Gene_ID

top_anno <- HeatmapAnnotation(
  Group = gene_group_label,
  col = list(Group = group_colors_with_counts),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 12),
  show_annotation_name = TRUE,
  simple_anno_size = unit(4, "mm"),
  
  # 关闭默认图注，我们将在后面手动绘制它
  show_legend = FALSE 
)


############################
# 10. Column split for family separation
############################

column_split_labels <- factor(
  gene_group_label,
  levels = group_counts$Group_label
)


############################
# 11. Heatmap color scale
############################

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2E8B57", "#FFFFCC", "#B2182B")
)


############################
# 12. Draw Heatmap (Base)
############################

ht <- Heatmap(
  heatmap_mat,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = column_split_labels,
  cluster_column_slices = FALSE,
  column_gap = unit(1.2, "mm"),
  column_title = NULL,
  top_annotation = top_anno,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 15),
  rect_gp = gpar(col = NA),
  border = FALSE,
  
  # 关闭默认图注，我们将在后面手动绘制它
  show_heatmap_legend = FALSE 
)


############################
# 13. Create Custom Legends with Dashed Boxes
############################

# -- a. 构建 Z-score 图注
leg_zscore <- Legend(
  title = "Z-score",
  col_fun = col_fun,
  at = c(-2, -1, 0, 1, 2),
  labels = c("-2", "-1", "0", "1", "2"),
  direction = "horizontal",
  title_gp = gpar(fontsize = 12),
  labels_gp = gpar(fontsize = 10)
)

# -- b. 构建 Group 图注
leg_group <- Legend(
  title = "Group",
  labels = names(group_colors_with_counts),
  legend_gp = gpar(fill = group_colors_with_counts),
  nrow = 3,
  title_gp = gpar(fontsize = 12),
  labels_gp = gpar(fontsize = 10)
)

# -- c. 编写包装函数：给图注外层加上虚线方框
wrap_in_dashed_box <- function(leg) {
  # 计算原始图注尺寸并增加 4mm 的内边距 padding
  w <- grobWidth(leg) + unit(4, "mm")
  h <- grobHeight(leg) + unit(4, "mm")
  
  # 组合方框(rectGrob)和图注内容(leg)
  gb <- grobTree(
    rectGrob(width = w, height = h, gp = gpar(col = "black", lty = 2, fill = NA)),
    leg,
    vp = viewport(width = w, height = h)
  )
  
  # 赋予自定义类名，以便 R 能够自动获取其包装后的尺寸
  class(gb) <- c("dashed_legend", class(gb))
  return(gb)
}

# 注册提取自定义组件尺寸的 S3 方法
grobWidth.dashed_legend <- function(x) x$vp$width
grobHeight.dashed_legend <- function(x) x$vp$height

# -- d. 将图注放进虚线盒子
leg_zscore_boxed <- wrap_in_dashed_box(leg_zscore)
leg_group_boxed <- wrap_in_dashed_box(leg_group)

# -- e. 将打包好的图注水平排列在同一水平线上 (利用 column_gap 设置间距)
packed_legends <- packLegend(
  leg_zscore_boxed, 
  leg_group_boxed, 
  direction = "horizontal", 
  column_gap = unit(15, "mm") 
)


############################
# 14. Save heatmap
############################

n_genes <- ncol(heatmap_mat)
pdf_width <- max(11, min(18, n_genes * 0.09 + 6))
pdf_height <- 4.5 

png_width <- as.integer(pdf_width * 300)
png_height <- as.integer(pdf_height * 300)

pdf(
  file.path(output_dir, "Fig5B_flavonoid_gene_heatmap_original_style.pdf"),
  width = pdf_width,
  height = pdf_height
)

draw(
  ht,
  annotation_legend_list = packed_legends, # 注入自定义的组合图注
  annotation_legend_side = "top"           # 强制放置在最上方
)

dev.off()

png(
  file.path(output_dir, "Fig5B_flavonoid_gene_heatmap_original_style.png"),
  width = png_width,
  height = png_height,
  res = 300
)

draw(
  ht,
  annotation_legend_list = packed_legends, 
  annotation_legend_side = "top"
)

dev.off()


############################
# 15. Save session information
############################

sink(file.path(output_dir, "sessionInfo_Fig5_heatmap.txt"))
sessionInfo()
sink()

message("Fig. 5 heatmap with custom boxed legends generated successfully.")
