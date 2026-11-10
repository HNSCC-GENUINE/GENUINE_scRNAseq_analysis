# Fig4-5 main panels: 7 CD8 T cell subsets (NK/MAIT/rdT used only in supplementary panels)
CD8_7 <- subset(CD8, subset = level2 %in% c("CD8_naive", "CD8_CM", "CD8_EM_GZMK+", "CD8_Tpex",
                                            "CD8_Tex", "CD8_ISG+", "Proliferating_CD8"))
CD8_7$level2 <- droplevels(CD8_7$level2)


# Fig 5. ============================================================================================================
## A. TCR counts per tissue ----
plot_df <- bind_rows(tcr$data) %>% filter(celltype %in% levels(CD8_7$level2)) %>%
  group_by(tissue, celltype) %>%
  dplyr::summarise(TCR_count = n(), .groups = "drop")

plot_df$celltype <- factor(plot_df$celltype, levels = levels(CD8_7$level2))

ggplot(plot_df, aes(x = tissue, y = TCR_count, fill = celltype)) +
  geom_col(position = "fill") +
  labs(x = "tissue", y = "Proportion of TCR count", fill = "celltype") +
  theme_classic() + scale_fill_manual(breaks = levels(CD8_7$level2), values = CD8@misc$cols$level2)


## B. CD8 TR prop change ----
tmp <- subset(CD8_7, Tissue != "post_LN") %>% subset(., level2_robust == "TR")
df <- tmp@meta.data %>%
  mutate(across(c(Patient, Response2, Tissue, level2), as.character)) %>%
  group_by(Patient, Response2, Tissue, level2) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(Patient, Response2, Tissue) %>%
  dplyr::mutate(sum = sum(n), ratio = n / sum) %>% ungroup() %>% dplyr::select(-c(n, sum)) %>%
  pivot_wider(names_from = "Tissue", values_from = "ratio") %>% as.data.frame() %>%
  filter(!Patient %in% c("P004", "P010", "P007")) %>%   
  replace(., is.na(.), 0) %>% dplyr::mutate(ratio = post_tumor / pre_tumor)
df$level2 <- factor(df$level2, levels = rev(levels(CD8_7$level2)))

ggplot(df, aes(level2, sqrt(ratio), fill = level2)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = FALSE, label.y.npc = 0.85) +
  coord_flip() +
  theme_classic() + scale_fill_manual(breaks = levels(CD8_7$level2), values = CD8@misc$cols$level2) +
  labs(x = "", fill = "celltype (n=8)")

ggplot(df, aes(level2, sqrt(ratio), fill = level2)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = FALSE, label.y.npc = 0.85) +
  coord_flip() +
  facet_grid(. ~ Response2) +
  theme_classic() + scale_fill_manual(breaks = levels(CD8_7$level2), values = CD8@misc$cols$level2) +
  labs(x = "", fill = "celltype (n=8)")


## C. Expanded TCR clones by timepoint and response ----
tcr_df <- bind_rows(tcr$data) %>% filter(tissue != "post_LN") %>%
  filter(celltype %in% levels(CD8_7$level2))

df <- tcr_df %>%
  filter(!CDR3.aa %in% c("None_None", "None_NA", "NA_None", "NA_NA")) %>%
  group_by(clinical_3, tissue, patient, CDR3.aa, celltype) %>%
  dplyr::summarise(cell = n(), .groups = "drop") %>%
  group_by(clinical_3, tissue, patient, CDR3.aa) %>%
  dplyr::mutate(clone_size = sum(cell)) %>%
  ungroup()

top5 <- c()
top10 <- c()
for (i in c("pre_tumor", "post_tumor")) {
  for (j in c("LR", "ER")) {
    df2 <- df %>%
      filter(tissue == i, clinical_3 == j) %>%
      distinct(patient, CDR3.aa, clone_size) %>%
      group_by(CDR3.aa) %>%
      dplyr::summarise(clone_size = sum(clone_size), .groups = "drop") %>%
      arrange(-clone_size)
    
    clones      <- df2$CDR3.aa
    top5_clone  <- clones[1:5]
    top10_clone <- clones[6:10]
    
    sub   <- tcr_df %>% filter(tissue == i, clinical_3 == j)
    top5  <- c(top5,  sub %>% filter(CDR3.aa %in% top5_clone)  %>% pull(cell_id))
    top10 <- c(top10, sub %>% filter(CDR3.aa %in% top10_clone) %>% pull(cell_id))
  }
}

CD8_7$top_clone <- case_when(colnames(CD8_7) %in% top5  ~ "Top5 expanded",
                             colnames(CD8_7) %in% top10 ~ "Top10 expanded",
                             TRUE ~ "Rest")

for (i in c("pre_tumor", "post_tumor")) {
  for (j in c("LR", "ER")) {
    tmp <- subset(CD8_7, subset = Tissue == i) %>% subset(., Response2 == j)
    p <- rasterize(DimPlot(tmp, reduction = "umap", raster = FALSE, pt.size = 0.5, group.by = "top_clone",
                           order = c("Top5 expanded", "Top10 expanded", "Rest"),
                           cols = c("Top5 expanded" = "darkred", "Top10 expanded" = "darkblue", "Rest" = "grey")), dpi = 500) +
      labs(title = paste(j, i), color = "", x = "umap_1", y = "umap_2") + NoAxes() + NoLegend()
    print(p)
  }
}

top_clones <- df %>%
  filter(tissue %in% c("pre_tumor", "post_tumor"), clinical_3 %in% c("LR", "ER")) %>%
  distinct(tissue, clinical_3, patient, CDR3.aa, clone_size) %>%
  group_by(tissue, clinical_3, CDR3.aa) %>%
  dplyr::summarise(clone_size = sum(clone_size), .groups = "drop") %>%
  group_by(tissue, clinical_3) %>%
  arrange(-clone_size, .by_group = TRUE) %>%
  dplyr::mutate(rank = row_number()) %>%
  filter(rank <= 10) %>%
  dplyr::mutate(rank_group = if_else(rank <= 5, "top1-5", "top6-10")) %>%
  ungroup()

plot_df_top10 <- tcr_df %>%
  inner_join(top_clones, by = c("tissue", "clinical_3", "CDR3.aa")) %>%
  dplyr::mutate(tissue_res = paste(tissue, clinical_3, sep = "-")) %>%
  group_by(tissue_res, celltype) %>%
  dplyr::summarise(cell = n(), .groups = "drop") %>%
  dplyr::mutate(celltype = factor(celltype, levels = levels(CD8_7$level2)))

ggplot(plot_df_top10, aes(tissue_res, cell, fill = celltype)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(breaks = levels(CD8_7$level2), values = CD8@misc$cols$level2) +
  labs(x = "", y = "Fraction of T cells", fill = "celltype",
       title = "Top 10 expanded clones") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## D. DEGs of pre-treatment CD8_EM_GZMK+ by response ----
tmp <- subset(CD8_7, subset = Tissue == "pre_tumor") %>% subset(., level2 == "CD8_EM_GZMK+")
Idents(tmp) <- "Response2"

deg <- FindMarkers(tmp, ident.1 = "ER", ident.2 = "LR", assay = "RNA")
df <- deg %>%
  mutate(gene = rownames(.),
         sig = case_when(p_val_adj <= 0.01 & avg_log2FC >  log2(1.3) ~ "ER",
                         p_val_adj <= 0.01 & avg_log2FC < -log2(1.3) ~ "LR",
                         TRUE ~ "NS"),
         ident1      = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "NS", TRUE ~ sig),
         gene_label1 = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "", sig == "NS" ~ "", TRUE ~ gene))

df2 <- df %>% filter(gene_label1 != "")

ggplot(df2, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(aes(col = ident1)) + geom_text_repel(aes(label = gene_label1), max.overlaps = 3000, size = 3) + theme_classic() +
  scale_color_manual(values = c(CD8@misc$cols$Response2, NS = "#ECECEC")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype = "dashed", color = "darkgrey") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 8)) + NoLegend() + xlim(c(-1.8, 2.7)) + ylim(c(0, 100))


## E. Associations related to GZMK+ CD8 T cells ----
df <- lunit %>% dplyr::select(Patient, Response2, tls_area_ca_ratio_post, tls_cluster_cnt_ca_ratio_post)

## (1) intra-TLS count ~ GZMK+/TR ratio
tmp <- CD8_7 %>% subset(., subset = Tissue == "post_tumor") %>% subset(., level2_robust == "TR")
df1 <- tmp@meta.data %>%
  mutate(across(c(Patient, Tissue, level2), as.character)) %>%
  group_by(Patient, Tissue, level2) %>% dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(Patient, Tissue) %>% dplyr::mutate(sum = sum(n), ratio = n / sum) %>% ungroup() %>%
  filter(!Patient %in% c("P010", "P007")) %>%         
  filter(level2 == "CD8_EM_GZMK+")
test <- merge(df, df1, by = "Patient")

## (2) CXCR5 ~ intra TLS area
tmp <- subset(CD8_7, Tissue == "post_tumor") %>% subset(., level2 == "CD8_EM_GZMK+")
score <- AddModuleScore(tmp, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Tumor_reactive1" = TR_gene1,
                                        "Tpex_score"      = pex_gene,
                                        "Exhaustion"      = exh_gene,
                                        "cytotoxicity"    = cyto_gene2,
                                        "CXCR5"           = "CXCR5"),
                        name = c("Stemness_Memory", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity", "CXCR5"))

df1 <- score@meta.data %>%
  mutate(across(c(Patient, level2), as.character)) %>%
  group_by(Patient, level2) %>%
  dplyr::summarise(CXCR5 = mean(CXCR56), .groups = "drop")

lun <- merge(df, df1, by = "Patient") %>% filter(!Patient %in% c("P010", "P007"))  

fit_ratio <- lm(ratio ~ tls_cluster_cnt_ca_ratio_post, data = test)
fit_cxcr5 <- lm(tls_area_ca_ratio_post ~ CXCR5,        data = lun)
.pf <- function(m) { f <- summary(m)$fstatistic; unname(pf(f[1], f[2], f[3], lower.tail = FALSE)) }
p_5E <- c(ratio = .pf(fit_ratio), CXCR5 = .pf(fit_cxcr5))
print(p_5E)

lab_5E <- function(m, key) sprintf("p = %s\nR2 = %.3f,  n = %d",
                                   format.pval(p_5E[[key]], digits = 2),
                                   summary(m)$r.squared, nobs(m))

## plot 1 : intra-TLS count ~ GZMK+/TR ratio
ggplot(test, aes(tls_cluster_cnt_ca_ratio_post, ratio)) +
  geom_point(aes(col = Response2)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = CD8@misc$cols$Response2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.2, color = "black",
           label = lab_5E(fit_ratio, "ratio")) +
  ylab("Tem_GZMK/TR ratio in post tumor") + NoLegend()

## plot 2 : CXCR5 ~ intra TLS area
ggplot(lun, aes(CXCR5, tls_area_ca_ratio_post)) +
  geom_point(aes(col = Response2)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = CD8@misc$cols$Response2) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.2, color = "black",
           label = lab_5E(fit_cxcr5, "CXCR5")) +
  NoLegend()


## F. Tumor-LN shared TCR count / circos ----
lv_cd8 <- levels(CD8_7$level2)
na_cdr3 <- c("NA_NA", "NA_None", "None_NA", "None_None")

tcr_df <- bind_rows(tcr$data) %>% filter(celltype %in% lv_cd8) %>%
  filter(tissue %in% c("post_LN", "post_tumor")) %>% distinct(.keep_all = TRUE)

overlap <- tcr_df %>% group_by(CDR3.aa, tissue) %>% dplyr::summarise(clone = sum(Clones), .groups = "drop") %>%
  filter(!CDR3.aa %in% na_cdr3) %>%
  pivot_wider(names_from = "tissue", values_from = "clone") %>%
  replace(., is.na(.), 0) %>% filter(post_LN > 0 & post_tumor > 0) %>% pull(CDR3.aa)

df <- tcr_df %>% group_by(CDR3.aa, tissue, celltype) %>% dplyr::summarise(total_clone = sum(Clones), .groups = "drop") %>%
  filter(!CDR3.aa %in% na_cdr3) %>% filter(CDR3.aa %in% overlap) %>% arrange(CDR3.aa) %>%
  dplyr::mutate(unique_clone = 1)

## bar plot
df3 <- df %>% group_by(tissue, celltype) %>% dplyr::summarise(uniq = n_distinct(CDR3.aa), .groups = "drop")
df3$celltype <- factor(df3$celltype, levels = rev(lv_cd8))
ggplot(df3, aes(tissue, uniq, fill = celltype)) + geom_col(position = "stack") +
  geom_text(position = position_stack(vjust = 0.5), aes(group = celltype, label = uniq)) + theme_linedraw() +
  scale_fill_manual(breaks = lv_cd8, values = CD8@misc$cols$level2) +
  ggtitle("unique overlapped clone")

## circos plot
matrix_list <- list()
for (cdr3 in unique(df$CDR3.aa)) {
  tmp <- df %>% filter(CDR3.aa == cdr3) %>% mutate(sample = paste0(celltype, "_", tissue)) %>% dplyr::select(sample, unique_clone)
  result <- t(tmp)
  colnames(result) <- result[1, ]
  matrix_list[[cdr3]] <- result %>% as.data.frame() %>% subset(rownames(result) == "unique_clone")
}

dff <- rbindlist(matrix_list, fill = TRUE) %>%
  replace(., is.na(.), 0) %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate_all(as.numeric) %>%
  dplyr::mutate(sum = rowSums(across(where(is.numeric)))) %>%
  mutate(across(starts_with("V"), ~ . / sum)) %>%
  dplyr::select(-sum)

patterns <- apply(dff, 2, function(col) paste(as.numeric(col != 0), collapse = ""))
grouped_cols <- lapply(unique(patterns), function(p) rowSums(dff[, which(patterns == p), drop = FALSE]))
names(grouped_cols) <- paste0("Group_", seq_along(grouped_cols))
new_dff <- as.data.frame(grouped_cols)

reorder_by_nonzero_pattern <- function(df, pattern_order) {
  pats <- sapply(colnames(df), function(col) paste(as.numeric(df[[col]] != 0), collapse = ""))
  sorted_col_names <- names(pats[order(pats, decreasing = TRUE)])
  ordered_cols <- unlist(lapply(pattern_order, function(p) grep(p, sorted_col_names, value = TRUE)))
  df[, c(ordered_cols, setdiff(sorted_col_names, ordered_cols)), drop = FALSE]
}
dff <- reorder_by_nonzero_pattern(new_dff, c("CM", "EM_GZMK", "Tex", "ISG+", "Proliferating"))

sectors_LN <- paste0(c("CD8_naive", "CD8_CM", "CD8_EM_GZMK+", "CD8_Tex", "CD8_ISG+"), "_post_LN")
sectors    <- c(sectors_LN, paste0(lv_cd8, "_post_tumor"))
stopifnot(all(sectors %in% rownames(dff)))     
dff <- dff[sectors, ]

pattern <- apply(as.data.frame(t(dff)), 1, function(x) paste(as.integer(x != 0), collapse = ""))
dff <- as.data.frame(t(dff)) %>%
  mutate(pattern = pattern) %>%
  group_by(pattern) %>% dplyr::summarise(across(everything(), sum)) %>%
  dplyr::select(-pattern) %>% as.data.frame() %>% t() %>% as.data.frame()

df$celltype <- factor(df$celltype, levels = lv_cd8)
sector_widths <- df %>% group_by(celltype, tissue) %>% dplyr::summarise(unique_clone = sum(unique_clone), .groups = "drop") %>%
  dplyr::mutate(ratio = unique_clone / sum(unique_clone)) %>%
  arrange(tissue, celltype) %>% pull(ratio)
stopifnot(length(sector_widths) == length(sectors))

colors    <- unname(CD8@misc$cols$level2[sub("_post_(LN|tumor)$", "", sectors)])
n_LN      <- length(sectors_LN)
gap.after <- c(rep(1, n_LN - 1), 10, rep(1, length(lv_cd8) - 1), 10)

circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.after = gap.after)
circos.initialize(factors = sectors, xlim = c(0, 1), sector.width = sector_widths / sum(sector_widths))

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim <- CELL_META$xlim; ylim <- CELL_META$ylim
  circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black", col = NA)
  circos.text(mean(xlim), mean(ylim), labels = CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, cex = 0.6)
}, bg.border = NA)

for (z in seq_len(ncol(dff))) {
  for (i in seq_len(n_LN)) {
    for (j in (n_LN + 1):length(sectors)) {
      if (dff[i, z] == 0 || dff[j, z] == 0) next
      lo_i <- if (z == 1) 0 else sum(dff[i, 1:(z - 1)])
      lo_j <- if (z == 1) 0 else sum(dff[j, 1:(z - 1)])
      circos.link(sector.index1 = sectors[i], point1 = c(lo_i, lo_i + dff[i, z]),
                  sector.index2 = sectors[j], point2 = c(lo_j, lo_j + dff[j, z]),
                  col = colors[j], border = "black")
    }
  }
}
circos.clear()


## G. intra tumor TCR overlap ----
### Intratumoral TCR overlap
## Note: the node layout was generated once and saved (see preset file)
lv_cd8  <- levels(CD8_7$level2)
na_cdr3 <- c("NA_NA", "NA_None", "None_NA", "None_None")

tcr_df <- bind_rows(tcr$data) %>% filter(celltype %in% lv_cd8) %>%
  filter(tissue == "post_tumor") %>% filter(clinical_3 == "LR")

overlap <- tcr_df %>% group_by(CDR3.aa, celltype) %>% dplyr::summarise(clone = sum(Clones), .groups = "drop") %>%
  filter(!CDR3.aa %in% na_cdr3) %>%
  pivot_wider(names_from = "celltype", values_from = "clone")
overlap <- overlap[rowSums(!is.na(overlap)) >= 3, ] %>% pull(CDR3.aa)

df <- tcr_df %>% group_by(CDR3.aa, celltype) %>%
  dplyr::summarise(total_clone = sum(Clones), unique_clone = 1, .groups = "drop") %>%
  filter(!CDR3.aa %in% na_cdr3) %>% filter(CDR3.aa %in% overlap) %>% arrange(CDR3.aa)

matrix_list <- list()
for (cdr3 in unique(df$CDR3.aa)) {
  tmp <- df %>% filter(CDR3.aa == cdr3) %>% dplyr::select(celltype, unique_clone)
  result <- t(tmp)
  colnames(result) <- result[1, ]
  matrix_list[[cdr3]] <- result %>% as.data.frame() %>% filter(rownames(.) == "unique_clone")
}
test <- rbindlist(matrix_list, fill = TRUE) %>%
  replace(., is.na(.), 0) %>% as.data.frame() %>% mutate_all(as.numeric)

vertex <- df %>% group_by(celltype) %>% dplyr::summarise(sum = n(), .groups = "drop")

df2 <- matrix(NA, nrow = length(lv_cd8), ncol = length(lv_cd8), dimnames = list(lv_cd8, lv_cd8))
for (i in seq_along(lv_cd8)) {
  for (j in seq_along(lv_cd8)) {
    if (i == j) next
    a <- lv_cd8[i]; b <- lv_cd8[j]
    df2[a, b] <- if (all(c(a, b) %in% colnames(test))) sum(test[[a]] == 1 & test[[b]] == 1) else 0
  }
}

dff <- expand.grid(from = lv_cd8, to = lv_cd8, stringsAsFactors = FALSE)
dff$importance <- df2[cbind(dff$from, dff$to)]

dff_unique <- dff %>%
  mutate(pair_min = pmin(from, to), pair_max = pmax(from, to)) %>%
  distinct(pair_min, pair_max, importance, .keep_all = TRUE) %>%
  dplyr::select(-pair_min, -pair_max) %>% filter(!is.na(importance))

g <- graph_from_data_frame(dff_unique, directed = FALSE)
E(g)$width  <- dff_unique$importance
V(g)$size   <- vertex$sum[match(V(g)$name, vertex$celltype)]
V(g)$color  <- CD8@misc$cols$level2[V(g)$name]

plot(g, edge.width = E(g)$width * 0.5,
     vertex.label.cex = 0.8,
     vertex.size = V(g)$size * 0.3,
     layout = layout_saved)


## H. Tumor-LN shared phenotype score ----
## Classify clones as LN-only, tumor-only, or shared 
lv_cd8 <- levels(CD8_7$level2)

clone_groups <- bind_rows(tcr$data) %>%
  filter(celltype %in% lv_cd8) %>%
  group_by(CDR3.aa) %>%
  dplyr::summarise(has_tumor = any(tissue == "post_tumor"),
                   has_ln    = any(tissue == "post_LN"), .groups = "drop") %>%
  mutate(clone_group = case_when(has_tumor & !has_ln ~ "Tumor_only",
                                 has_tumor &  has_ln ~ "Tumor_LN_share",
                                 !has_tumor & has_ln ~ "LN_only",
                                 TRUE                ~ "Other_or_missing")) %>%
  filter(!CDR3.aa %in% c("None_NA", "None_None", "NA_None"))
table(clone_groups$clone_group)

tcr_clone_map <- bind_rows(tcr$data) %>%
  filter(celltype %in% lv_cd8, tissue %in% c("post_tumor", "post_LN")) %>%
  distinct(.keep_all = TRUE) %>%
  left_join(clone_groups %>% dplyr::select(CDR3.aa, clone_group), by = "CDR3.aa") %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  dplyr::select(cell_id, clone_group)

CD8_7$clone_group <- tcr_clone_map$clone_group[match(colnames(CD8_7), tcr_clone_map$cell_id)]
CD8_7$clone_group[is.na(CD8_7$clone_group)] <- "NA"

tmp <- subset(CD8_7, subset = Tissue != "pre_tumor") %>% subset(., subset = clone_group != "NA") %>%
  subset(., subset = level2 == "CD8_EM_GZMK+")
tmp <- AddModuleScore(tmp, assay = "RNA",
                      features = list("Stemness_Memory" = Tscm.cm_gene,
                                      "Exhaustion"      = exh_gene,
                                      "cytotoxicity"    = cyto_gene2,
                                      "Tumor_reactive1" = TR_gene1),
                      name = c("Stemness_Memory", "Exhaustion", "cytotoxicity", "Tumor_reactive1"))
table(tmp$clone_group)

my_comparisons <- list(c("LN_only", "Tumor_LN_share"),
                       c("Tumor_LN_share", "Tumor_only"),
                       c("LN_only", "Tumor_only"))

md_5H <- tmp@meta.data %>%
  dplyr::select(clone_group, Stemness_Memory1, Exhaustion2, cytotoxicity3, Tumor_reactive14) %>%
  tidyr::pivot_longer(-clone_group, names_to = "score_name", values_to = "value")

stat_5H <- adj_list(md_5H, y = "value", x = "clone_group",
                    comparisons = my_comparisons, facet_by = "score_name",
                    method = "t.test", padj = "BH") %>% add_p_lab()
print(stat_5H %>% dplyr::select(score_name, group1, group2, p, p.adj, n_test, p.adj.sig), n = 20)

ggplot(md_5H, aes(clone_group, value)) +
  geom_violin() +
  facet_wrap(~ score_name, scales = "free_y", nrow = 2) +
  theme_classic() +
  ggpubr::stat_pvalue_manual(stat_5H, label = "p.adj.lab", size = 2.5, tip.length = 0.01) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred") + NoLegend()


# Supple Fig12. ================================================================================================
## A. GZMK+ GSEA in post tumor ----
tmp <- subset(CD8_7, subset = Tissue == "post_tumor")
Idents(tmp) <- "level2"

deg <- FindMarkers(tmp, ident.1 = "CD8_EM_GZMK+", assay = "RNA")  
a <- deg %>% mutate(gene = rownames(.)) %>% arrange(desc(avg_log2FC))   # top = GZMK+

ranks <- a %>% filter(!grepl("^MT-|RPS|RPL|KRT|IGL|IGK|IGH", gene)) %>%
  filter(p_val < 0.01) %>%
  dplyr::select(gene, avg_log2FC) %>% deframe()

set.seed(42)
fgseaRes_G <- fgsea(pathways = gdbr_list,
                    stats    = ranks,
                    eps = 0.0, nPermSimple = 10000,
                    minSize  = 3, maxSize  = 500,
                    nproc = 1)
fgseaResTidy_G <- fgseaRes_G %>% as_tibble() %>% arrange(desc(NES)) %>% filter(!is.na(NES))

print(fgseaResTidy_G %>% dplyr::filter(grepl(T_pat, pathway)) %>%
        dplyr::summarise(n_pathway   = dplyr::n(),
                         n_pval_0.01 = sum(pval <= 0.01),
                         n_padj_0.05 = sum(padj <= 0.05)))

df <- fgseaResTidy_G %>% filter(pval <= 0.01) %>% filter(grepl(T_pat, pathway)) %>%
  filter(!grepl("_THYMIC|REELIN|MONOCYTE|SYNAPSE|B_CELL|NEURO|NECRO", pathway)) %>%
  filter(NES > 0)

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x = "", y = "Normalized Enrichment Score\n(Rest <-> CD8_EM_GZMK+)",
       title = "GOBP pathways NES from GSEA \nin post-tumor CD8_EM_GZMK+", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")


## B. Tumor-LN TCR similarity ----
tmp <- tcr
tmp <- repFilter(tmp, .method = "by.meta", .query = list(celltype   = include(levels(CD8_7$level2))))
tmp <- repFilter(tmp, .method = "by.meta", .query = list(tissue     = include(c("post_LN", "post_tumor"))))
tmp <- repFilter(tmp, .method = "by.meta", .query = list(clinical_3 = include(c("LR", "ER"))))

df <- repOverlap(tmp$data, .method = "tversky", .verbose = FALSE)
ms <- tmp$meta %>% dplyr::select(Sample, patient, tissue, celltype)

df2 <- df %>% as.data.frame() %>% rownames_to_column(var = "source") %>%
  pivot_longer(cols = -source, names_to = "target", values_to = "Similarity") %>%
  filter(!is.na(Similarity)) %>%
  left_join(ms %>% rename_with(~ paste0("source_", .x), -Sample), by = c("source" = "Sample")) %>%
  left_join(ms %>% rename_with(~ paste0("target_", .x), -Sample), by = c("target" = "Sample")) %>%
  filter(source_patient == target_patient, source_tissue != target_tissue) %>%
  dplyr::mutate(cells = paste(pmin(source_celltype, target_celltype),
                              pmax(source_celltype, target_celltype), sep = "|"))

a <- df2 %>% group_by(cells) %>% dplyr::summarise(Mean = median(Similarity)) %>% arrange(-Mean) %>% pull(cells)
df2$cells <- factor(df2$cells, levels = a)

stat_12B <- adj_refall(df2, y = "Similarity", x = "cells", padj = "BH") %>% add_p_lab()
print(stat_12B %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig), n = 30)

ggplot(df2, aes(cells, Similarity)) +
  geom_hline(yintercept = mean(df2$Similarity), lty = "dashed", color = "darkgrey") +
  geom_boxplot() + geom_point() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_12B %>% dplyr::filter(p.adj.sig != "ns"),
                             x = "group2", label = "p.adj.sig", size = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "TCR similarity (tversky)")

