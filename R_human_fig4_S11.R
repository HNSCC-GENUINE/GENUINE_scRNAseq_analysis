# Fig4-5 main panels: 7 CD8 T cell subsets (NK/MAIT/rdT used only in supplementary panels)
CD8_7 <- subset(CD8, subset = level2 %in% c("CD8_naive", "CD8_CM", "CD8_EM_GZMK+", "CD8_Tpex",
                                            "CD8_Tex", "CD8_ISG+", "Proliferating_CD8"))
CD8_7$level2 <- droplevels(CD8_7$level2)

# Fig 4. ============================================================================================================================
## A. CD8 T UMAP ----
rasterise(DimPlot(CD8_7, reduction = "umap", raster = FALSE, group.by = "level2", pt.size = 0.3,
                  order = rev(levels(CD8_7$level2)), cols = CD8@misc$cols$level2), dpi = 500) +
  labs(title = "") & NoAxes()

## B. CD8 signature scores ----
tmp2 <- CD8_7 %>% subset(., subset = Tissue != "post_LN")
score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory"   = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1"   = TR_gene1, # Oliveira et al.
                                        "Tpex_score"        = pex_gene,
                                        "Exhaustion"        = exh_gene,
                                        "cytotoxicity"      = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>%
  mutate(across(c(Patient, Response2, level2, Tissue), as.character)) %>%
  group_by(Patient, Response2, level2, Tissue) %>%
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1),
                   Tumor_reactive13 = mean(Tumor_reactive13),
                   Exhaustion5      = mean(Exhaustion5),
                   cytotoxicity6    = mean(cytotoxicity6), .groups = "drop")

df2 <- df_final %>% pivot_longer(names_to = "func", values_to = "score",
                                 cols = c("Stemness_Memory1", "cytotoxicity6", "Tumor_reactive13", "Exhaustion5"))

df2$func   <- factor(df2$func, levels = c("Stemness_Memory1", "cytotoxicity6", "Tumor_reactive13", "Exhaustion5"))
df2$level2 <- factor(df2$level2, levels = rev(levels(CD8_7$level2)))

stat_4B <- adj_refall(as.data.frame(df2), y = "score", x = "level2",
                      facet_by = "func", method = "wilcox.test", padj = "BH") %>% add_p_lab()
print(stat_4B %>% dplyr::select(func, group2, p, p.adj, n_test, p.adj.sig), n = 100)

ggplot(df2, aes(level2, score, fill = level2)) +
  geom_boxplot() + geom_point(size = 1) +
  facet_wrap(func ~ ., scales = "free", nrow = 2) +
  theme_classic() +
  scale_fill_manual(values = CD8@misc$cols$level2) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_4B, x = "group2", label = "p.adj.lab", size = 2.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2,
               color = "darkred", fill = "darkred", position = position_dodge(.75)) +
  NoLegend()


## C. Top 10 TCR clones on UMAP ----
tcr_df <- bind_rows(tcr$data) %>% filter(tissue != "post_LN") %>%
  filter(celltype %in% levels(CD8_7$level2))

df <- tcr_df %>%
  group_by(clinical_3, tissue, CDR3.aa, celltype) %>%
  dplyr::summarise(clone = sum(Clones), .groups = "drop") %>%
  group_by(clinical_3, tissue, CDR3.aa) %>%
  dplyr::mutate(sum_clone = sum(clone)) %>% ungroup() %>%
  filter(!CDR3.aa %in% c("None_None", "None_NA", "NA_None", "NA_NA"))

top5 <- c()
top10 <- c()
for (i in c("pre_tumor", "post_tumor")) {
  for (j in c("ER", "LR")) {
    df2 <- df %>% filter(tissue == i) %>% filter(clinical_3 == j) %>% arrange(-sum_clone)
    clones <- df2 %>% pull(CDR3.aa) %>% unique()
    top5_clone  <- clones[1:5]
    top10_clone <- clones[6:10]
    
    sub <- tcr_df %>% filter(tissue == i) %>% filter(clinical_3 == j)
    top5  <- c(top5,  sub %>% filter(CDR3.aa %in% top5_clone)  %>% pull(cell_id))
    top10 <- c(top10, sub %>% filter(CDR3.aa %in% top10_clone) %>% pull(cell_id))
  }
}

CD8_7$top_clone <- case_when(colnames(CD8_7) %in% top5  ~ "Top5 expanded",
                             colnames(CD8_7) %in% top10 ~ "Top10 expanded",
                             TRUE ~ "Rest")

tmp <- subset(CD8_7, subset = Tissue != "post_LN")
rasterize(DimPlot(tmp, reduction = "umap", raster = FALSE, pt.size = 0.2, group.by = "top_clone",
                  order = c("Top5 expanded", "Top10 expanded", "Rest"),
                  cols = c("Top5 expanded" = "darkred", "Top10 expanded" = "darkblue", "Rest" = "grey")), dpi = 500) +
  labs(title = "Top clones in tissue", color = "", x = "umap_1", y = "umap_2") + NoAxes()


## D. TCR clonality ----
tcr_df <- bind_rows(tcr$data) %>% filter(tissue != "post_LN") %>% filter(celltype %in% levels(CD8_7$level2))

df <- tcr_df %>% group_by(celltype, CDR3.aa) %>% dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(celltype, n) %>% dplyr::summarise(count = sum(n), .groups = "drop") %>%
  dplyr::mutate(clonality = case_when(n == 1 ~ "1",
                                      n == 2 ~ "2",
                                      n > 2 ~ ">=3")) %>%
  group_by(celltype, clonality) %>% dplyr::summarise(count = sum(count), .groups = "drop")

df$clonality <- factor(df$clonality, levels = rev(c("1", "2", ">=3")))
df$celltype  <- factor(df$celltype, levels = levels(CD8_7$level2))

p1 <- ggplot(df, aes(count, clonality, fill = celltype)) + geom_col() +
  theme_classic() + scale_fill_manual(breaks = levels(CD8_7$level2), values = CD8@misc$cols$level2) +
  theme(axis.text = element_text(size = 10))

p1


## E. TR pre/post volcano plot ----
tmp <- subset(CD8_7, Tissue != "post_LN") %>% subset(., level2_robust == "TR")
Idents(tmp) <- "Tissue"

deg <- FindMarkers(tmp, ident.1 = "post_tumor", ident.2 = "pre_tumor", assay = "RNA")

df <- deg %>%
  mutate(gene = rownames(.),
         sig = case_when(p_val_adj <= 0.01 & avg_log2FC >  log2(1.3) ~ "post_tumor",
                         p_val_adj <= 0.01 & avg_log2FC < -log2(1.3) ~ "pre_tumor",
                         TRUE ~ "NS"),
         ident1       = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "NS", TRUE ~ sig),
         gene_label1  = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "", sig == "NS" ~ "", TRUE ~ gene),
         func         = gene2func(gene = gene_label1, celltype = "T"),
         gene_label_2 = case_when(gene_label1 != "" & func == "functional" ~ gene, TRUE ~ ""),
         ident2       = case_when(gene_label_2 == "" ~ "NS", TRUE ~ ident1))
df2 <- df %>% filter(gene_label1 != "")

ggplot(df2, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(aes(col = ident2)) + geom_text_repel(aes(label = gene_label_2), max.overlaps = 3000, size = 3) + theme_classic() +
  scale_color_manual(values = c("pre_tumor" = "black", "post_tumor" = unname(CD8@misc$cols$Tissue["post_tumor"]), NS = "#ECECEC")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype = "dashed", color = "darkgrey") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 8)) + NoLegend() + xlim(c(-1.8, 2.7)) + ylim(c(0, 100))


## F. TR pre/post GSEA GOBP ----
a <- deg %>% mutate(gene = rownames(.)) %>% arrange(desc(avg_log2FC))   # top = post

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
                         n_pval_0.05 = sum(pval <= 0.05),
                         n_padj_0.05 = sum(padj <= 0.05)))

df <- fgseaResTidy_G %>% filter(pval <= 0.05) %>% filter(grepl(T_pat, pathway)) %>%
  filter(!grepl("_THYMIC|REELIN|MONOCYTE|SYNAPSE|B_CELL|HUMORAL", pathway))

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x = "", y = "Normalized Enrichment Score\n(Pre <-> Post)",
       title = "GOBP pathways NES from GSEA \nin TR CD8 T cells", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")


## G. TR CD8 signature scores, pre vs post ----
tmp2 <- CD8_7 %>% subset(., subset = Tissue != "post_LN") %>% subset(., level2_robust == "TR")
score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory"   = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1"   = TR_gene1, # Oliveira et al.
                                        "Tpex_score"        = pex_gene,
                                        "Exhaustion"        = exh_gene,
                                        "cytotoxicity"      = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>%
  mutate(across(c(Patient, Response2, Tissue, level2), as.character)) %>%
  group_by(Patient, Response2, Tissue, level2) %>%
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1),
                   Tumor_reactive13 = mean(Tumor_reactive13),
                   Exhaustion5      = mean(Exhaustion5),
                   cytotoxicity6    = mean(cytotoxicity6), .groups = "drop")

df2 <- df_final %>%
  pivot_longer(names_to = "func", values_to = "score",
               cols = c("Stemness_Memory1", "Tumor_reactive13", "Exhaustion5", "cytotoxicity6")) %>%
  mutate(func2 = case_when(func == "Stemness_Memory1" ~ "Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                           func == "Tumor_reactive13" ~ "Tumor reactivity\n(Oliveira,Nat,2021)",
                           func == "Exhaustion5"      ~ "T cell Exhaustion\n(Beltra,Immunity,2020)",
                           func == "cytotoxicity6"    ~ "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)"))
df2$func2 <- factor(df2$func2, levels = c("Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                                          "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)",
                                          "Tumor reactivity\n(Oliveira,Nat,2021)",
                                          "T cell Exhaustion\n(Beltra,Immunity,2020)"))

df3 <- df2 %>%
  filter(!Patient %in% c("P010", "P007")) %>%      
  group_by(func2, Patient, level2) %>%
  filter(n_distinct(Tissue) == 2, n() == 2) %>%
  ungroup() %>%
  mutate(Tissue   = factor(Tissue, levels = c("pre_tumor", "post_tumor")),
         .pair_id = paste0(Patient, "_", level2))

stat_4G <- adj_two(as.data.frame(df3), y = "score", x = "Tissue",
                   facet_by = "func2", method = "wilcox.test",
                   paired = TRUE, id = ".pair_id", padj = "BH") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_4G %>% dplyr::select(func2, n, p, p.adj, n_test, p.adj.sig))   # n = number of pairs

ggplot(df3, aes(Tissue, score, fill = Tissue)) +
  geom_boxplot() + geom_point(size = 1) +
  facet_wrap(func2 ~ ., scales = "free", nrow = 2) +
  geom_line(aes(group = interaction(Patient, level2)), col = "darkgrey", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(values = CD8@misc$cols$Tissue) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_4G, label = "p.lab", size = 2.8, tip.length = 0.01) +
  stat_summary(aes(group = Tissue), fun = mean, geom = "point", shape = 18, size = 2, color = "darkred", fill = "darkred",
               position = position_dodge(.75)) + NoLegend()


## H. CD8 heatmap ----
marker <- c("CXCR3", "ITGAL", "ITGA4", "CXCR4", "CXCR5", "CCR4", "CCR6", "CXCL13", "CCL3", "CCL4", "CX3CR1", "IL21R", "CD40LG",
            "EOMES", "KLF2", "SELL", "S1PR5", "TCF7", "BCL6", "BHLHE40", "BACH2", "ID3", "BATF3", "LEF1", "CCR7", "IL7R", "FOXO1", "FOXP1", "S1PR1", "RUNX3",
            "SLAMF7", "FASLG", "CD27", "GZMK", "CD28", "CX3CR1", "IFNGR1", "KLRG1", "KLRB1", "KLF4", "XCL1", "BTLA", "IL2", "TNFSF14", "TNFRSF4", "TNFRSF9",
            "STAT3", "IL2RA", "MYC", "IRF8",
            "IFNG", "IL15RA", "TNFRSF18", "IL4", "LCK", "GZMB", "KLRC1", "KLRD1", "SAMD9", "TNFSF10", "ISG15", "MX1", "IFIT1", "IFIT3", "IFIT2", "STAT1",
            "STAT2", "IRF7", "TBX21", "TNF", "SOCS3", "STAT4", "IL6",
            "NFATC2", "TOX", "STAT5A", "CSF1", "TIGIT", "ID2", "ENTPD1", "IRF2", "LAG3", "PDCD1", "HAVCR2", "CD38", "PRDM1", "TOP2A", "MKI67")

r_split <- factor(c(rep("Tissue Homing", 11),
                    rep("B cell interaction", 2),
                    rep("Memory formation\nCirculation", 17),
                    rep("Activation\nEffector Differenciation", 43),
                    rep("Exhaustion", 13),
                    rep("Proliferation", 2)),
                  levels = c("Tissue Homing", "B cell interaction", "Memory formation\nCirculation",
                             "Activation\nEffector Differenciation", "Exhaustion", "Proliferation"))

tmp <- subset(CD8_7, Tissue != "post_LN") %>% subset(., level2_robust == "TR")
tmp$level2 <- as.character(tmp$level2)     
tmp$Tissue <- as.character(tmp$Tissue)

mat <- AverageExpression(tmp, assays = "RNA", features = marker,
                         group.by = c("level2", "Tissue"), verbose = FALSE)$RNA %>%
  as.matrix() %>% t() %>% scale() %>% t()

ct_order  <- c("CD8_EM_GZMK+", "CD8_Tpex", "CD8_Tex", "CD8_ISG+", "Proliferating_CD8")
col_order <- paste0(rep(ct_order, each = 2), "_", c("pre_tumor", "post_tumor"))
stopifnot(all(col_order %in% colnames(mat)))  
mat <- mat[, col_order]

keep    <- !duplicated(marker)
marker  <- marker[keep]
r_split <- r_split[keep]
r_split <- r_split[marker %in% rownames(mat)]
marker  <- marker[marker %in% rownames(mat)]
mat     <- mat[marker, ]

split <- factor(rep(ct_order, each = 2), levels = ct_order)

colHA <- HeatmapAnnotation(which = "column", show_annotation_name = FALSE, show_legend = TRUE,
                           annotation_name_gp = gpar(fontsize = 10), gp = gpar(col = "white"),
                           Celltype = split,
                           Tissue   = factor(rep(c("pre_tumor", "post_tumor"), length(ct_order)),
                                             levels = c("pre_tumor", "post_tumor")),
                           col = list(Celltype = CD8@misc$cols$level2[ct_order],
                                      Tissue   = CD8@misc$cols$Tissue[c("pre_tumor", "post_tumor")]))

p1 <- Heatmap(mat, name = "Expression",
              col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B0055FF", "#00588BFF", "#009B95FF", "#53CC67FF", "#FDE333FF")),
              show_column_dend = FALSE, show_row_dend = FALSE, cluster_columns = FALSE,
              cluster_row_slices = FALSE, cluster_column_slices = FALSE,
              row_order = marker, row_names_gp = gpar(fontsize = 10),
              row_split = r_split,
              row_title_rot = 90, row_title_gp = gpar(fontsize = 10), row_gap = unit(3, "mm"),
              column_order = col_order,
              column_title_side = "top",
              column_names_side = "top", column_names_rot = 90, column_names_centered = FALSE, column_title_gp = gpar(fontsize = 7),
              row_names_side = "left", show_heatmap_legend = TRUE,
              top_annotation = colHA,
              heatmap_height = unit(30, "cm"), heatmap_width = unit(10, "cm"),
              rect_gp = gpar(col = "black"), show_column_names = FALSE, show_row_names = TRUE)
print(p1)


## I. Association of intra-TLS score with Tex/TR ratio and exhaustion ----
## TLS - Tex/TR ratio correlation
df <- lunit %>% dplyr::select(Patient, Response2, tls_area_ca_ratio_post, tls_cluster_cnt_ca_ratio_post, IS_post)

tmp <- CD8_7 %>% subset(., subset = Tissue != "post_LN")
df1 <- tmp@meta.data %>% filter(level2_robust == "TR") %>%
  mutate(across(c(Patient, Tissue, level2), as.character)) %>%
  group_by(Patient, Tissue, level2) %>% dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(Patient, Tissue) %>% dplyr::mutate(sum = sum(n), ratio = n / sum) %>% ungroup() %>%
  filter(!Patient %in% c("P010", "P007")) %>%        
  filter(level2 == "CD8_Tex") %>% filter(Tissue == "post_tumor")

test <- merge(df, df1, by = "Patient")
print(lm_p(test, x = "tls_area_ca_ratio_post", y = "ratio"))

ggplot(test, aes(tls_area_ca_ratio_post, ratio)) + geom_point(aes(col = Response2)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = CD8@misc$cols$Response2) +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") +
  ylab("Tex/TR ratio in post tumor") + NoLegend()

## TLS - exhaustion score correlation
tmp <- subset(CD8_7, Tissue == "post_tumor") %>% subset(., level2 == "CD8_Tex")
score <- AddModuleScore(tmp, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Tumor_reactive1" = TR_gene1,
                                        "Tpex_score"      = pex_gene,
                                        "Exhaustion"      = exh_gene,
                                        "cytotoxicity"    = cyto_gene2),
                        name = c("Stemness_Memory", "Tumor_reactive1", "Tpex_score", "Exhaustion", "cytotoxicity"))

df1 <- score@meta.data %>%
  mutate(across(c(Patient, level2), as.character)) %>%
  filter(!Patient %in% c("P010", "P007")) %>%         
  group_by(Patient, level2) %>%
  dplyr::summarise(Stemness_Memory = mean(Stemness_Memory1),
                   Tumor_reactive1 = mean(Tumor_reactive12),
                   Tpex_score      = mean(Tpex_score3),
                   Exhaustion      = mean(Exhaustion4),
                   cytotoxicity    = mean(cytotoxicity5), .groups = "drop")
lun <- merge(df, df1, by = "Patient")

print(lm_p(lun, x = "tls_area_ca_ratio_post", y = "Exhaustion"))

ggplot(lun, aes(tls_area_ca_ratio_post, Exhaustion)) + geom_point(aes(col = Response2)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = CD8@misc$cols$Response2) +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") + NoLegend()


# Supple Fig11. ========================================================================================================================
## A. T/NK feature plot ----
list.rna <- rev(c("S1PR1", "CCR7", "TCF7", "BACH2", "SELL",
                  "CD28", "FOSB", "FAS",
                  "GZMK", "EOMES", "CD27", "PDCD1", "GZMM",
                  "XCL1", "STAT5A", "TNFRSF9",
                  "CXCL13", "TOX",
                  "ZNF683", "HOPX", "KLRD1", "KIR2DL4",
                  "IFIT1", "OAS1", "RSAD2",
                  "MKI67", "CDC25A"))

p1 <- FeaturePlot(CD8, features = rev(list.rna), min.cutoff = "q10", max.cutoff = "q90",
                  reduction = "umap", ncol = 8, raster = TRUE) & NoAxes()
rasterize(p1, dpi = 500)


## B. T/NK UMAP ----
rasterise(DimPlot(CD8, reduction = "umap", raster = FALSE, group.by = "level2", pt.size = 0.3,
                  order = rev(levels(CD8$level2)), cols = CD8@misc$cols$level2), dpi = 500) +
  labs(title = "") & NoAxes()


## C. CD8 dot plot ----
p1 <- DotPlot(CD8_7, features = list.rna, group.by = "level2", assay = "RNA", dot.scale = 5, col.max = 2) +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = levels(CD8_7$level2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        legend.position = "top",
        panel.grid = element_line(linewidth = 0.1, linetype = 2)) + NoLegend()
p1


## D. TR CD8 signature scores, pre vs post per cell type ----
tmp2 <- CD8_7 %>% subset(., subset = Tissue != "post_LN") %>% subset(., level2_robust == "TR")

score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory"   = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1"   = TR_gene1, # Oliveira et al.
                                        "Tpex_score"        = pex_gene,
                                        "Exhaustion"        = exh_gene,
                                        "cytotoxicity"      = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>%
  mutate(across(c(Patient, Response2, Tissue, level2), as.character)) %>%
  group_by(Patient, Response2, Tissue, level2) %>%
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1),
                   Tumor_reactive13 = mean(Tumor_reactive13),
                   Exhaustion5      = mean(Exhaustion5),
                   cytotoxicity6    = mean(cytotoxicity6), .groups = "drop")

df2 <- df_final %>%
  pivot_longer(names_to = "func", values_to = "score",
               cols = c("Stemness_Memory1", "Tumor_reactive13", "Exhaustion5", "cytotoxicity6")) %>%
  mutate(func2 = case_when(func == "Stemness_Memory1" ~ "Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                           func == "Tumor_reactive13" ~ "Tumor reactivity\n(Oliveira,Nat,2021)",
                           func == "Exhaustion5"      ~ "T cell Exhaustion\n(Beltra,Immunity,2020)",
                           func == "cytotoxicity6"    ~ "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)"))
df2$func2 <- factor(df2$func2, levels = c("Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                                          "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)",
                                          "Tumor reactivity\n(Oliveira,Nat,2021)",
                                          "T cell Exhaustion\n(Beltra,Immunity,2020)"))

ct_order  <- c("CD8_EM_GZMK+", "CD8_Tpex", "CD8_Tex", "CD8_ISG+", "Proliferating_CD8")
lv_sample <- paste0(rep(ct_order, each = 2), "_", c("pre_tumor", "post_tumor"))

df3 <- df2 %>%
  filter(!Patient %in% c("P010", "P007")) %>%        
  dplyr::mutate(sample = paste0(level2, "_", Tissue)) %>%
  group_by(func2, Patient, level2) %>%
  filter(n_distinct(Tissue) == 2, n() == 2) %>%
  ungroup() %>%
  dplyr::mutate(sample = factor(sample, levels = lv_sample),
                Tissue = factor(Tissue, levels = c("pre_tumor", "post_tumor")),
                level2 = factor(level2, levels = rev(levels(CD8_7$level2)))) %>%
  filter(!is.na(sample))

stat_paired <- df3 %>%
  dplyr::select(func2, level2, Patient, Tissue, score) %>%
  tidyr::pivot_wider(names_from = Tissue, values_from = score) %>%
  tidyr::drop_na(pre_tumor, post_tumor) %>%
  group_by(func2, level2) %>%
  dplyr::summarise(
    n_pair = dplyr::n(),
    p      = tryCatch(stats::wilcox.test(pre_tumor, post_tumor, paired = TRUE, exact = FALSE)$p.value,
                      error = function(e) NA_real_),
    .groups = "drop") %>%
  group_by(func2) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH"), n_test = dplyr::n()) %>%
  ungroup() %>%
  dplyr::mutate(group1    = paste0(level2, "_pre_tumor"),
                group2    = paste0(level2, "_post_tumor"),
                p.adj.sig = .sig(p.adj),
                label     = ifelse(is.na(p), "n.s.", paste0("p = ", signif(p, 2))))

print(stat_paired, n = 50)   # n_pair = number of patients per comparison

stat_paired <- df3 %>%
  group_by(func2) %>%
  dplyr::summarise(.lo = min(score, na.rm = TRUE),
                   .hi = max(score, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(.rng = .hi - .lo) %>%
  right_join(stat_paired, by = "func2") %>%
  group_by(func2) %>%
  dplyr::mutate(y.position = .hi + .rng * (0.06 + 0.08 * (dplyr::row_number() - 1))) %>%
  ungroup()

ggplot(df3, aes(sample, score, fill = Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1) +
  facet_wrap(func2 ~ ., scales = "free", nrow = 2) +
  geom_line(aes(group = interaction(Patient, level2)), col = "darkgrey", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(values = CD8@misc$cols$Tissue) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 7)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_paired, label = "label",
                             tip.length = 0.01, size = 2.5, bracket.size = 0.3) +
  stat_summary(aes(group = Tissue), fun = mean, geom = "point",
               shape = 18, size = 2, color = "darkred", fill = "darkred",
               position = position_dodge(.75)) +
  NoLegend()


## E. ISG score ----
tmp2 <- CD8_7 %>% subset(., subset = Tissue != "post_LN")
score <- AddModuleScore(tmp2, assay = "RNA", slot = "data",
                        features = list("good" = IFNG.GS,
                                        "bad"  = ISG.RS),
                        name = c("good", "bad"))

md <- score@meta.data %>%
  filter(level2_robust == "TR") %>%
  mutate(good_z = as.numeric(scale(good1)),
         bad_z  = as.numeric(scale(bad2)),
         delta  = good_z - bad_z,
         level2 = factor(as.character(level2), levels = intersect(levels(CD8_7$level2), unique(as.character(level2)))),
         Patient = as.character(Patient), Response2 = as.character(Response2))

df_final <- md %>%
  group_by(Patient, Tissue, level2, Response2) %>%
  dplyr::summarise(good_z = mean(good_z),
                   bad_z  = mean(bad_z),
                   ratio  = mean(delta),
                   n_cell = n(), .groups = "drop")

inter <- mean(df_final$ratio)
stat_11E1 <- adj_refall(df_final, y = "ratio", x = "level2", padj = "BH") %>% add_p_lab()
print(stat_11E1 %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

p1 <- ggplot(df_final, aes(level2, ratio, fill = level2)) +
  geom_boxplot() + ylab("delta good/bad IFNG effect") + xlab("") +
  theme_classic() + geom_point(position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = inter) +
  scale_fill_manual(values = CD8@misc$cols$level2) +
  ggpubr::stat_pvalue_manual(stat_11E1, x = "group2", label = "p.adj.lab", size = 2.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred") +
  NoLegend()

##
df_final$Tissue <- factor(as.character(df_final$Tissue), levels = c("pre_tumor", "post_tumor"))
stat_11E2 <- adj_two(as.data.frame(df_final), y = "ratio", x = "Tissue",
                     facet_by = "level2", method = "wilcox.test", padj = "BH") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_11E2 %>% dplyr::select(level2, n, p, p.adj, n_test, p.adj.sig))

p2 <- ggplot(df_final, aes(Tissue, ratio, fill = Tissue)) +
  geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) +
  facet_wrap(~ level2, nrow = 1) +
  theme_classic() + theme(legend.position = "top") +
  ggpubr::stat_pvalue_manual(stat_11E2, label = "p.lab", size = 2.2, tip.length = 0.01) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred") +
  scale_fill_manual(values = CD8@misc$cols$Tissue) + xlab("")

p1 / p2



