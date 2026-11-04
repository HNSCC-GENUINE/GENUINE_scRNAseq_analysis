# Fig3. ================================================================================================================================
## B. whole immune GSEA ----
tmp <- subset(adata_filtered, ident_tissue != "post_LN") %>% subset(., immune == "immune")
Idents(tmp) <- "ident_tissue"
deg <- FindMarkers(tmp, ident.1 = "post_tumor", ident.2 = "pre_tumor", assay = "RNA", logfc.threshold = 0.2, min.pct = 0.1)

a <- deg %>% mutate(gene = rownames(.))
a <- a[order(a$avg_log2FC, decreasing = T), ]

ranks <- a %>% filter(!grepl("^MT-|RPS|RPL|KRT|IGL|IGK|IGH", gene)) %>% 
  filter(p_val < 0.01) %>%
  as.data.frame() %>%
  dplyr::select(gene, avg_log2FC)  %>% deframe()

set.seed(42)
fgseaRes_G <- fgsea(pathways = gdbr_list , 
                    stats    = ranks,
                    eps = 0.0, nPermSimple = 10000,
                    minSize  = 3, maxSize  = 500)
fgseaResTidy_G <- fgseaRes_G %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy_G <- fgseaResTidy_G %>% filter(!is.na(.$NES))

df <- fgseaResTidy_G %>% filter(pval <= 0.05) %>% filter(grepl(immune_pat, pathway)) %>%
  filter(row_number() %in% c(1:10, (n() - 10):n()))

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score\n(Pre <-> Post)",
       title="GOBP pathways NES from GSEA \nin Whole tissue", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")


## C. B cell annotation, monocle3 ----
rasterize(DimPlot(B_subset, reduction = "standardintegrationharmony_umap", group.by = "celltype_B2", order = rev(levels_B2), cols = values_B2, 
                  size = 0.2), dpi = 500) + NoAxes() + labs(title = "") + NoLegend()

p1 <- plot_cells(cds.tissue, show_trajectory_graph = T, trajectory_graph_segment_size = 1,
                 label_principal_points = F, label_branch_points = F,  labels_per_group = T,
                 label_cell_groups = F, label_groups_by_cluster = F, label_leaves = F, label_roots = T,
                 rasterize = F, alpha = 0.5,
                 color_cells_by = "pathway", trajectory_graph_color = "black",
                 norm_method = "log", cell_size = 0.5, graph_label_size = 3) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15), legend.position = "top") +
  scale_color_manual(breaks = c("Path1", "Path2", "Mixed", "None"), values = c("darkred", "darkblue", "darkgreen", "grey"))+
  labs(x= "", y = "") + NoAxes() + NoLegend()

rasterise(p1, dpi = 500)


## D. EF/GC score feature plot ----
tmp <- subset(B_subset, subset = ident_tissue != "post_LN")

GC <- EFGC_ASC %>% filter(Group == "GC") %>% arrange(LogFC) %>% 
  filter(LogFC <= -1) %>% filter(P_adj <= 0.05) %>% pull(Gene)
EF <- EFGC_ASC %>% filter(Group == "EF") %>% arrange(-LogFC) %>%
  filter(LogFC >= 1) %>% filter(P_adj <= 0.05) %>% pull(Gene)

score <- AddModuleScore(tmp, assay = "RNA", 
                        features = list("GC_related" = GC,
                                        "EF_related" = EF,
                                        "rafi_TIL_specific" = TILB_gene),
                        name = c("GC_related", "EF_related", "rafi_TIL_specific"))

score$ident_tissue <- factor(score$ident_tissue, levels = levels_ident_tissue)
Minval <- min(c(score$GC_related1, score$EF_related2, score$rafi_TIL_specific3))

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GC_related2 = GC_related1 - Minval,
                EF_related3 = EF_related2 - Minval,
                rafi_TIL_specific4 = rafi_TIL_specific3 - Minval)

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GCEF_ratio = GC_related2/(EF_related3 + 0.001))

score$pseudotime <- score$pseudotime %>% unlist()

tmp2 <- subset(score, celltype_B2 %in% c("B_naive_1", "B_naive_2", "B_EGR1_activated", "B_NR4A2_activated",  
                                         "B_memory", "B_AtM", "B_GCLZ" , "B_GCDZ"))
rasterize(FeaturePlot(tmp2, features = c("GC_related1", "EF_related2"), reduction = "standardintegrationharmony_umap", 
                      max.cutoff = "q75", min.cutoff = "q25"), dpi = 500) & NoAxes()


## E. GC/EF ratio by path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))
df$log2_GCEF <- log2(df$GCEF_ratio)

stat_3E <- adj_two(df, y = "log2_GCEF", x = "pathway", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_3E %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(pathway, log2_GCEF, fill = pathway), col = "black") +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_3E, label = "p.lab", size = 3.2, tip.length = 0.01) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(breaks = c("Path1", "Path2"),
                    values = c("Path1" = "darkred", "Path2" = "darkblue")) +
  labs(y = "log2(GC/EF ratio)")


## F. HPV specificity along path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))

df_3F <- df %>% dplyr::mutate(log2_HPV = log2(rafi_TIL_specific4)) %>%
  dplyr::filter(is.finite(log2_HPV), is.finite(pseudotime))
stat_3F <- dplyr::bind_rows(
  lm_p(df_3F %>% dplyr::filter(pathway == "Path1"), x = "pseudotime", y = "log2_HPV") %>%
    dplyr::mutate(pathway = "Path1"),
  lm_p(df_3F %>% dplyr::filter(pathway == "Path2"), x = "pseudotime", y = "log2_HPV") %>%
    dplyr::mutate(pathway = "Path2"))
print(stat_3F %>% dplyr::select(pathway, n, p, r2))

ggplot(df, aes(pseudotime, log2(rafi_TIL_specific4), col = pathway)) + rasterize(geom_point(alpha = 0.2, size = 0.3), dpi = 500) + geom_smooth(method = "lm") +
  scale_color_manual(breaks = c('Path1', "Path2"), values = c("Path1" = "darkred", "Path2" = "darkblue")) +
  stat_poly_eq(aes(label = paste("p =", format.pval(after_stat(p.value), digits = 2), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, label.y = c(0.98, 0.88), size = 3) +
  theme_classic()


## G. cell prop along pseudotime ----
tmp <- subset(B_subset, subset = ident_tissue != "post_LN")

df <- tmp@meta.data %>% filter(pathway %in% c("Path1")) %>% 
  filter(!celltype_B2 %in% c("Plasma", "B_GCDZ"))
df$celltype_B2 <- factor(df$celltype_B2, levels = rev(levels_B2))

ggplot(df, aes(pseudotime, group = celltype_B2, fill = celltype_B2, col = celltype_B2)) +
  geom_density(position = "fill", alpha = 0.7) +
  facet_wrap(pathway~., scales = "free_x", ncol = 1) +
  theme_minimal() +
  scale_color_manual(breaks = levels_B2, values = values_B2) +
  scale_fill_manual(breaks = levels_B2, values = values_B2)  + NoLegend()


## H. HPV specificity by cell type ----
tmp <- subset(score, subset = pathway == 'Path1')

df <- tmp@meta.data %>% filter(celltype_B2 != 'B_GCDZ') %>% 
  group_by(patient_name, ident_tissue, celltype_B2) %>% dplyr::summarise(HPV = mean(rafi_TIL_specific4))

order <- df %>% group_by(celltype_B2) %>% dplyr::summarise(mean = mean(HPV)) %>% arrange(-mean) %>% pull(celltype_B2)
df$celltype_B2 <- factor(df$celltype_B2, levels = order)

stat_3H <- adj_refall(df, y = "HPV", x = "celltype_B2", padj = "BH") %>% add_p_lab()
print(stat_3H %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(celltype_B2, HPV, fill = celltype_B2)) +
  geom_hline(yintercept = mean(df$HPV), lty = "dashed") +
  geom_boxplot() + geom_point() +
  scale_fill_manual(breaks = order, values = values_B2) + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_3H, x = "group2", label = "p.adj.lab", size = 2.5)


## I. SHM of B-AtM ----
df <- vdj %>% filter(ident_tissue != "post_LN") %>% filter(pathway == "Path1") %>% filter(celltype_B2 == "B_AtM")

stat_3I <- adj_two(df, y = "mu_freq", x = "ident_tissue", facet_by = "isotype_switch",
                   method = "t.test", padj = "BH") %>% add_p_lab(prefix = "t.test p = ")
print(stat_3I %>% dplyr::select(isotype_switch, group1, group2, n, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(ident_tissue, mu_freq, fill = ident_tissue)) + geom_boxplot() +
  facet_grid(~isotype_switch) + theme_bw() +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue) +
  ggpubr::stat_pvalue_manual(stat_3I, label = "p.lab", size = 3, tip.length = 0.01)


## K. intra-TLS score ----
df <- lunit %>% dplyr::select(patient_name, ERLR, tls_area_ca_ratio_post, tls_area_ca_ratio_pre, `Pathologic Tumor regression`)

df2 <- df %>% pivot_longer(cols = c("tls_area_ca_ratio_post", "tls_area_ca_ratio_pre"), names_to = "tissue", values_to = "tls_score") %>% 
  dplyr::mutate(tissue = case_when(tissue == "tls_area_ca_ratio_post" ~ "post", T ~ "pre")) %>% 
  filter(!is.na(tls_score)) %>% filter(!patient_name %in% c("P07", "P06"))

df2$tissue <- factor(df2$tissue, levels = c("pre", "post"))

stat_3K <- adj_two(df2, y = "tls_score", x = "tissue",
                   method = "wilcox.test", paired = TRUE, id = "patient_name") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_3K %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df2, aes(tissue, tls_score)) + geom_boxplot() +
  ggpubr::stat_pvalue_manual(stat_3K, label = "p.lab", size = 3.2, tip.length = 0.01) +
  theme_classic() + geom_point()


## L. association intra-TLS and p-tumor regression ----
df2 <- df %>% filter(!patient_name %in% c("P07", "P08"))

stat_3L <- lm_p(df2, x = "tls_area_ca_ratio_post", y = "Pathologic Tumor regression")
print(stat_3L)

ggplot(df2, aes(tls_area_ca_ratio_post, `Pathologic Tumor regression`)) + geom_point(aes(col = ERLR)) +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = values_response_3) +
  theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 8), ",\n",after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") + NoLegend()


# Supple Fig7.================================================================================================================================
## A. whole cell annotation ----
DimPlot(adata_filtered, reduction = "standardintegrationharmony_umap",  raster = F, group.by = "major_anno_yj_4", pt.size = 0.2, label = F, 
        order = rev(levels_major),cols = values_major) + NoAxes() +labs(title = "")

## B. level1 dot plot ----
list.rna <- c("CPA3",  "KIT", "FCGR3A", 'CD14', "CD68", "APOE", "C1QC", "LAMP3", "CD1C", 
              "TRDC", "NCAM1", "CD3E", "CD8A", "PDCD1", "CD4", "FOXP3", "IL2RA", "CD19", "MS4A1", 
              "MZB1", "IGHG1", "IGHA1", "EPCAM", "MUC1", "PECAM1", "VWF", "COL1A1", "MKI67", "TOP2A")

DotPlot(adata_filtered, features = list.rna, group.by = "major_anno_yj_4", dot.scale = 7, assay = "RNA") +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = rev(levels_major_yj2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_line(size = 0.1, linetype = 2)) + NoLegend()

## C. level1 feature plot ----
rasterize(FeaturePlot(adata_filtered, reduction = "standardintegrationharmony_umap",  raster = T,  
                      min.cutoff = "q5", max.cutoff = "q95",
                      features = list.rna, ncol = 7), dpi = 500) & NoAxes() +
  theme(text = element_text(size = 5),
        legend.key.height = unit(0.2, 'cm'), 
        legend.key.width = unit(0.1, 'cm'))

## D. level1 prop change ----
tmp <- subset(adata_filtered, ident_tissue != "post_LN") %>% subset(., immune == "immune")
df <- tmp@meta.data %>% group_by(patient_name, ident.clinical_response_3, ident_tissue, major_anno_yj_4) %>% 
  dplyr::summarise(n = n()) %>% ungroup() %>% group_by(patient_name, ident.clinical_response_3, ident_tissue) %>% 
  dplyr::mutate(sum = sum(n), ratio = n/sum) %>% ungroup() %>% dplyr::select(-c(n, sum)) %>% 
  pivot_wider(names_from = "ident_tissue", values_from = "ratio") %>% as.data.frame() %>% 
  filter(!patient_name %in% c("P11", "P07", "P08")) %>% replace(., is.na(.), 0) %>% dplyr::mutate(ratio = post_tumor/pre_tumor) %>% 
  filter(ident.clinical_response_3 %in% c("ER", "LR"))
df$major_anno_yj_4 <- factor(df$major_anno_yj_4, levels = rev(levels_major))

ggplot(df, aes(major_anno_yj_4, sqrt(ratio), fill = major_anno_yj_4)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = F, label.y.npc = 0.85) +
  coord_flip() +
  theme_classic() + scale_fill_manual(breaks = levels_major, values = values_major) +
  labs(x = "", fill = "celltype (n=8)")

ggplot(df, aes(major_anno_yj_4, sqrt(ratio), fill = major_anno_yj_4)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = F, label.y.npc = 0.85) +
  coord_flip() +
  facet_grid(.~ident.clinical_response_3) +
  theme_classic() + scale_fill_manual(breaks = levels_major, values = values_major) +
  labs(x = "", fill = "celltype (n=8)")

# Supple Fig8. ================================================================================================================================
## A. B feature plot ----
list.rna <- c("IGHD", "FCER2", "TCL1A", "YBX3",  "ZNF331", "CCR7", "KLF2",
              "HSPA1A","EGR1", "NR4A2", 
              "CD69", "NFKB1",  "VPS37B", 
              "CRIP1",  "ITGB1", 
              "FCRL5", "IFI44L", "TBX21", "ITGAX", "TLR7",  "STAT1",  "ZEB2",
              "DUSP4", "SYK", "HHEX", "MYC", "FOXO1", "SOX4",
              "CXCL13", "CXCR5", "CXCR3", "CCL4", "CCL5", "CD40",  "IL21R",  "AICDA",  "MIR155HG","EBI3", "TRAF1", 
              "BCL6", "IRF4", "NEIL1",    "SUGCT", "MEF2B", "STMN1", "MKI67",
              "CD27",  "XBP1", "PRDM1", 
              "IGHG1", "IGHA1")

DefaultAssay(B_subset) <- "RNA"
rasterize(FeaturePlot(B_subset, reduction = "standardintegrationharmony_umap",  raster = T, min.cutoff = "q5", max.cutoff = "q95", pt.size = 1.5,
                      features = list.rna, ncol = 6), dpi = 500) & NoAxes() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        legend.key.height = unit(0.2, 'cm'), 
        legend.key.width = unit(0.1, 'cm'))

## B. B dot plot ----
DotPlot(B_subset, features = list.rna, group.by = "celltype_B2", dot.scale = 5) +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = rev(levels_B2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_line(size = 0.1, linetype = 2)) + NoLegend()

# Supple Fig9. ================================================================================================================================
## A. HPV specificity per path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))
df$log2_HPV <- log2(df$rafi_TIL_specific4)

stat_9A <- adj_two(df, y = "log2_HPV", x = "pathway", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9A %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(pathway, log2_HPV, fill = pathway), col = "black") +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_9A, label = "p.lab", size = 3.2, tip.length = 0.01) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(breaks = c("Path1", "Path2"),
                    values = c("Path1" = "darkred", "Path2" = "darkblue")) +
  labs(y = "log2(TIL-B specificity score)")


## B. association GC - HPV or SHM ----
tmp <- subset(B_subset, subset = ident_tissue != "post_LN") %>% subset(., pathway == "Path1")

GC <- EFGC_ASC %>% filter(Group == "GC") %>% arrange(LogFC) %>% 
  filter(LogFC <= -1) %>%
  filter(P_adj <= 0.05) %>% pull(Gene)
EF <- EFGC_ASC %>% filter(Group == "EF") %>% arrange(-LogFC) %>%
  filter(LogFC >= 1) %>%
  filter(P_adj <= 0.05) %>% pull(Gene)

score <- AddModuleScore(tmp, assay = "RNA", 
                        features = list("GC_related" = GC,
                                        "EF_related" = EF,
                                        "rafi_TIL_specific" = TILB_gene),
                        name = c("GC_related", "EF_related", "rafi_TIL_specific"))


score$ident_tissue <- factor(score$ident_tissue, levels = levels_ident_tissue)
score$celltype_B <- factor(score$celltype_B, levels = levels_B)

Minval <- min(c(score$GC_related1, score$EF_related2, score$rafi_TIL_specific3))

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GC_related2 = GC_related1 - Minval,
                EF_related3 = EF_related2 - Minval,
                rafi_TIL_specific4 = rafi_TIL_specific3 - Minval)

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GCEF_ratio = GC_related2/(EF_related3 + 0.001))
score$pseudotime <- score$pseudotime %>% unlist()
meta <- score@meta.data %>%
  dplyr::select(barcode4merge2, GC_related2, EF_related3, rafi_TIL_specific4, GCEF_ratio)


## GC - HPV
meta2 <- vdj %>% filter(ident_tissue != "post_LN") %>% 
  dplyr::select(c(barcode4merge2, patient_name, celltype_B2, ident_tissue, ident.clinical_response_3, mu_freq))%>% left_join(., meta, by = "barcode4merge2")

df <- meta2 %>% filter(!is.na(GC_related2)) %>% filter(celltype_B2 != "B_GCDZ") 
df <- df %>% group_by(patient_name, celltype_B2) %>%
  dplyr::summarise(GC = mean(GC_related2), Mu = mean(mu_freq), HPV = mean(rafi_TIL_specific4), GCEF = mean(GCEF_ratio))

df_9B1  <- df
stat_9B1 <- lm_p(df_9B1, x = "GC", y = "HPV")
print(stat_9B1)

ggplot(df, aes(GC, HPV)) + geom_point() + geom_smooth(method = "lm", color = "darkred") + theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n",after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black")

## GC - Mu
meta2 <- vdj %>% filter(ident_tissue != "post_LN") %>% 
  dplyr::select(c(barcode4merge2, patient_name, celltype_B2, ident_tissue, ident.clinical_response_3, mu_freq, isotype_switch))%>% left_join(., meta, by = "barcode4merge2")

df <- meta2 %>% filter(!is.na(GC_related2)) %>% filter(celltype_B2 != "B_GCDZ") 
df <- df %>% group_by(patient_name, celltype_B2, isotype_switch) %>%
  dplyr::summarise(GC = mean(GC_related2), Mu = mean(mu_freq), HPV = mean(rafi_TIL_specific4))

df_9B2  <- df
stat_9B2 <- lm_p(df_9B2, x = "GC", y = "Mu")
print(stat_9B2)

ggplot(df, aes(GC, Mu)) + geom_point() + geom_smooth(method = "lm", color = "darkred") + theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n",after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black")

## C. GC/EF ratio by response ----
df <- score@meta.data %>% filter(ident_tissue != "post_LN") %>% filter(pathway == "Path1") %>% filter(celltype_B2 == "B_AtM")
df$ident_tissue <- factor(as.character(df$ident_tissue), levels = c("pre_tumor", "post_tumor"))

stat_9C <- adj_two(df, y = "GCEF_ratio", x = "ident_tissue",
                   facet_by = "ident.clinical_response_3", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9C %>% dplyr::select(ident.clinical_response_3, n, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(ident_tissue, GCEF_ratio, fill = ident_tissue)) + geom_violin() + geom_boxplot() +
  facet_grid(.~ident.clinical_response_3) + theme_classic() +
  stat_summary(fun = "mean", geom = "point", shape = 18, color = "darkred", size = 3) +
  ggpubr::stat_pvalue_manual(stat_9C, label = "p.lab", size = 3, tip.length = 0.01) +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue)


## D. path2 cell prop by pseudotime ----
tmp <- subset(B_subset, subset = ident_tissue != "post_LN")

df <- tmp@meta.data %>% filter(pathway %in% c("Path2")) %>% 
  filter(!celltype_B2 %in% c("Plasma", "B_GCDZ"))
df$celltype_B2 <- factor(df$celltype_B2, levels = rev(levels_B2))

ggplot(df, aes(pseudotime, group = celltype_B2, fill = celltype_B2, col = celltype_B2)) +
  geom_density(position = "fill", alpha = 0.7) +
  facet_wrap(pathway~., scales = "free_x", ncol = 1) +
  theme_minimal() +
  scale_color_manual(breaks = levels_B2, values = values_B2) +
  scale_fill_manual(breaks = levels_B2, values = values_B2)  + NoLegend()

## E. GC score in path1 by celltype ----
df <- score@meta.data %>% filter(pathway %in% c("Path1")) %>% filter(celltype_B2 != "B_GCDZ")
df1 <- df %>% group_by(patient_name, ident_tissue, celltype_B2) %>% dplyr::summarise(Mean = mean(GC_related2))
order <- df1 %>% group_by(celltype_B2) %>% dplyr::summarise(mean = mean(Mean)) %>% arrange(-mean) %>% pull(celltype_B2) 

stat_9E <- adj_refall(df1, y = "Mean", x = "celltype_B2", padj = "BH") %>% add_p_lab()
print(stat_9E %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

ggplot(df1, aes(celltype_B2, Mean, fill = celltype_B2), col = "black") +
  geom_boxplot() + geom_point() +
  ggpubr::stat_pvalue_manual(stat_9E, x = "group2", label = "p.adj.lab", size = 2.5) +
  theme_classic() +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(values = values_B2) + NoLegend()


## F. BCR diversity in path1 per celltype ----
row_name <- c("UnSw_pre_tumor", "UnSw_post_tumor", "Sw_pre_tumor", "Sw_post_tumor")
col_name <- levels_B2

vdj_tissue <- vdj %>% filter(ident_tissue != "post_LN") %>% mutate(iso_tissue = paste0(isotype_switch, "_", ident_tissue))
tmp2 <- vdj_tissue %>% filter(pathway == "Path1")  %>% filter(celltype_B2 != "B_GCDZ") 

mat_list <- list()
for (k in unique(tmp2$patient_id)) {
  vdj <- tmp2
  tmp <- vdj %>% filter(patient_id == k)
  val <- vdj %>% filter(patient_name == unique(tmp$patient_name)) %>% pull(ident_tissue) %>% unique() %>% length()
  print(val)
  
  if (val > 0) {
    print(unique(tmp$patient_id))
    mat <- matrix(nrow = length(row_name), ncol = length(col_name), data = 0)
    
    for (i in 1:length(row_name)) {
      for (j in 1:length(col_name)){
        
        clone_counts <- tmp %>% filter(iso_tissue == row_name[[i]]) %>% filter(celltype_B2 == col_name[[j]]) %>% 
          group_by(clone_id.x) %>%summarise(n= n()) %>% pull(n)
        
        shannon <- diversity(clone_counts, index = "invsimpson") 
        mat[i, j] <- shannon}}
    
    rownames(mat) <- row_name
    colnames(mat) <- col_name
    
    mat <- as.data.frame(mat) %>% dplyr::mutate(patient = k, iso_tissue = rownames(.))
    mat <- pivot_longer(mat, cols = levels_B2, names_to = "celltype", values_to = "div") %>% as.data.frame()
    
    mat_list[[k]] <- mat
  } else {next}
}

patient_key = c(c("P01", "P02", "P03", "P04"), c("P05", "P06", "P07"), c("P08", "P09", "P10", "P11"))

div_LR <- mat_list[grep("P01|P02|P03|P04", names(mat_list))] 
names(div_LR) <- NULL
div_LR <- do.call("rbind", div_LR) %>% mutate(response = "LR")

div_ER <- mat_list[grep("P05|P06|P07|P08|P09|P10|P11", names(mat_list))]
names(div_ER) <- NULL
div_ER <- do.call("rbind", div_ER) %>% mutate(response = "ER")

div <- do.call("rbind", list(div_ER, div_LR)) %>% filter(div != Inf)
div <- div %>% dplyr::mutate(patient_name = str_extract(patient, paste(patient_key, collapse = "|")),
                             isotype_switch = str_extract(iso_tissue, "UnSw|Sw"),
                             ident_tissue = str_extract(iso_tissue, "pre_tumor|post_tumor"))

div$response <- factor(div$response, levels = levels_response_3)
div$celltype <- factor(div$celltype, levels = levels_B2)
div$ident_tissue <- factor(div$ident_tissue, levels = levels_ident_tissue)

sort <- div %>% group_by(celltype) %>% dplyr::summarise(Mean = mean(div)) %>% dplyr::arrange(-Mean) %>% pull(celltype)
div$celltype <- factor(div$celltype, levels = sort)


div$div_sqrt <- sqrt(div$div)
stat_9F <- adj_refall(div, y = "div_sqrt", x = "celltype", padj = "BH") %>% add_p_lab()
print(stat_9F %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

ggplot(div, aes(celltype, div_sqrt, fill = celltype)) +
  geom_boxplot() + geom_point() +
  theme_classic() + ylab("sqrt(inv.simpson)") + xlab("") +
  geom_hline(yintercept = sqrt(mean(div$div)), lty = "dashed") +
  ggpubr::stat_pvalue_manual(stat_9F, x = "group2", label = "p.adj.lab", size = 2.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred", position = position_dodge(.75)) +
  scale_fill_manual(breaks = levels_B, values = values_B) + NoLegend() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


## G. GC score of B-AtM per timepoint ----
df <- score@meta.data %>% filter(pathway == "Path1") %>% filter(celltype_B2 == "B_AtM")
df$ident_tissue <- factor(as.character(df$ident_tissue), levels = c("pre_tumor", "post_tumor"))

stat_9G <- adj_two(df, y = "GC_related2", x = "ident_tissue", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9G %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(ident_tissue, GC_related2, fill = ident_tissue)) +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_9G, label = "p.lab", size = 3.2, tip.length = 0.01) +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred", position = position_dodge(.75))

# Fig 4. ============================================================================================================================
## A. T UMAP ----
rasterise(DimPlot(CD8, reduction = "standardintegrationharmony_umap",  raster = F, group.by = "celltype_T", pt.size = 0.3, 
                  order = rev(levels_CD8), cols = rev(values_CD8)), dpi = 500) + 
  labs(title = "") & NoAxes()

## B. CD8 signature scores ----
tmp2 <- CD8 %>% subset(., subset = ident_tissue != "post_LN")
score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1" = TR_gene1, # oliveira
                                        "Tpex_score" = pex_gene,
                                        "Exhaustion" = exh_gene,
                                        "cytotoxicity" = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>% group_by(patient_ID, patient_name, ident.clinical_response_3, celltype_T, ident_tissue) %>% 
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1), 
                   Tumor_reactive13 = mean(Tumor_reactive13), 
                   Exhaustion5 = mean(Exhaustion5), 
                   cytotoxicity6 = mean(cytotoxicity6)) 

df2 <- df_final %>% pivot_longer(., names_to = "func", values_to = "score", 
                                 cols = c("Stemness_Memory1", "cytotoxicity6", "Tumor_reactive13", "Exhaustion5"))

df2$func <- factor(df2$func, levels = c("Stemness_Memory1", "cytotoxicity6", "Tumor_reactive13", "Exhaustion5"))
df2$celltype_T <- factor(df2$celltype_T, levels = rev(levels_CD8))

stat_4B <- adj_refall(as.data.frame(df2), y = "score", x = "celltype_T",
                      facet_by = "func", method = "wilcox.test", padj = "BH") %>% add_p_lab()
print(stat_4B %>% dplyr::select(func, group2, p, p.adj, n_test, p.adj.sig), n = 100)

ggplot(df2, aes(celltype_T, score, fill = celltype_T)) +
  geom_boxplot() + geom_point(size = 1) +
  facet_wrap(func ~ ., scales = "free", nrow = 2) +
  theme_classic() +
  scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_4B, x = "group2", label = "p.adj.lab", size = 2.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2,
               color = "darkred", fill = "darkred", position = position_dodge(.75)) +
  NoLegend()



## C. top10 TCR umap ----
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(tissue != "post_LN") %>% 
  filter(celltype %in% levels_CD8_abb)

df <- tcr_df %>% filter(celltype %in% levels_CD8_abb) %>% 
  group_by(clinical_3, tissue, CDR3.aa, celltype) %>% 
  dplyr::mutate(cell= n(), 
                clone_size_max2 = case_when(cell == 1 ~ "1", cell > 1 ~ ">=2"), 
                clone_size_max3 = case_when(cell == 1 ~ "1", cell == 2 ~ "2", cell > 2 ~ ">=3")) %>% 
  dplyr::summarise(clone = sum(Clones)) %>% group_by(clinical_3, tissue, CDR3.aa) %>% 
  dplyr::mutate(sum_clone = sum(clone)) %>% filter(!CDR3.aa %in% c("None_None", "None_NA", "NA_None", "NA_NA"))

top5 <- c()
top10 <- c()
for (i in c("pre_tumor", "post_tumor")) {
  for (j in c("ER", 'LR')) {
    df2 <- df %>% filter(tissue == i) %>% filter(clinical_3 == j) %>% arrange(-sum_clone) 
    clones <- df2 %>% pull(CDR3.aa) %>% unique()
    top5_clone <- clones[1:5]
    top10_clone <- clones[6:10]
    
    top5_cell <- tcr_df %>% filter(tissue == i) %>% filter(clinical_3 == j) %>% filter(CDR3.aa %in% top5_clone) %>% 
      dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(Barcode, "-"), `[`, 1), "_", patient, "_", tissue)) %>% pull(barcode4merge2)
    top10_cell <- tcr_df %>% filter(tissue == i) %>% filter(clinical_3 == j) %>% filter(CDR3.aa %in% top10_clone) %>% 
      dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(Barcode, "-"), `[`, 1), "_", patient, "_", tissue)) %>% pull(barcode4merge2)
    
    top5 <- c(top5, top5_cell)
    top10 <- c(top10, top10_cell)
  }
}

CD8@meta.data <- CD8@meta.data %>% 
  dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(barcode, "-"), `[`, 1), "_", patient_name, "_", ident_tissue),
                top_clone = case_when(barcode4merge2 %in% top5 ~ "Top5 expanded",
                                      barcode4merge2 %in% top10 ~ "Top10 expanded",
                                      T ~ "Rest"))

tmp <- subset(CD8, subset = ident_tissue != "post_LN")
tmp@meta.data <- tmp@meta.data %>% dplyr::mutate(res_tissue = paste0(ident.clinical_response_3, "_", ident_tissue))
rasterize(DimPlot(tmp, reduction = "standardintegrationharmony_umap",  raster = F,  pt.size = 0.2, group.by = "top_clone", 
                  order = c("Top5 expanded", "Top10 expanded", "Rest"),
                  cols = c("Top5 expanded" = "darkred", "Top10 expanded" = "darkblue", "Rest" = "grey")), dpi = 500) +
  labs(title = "Top clones in tissue", color = "", x = "umap_1", y = 'umap_2') + NoAxes()


## D. TCR clonality ----
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(tissue != "post_LN") %>% filter(celltype %in% levels_CD8_abb)

df <- tcr_df %>% group_by(celltype, CDR3.aa) %>% dplyr::summarise(n = n()) %>% ungroup() %>% 
  group_by(celltype, n) %>% dplyr::summarise(count= sum(n)) %>% 
  dplyr::mutate(clonality = case_when(n == 1 ~ "1",
                                      n == 2 ~ "2",
                                      n > 2 ~ ">=3")) %>% 
  group_by(celltype, clonality) %>% dplyr::summarise(count = sum(count))

df$clonality <- factor(df$clonality, levels = rev(c("1", "2", ">=3")))
df$celltype <- factor(df$celltype, levels = levels_CD8)

p1 <- ggplot(df, aes(count, clonality, fill = celltype)) + geom_col() +
  theme_classic() + scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  theme(axis.text = element_text(size = 5)) 

p1

## E. TR pre/post volcano plot ----
tmp <- subset(CD8, ident_tissue != "post_LN") %>% subset(., celltype_robust == "TR")
Idents(tmp) <- "ident_tissue"
DefaultAssay(tmp) <- "RNA"

deg <- FindMarkers(tmp, ident.1 = "post_tumor", ident.2 =  "pre_tumor", assay = "RNA")

df <- deg %>% 
  mutate(gene = rownames(.),
         sig = case_when(p_val_adj <= 0.01 & avg_log2FC > log2(1.3) ~ "post_tumor",
                         p_val_adj <= 0.01 & avg_log2FC < -log2(1.3) ~ "pre_tumor",
                         T ~ "NS"),
         ident1 = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "NS", T ~ sig),
         gene_label1 = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "", sig == "NS" ~ "", T ~ gene),
         func = gene2func(gene = gene_label1, celltype = "T"),
         gene_label_2 = case_when(gene_label1 != "" & func == 'functional' ~ gene, T ~ ""),
         ident2 = case_when(gene_label_2 == "" ~ "NS", T ~ ident1))
df2 <- df %>% filter(gene_label1 != "")

ggplot(df2, aes(avg_log2FC, -log10(p_val_adj))) + 
  geom_point(aes(col = ident2)) + geom_text_repel(aes(label = gene_label_2), max.overlaps = 3000, size = 3) + theme_classic() + 
  scale_color_manual(values = c("pre_tumor" = "black", "post_tumor" = "#5F967C", NS = "#ECECEC")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") + 
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype = "dashed", color = "darkgrey") +
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 8)) +  NoLegend() + xlim(c(-1.8, 2.7)) + ylim(c(0, 100))


## F. TR pre/post GSEA GOBP ----
a <- deg %>% mutate(gene = rownames(.))
a <- a[order(a$avg_log2FC, decreasing = T), ] # top = post

ranks <- a %>% filter(!grepl("^MT-|RPS|RPL|KRT|IGL|IGK|IGH", gene)) %>% 
  filter(p_val < 0.01) %>%
  as.data.frame() %>%
  dplyr::select(gene, avg_log2FC)  %>% deframe()


#GOBP gdbr_list
set.seed(42)
fgseaRes_G <- fgsea(pathways = gdbr_list , 
                    stats    = ranks,
                    eps = 0.0, nPermSimple = 10000,
                    minSize  = 3, maxSize  = 500)
fgseaResTidy_G <- fgseaRes_G %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy_G <- fgseaResTidy_G %>% filter(!is.na(.$NES))

min(fgseaResTidy_G$NES)
max(fgseaResTidy_G$NES)

df <- fgseaResTidy_G %>% filter(pval <= 0.05) %>% arrange(-NES)
min(df$NES)
max(df$NES)

print(fgseaResTidy_G %>% dplyr::filter(grepl(T_pat, pathway)) %>%
        dplyr::summarise(n_pathway = dplyr::n(),
                         n_pval_0.05 = sum(pval <= 0.05),
                         n_padj_0.05 = sum(padj <= 0.05))) 

df <- fgseaResTidy_G %>% filter(pval <= 0.05) %>% filter(grepl(T_pat, pathway)) %>%
  filter(!grepl("_THYMIC|REELIN|MONOCYTE|SYNAPSE|B_CELL|HUMORAL", pathway))

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score\n(Pre <-> Post)",
       title="GOBP pathways NES from GSEA \nin Whole tissue", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")


## G. TR CD8 signature scores, pre vs post ----
tmp2 <- CD8 %>% subset(., subset = ident_tissue != "post_LN") %>% subset(., celltype_robust == 'TR')
score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1" = TR_gene1, # oliveira
                                        "Tpex_score" = pex_gene,
                                        "Exhaustion" = exh_gene,
                                        "cytotoxicity" = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>% group_by(patient_ID, patient_name, ident.clinical_response_3, ident_tissue, celltype_T) %>% 
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1), 
                   Tumor_reactive13 = mean(Tumor_reactive13), 
                   Exhaustion5 = mean(Exhaustion5), 
                   cytotoxicity6 = mean(cytotoxicity6))

df2 <- df_final %>% pivot_longer(., names_to = "func", values_to = "score", 
                                 cols = c("Stemness_Memory1", "Tumor_reactive13", "Exhaustion5", "cytotoxicity6")) %>% 
  mutate(func2 = case_when(func == "Stemness_Memory1" ~ "Stemness and memory\n(Galletti,Nat.Immunol,2020)", 
                           func == "Tumor_reactive13"~ "Tumor reactivity\n(Oliveira,Nat,2021)", 
                           func ==  "Exhaustion5" ~ "T cell Exhaustion\n(Beltra,Immunity,2020)", 
                           func == "cytotoxicity6" ~ "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)"))
df2$func2 <- factor(df2$func2, levels = c("Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                                          "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)", 
                                          "Tumor reactivity\n(Oliveira,Nat,2021)",
                                          "T cell Exhaustion\n(Beltra,Immunity,2020)"))

df3 <- df2 %>% 
  filter(!patient_name %in% c("P07", "P08")) %>%
  group_by(func2, patient_name, celltype_T) %>%
  filter(n_distinct(ident_tissue) == 2, n() == 2) %>%
  ungroup()

df3$.pair_id <- paste0(df3$patient_name, "_", df3$celltype_T)

stat_4G <- adj_two(as.data.frame(df3), y = "score", x = "ident_tissue",
                   facet_by = "func2", method = "wilcox.test",
                   paired = TRUE, id = ".pair_id", padj = "BH") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_4G %>% dplyr::select(func2, n, p, p.adj, n_test, p.adj.sig))   # n = number of pairs

ggplot(df3, aes(ident_tissue, score, fill = ident_tissue)) +
  geom_boxplot() + geom_point(size = 1) +
  facet_wrap(func2~., scales = "free", nrow = 2) +
  geom_line(aes(group = interaction(patient_name, celltype_T)), col = "darkgrey", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_4G, label = "p.lab", size = 2.8, tip.length = 0.01) +
  stat_summary(aes(group = ident_tissue), fun=mean, geom="point",  shape = 18, size = 2, color="darkred", fill="darkred",
               position = position_dodge(.75)) + NoLegend()

## H. CD8 heatmap ----
for (i in 1) {
  print(i)
  marker <- c('CXCR3', 'ITGAL', "ITGA4", 'CXCR4', 'CXCR5',  'CCR4', 'CCR6', 'CXCL13', "CCL3", "CCL4", "CX3CR1",'IL21R', 'CD40LG',  
              'EOMES', 'KLF2', 'SELL', "S1PR5", 'TCF7', 'BCL6', 'BHLHE40',  'BACH2', 'ID3', 'BATF3','LEF1', 'CCR7', 'IL7R',  'FOXO1', "FOXP1", "S1PR1",  'RUNX3',  
              'SLAMF7', 'FASLG','CD27',  'GZMK', 'CD28',  'CX3CR1','IFNGR1', 'KLRG1', 'KLRB1','KLF4', 'XCL1','BTLA', 'IL2',  'TNFSF14', 'TNFRSF4', 'TNFRSF9', 
              'STAT3', "IL2RA", 'MYC', 'IRF8',
              'IFNG', 'IL15RA',  'TNFRSF18', 'IL4', "LCK", "GZMB", 'KLRC1', 'KLRD1','SAMD9', "TNFSF10", 'ISG15', 'MX1',  'IFIT1', 'IFIT3', 'IFIT2',  'STAT1', 
              'STAT2', 'IRF7',"TBX21", 'TNF','SOCS3','STAT4', 'IL6', 
              'NFATC2', 'TOX', "STAT5A", 'CSF1',  'TIGIT', 'ID2', 'ENTPD1',  'IRF2',"LAG3", 'PDCD1', 'HAVCR2', 'CD38',  "PRDM1", 'TOP2A', 'MKI67'
  )
  
  r_split <- factor(c(rep("Tissue Homing", 11), 
                      rep("B cell interaction", 2), 
                      rep("Memory formation\nCirculation", 17), 
                      rep("Activation\nEffector Differenciation", 43),
                      rep("Exhaustion", 13),
                      rep("Proliferation", 2)),
                    levels = c("Tissue Homing", "B cell interaction","Memory formation\nCirculation", 
                               "Activation\nEffector Differenciation", "Exhaustion",  "Proliferation"))

  
  tmp <- subset(CD8, ident_tissue != "post_LN") %>% subset(., celltype_robust == "TR")
  mat <- AverageExpression(tmp, assays = "RNA", features = marker, 
                           group.by = c("celltype_T", "ident_tissue"), verbose = F)$RNA %>% 
    as.matrix() %>% t() %>% scale() %>% t()
  
  
  col_order <- c("CD8_EM_GZMK+_pre_tumor", "CD8_EM_GZMK+_post_tumor",
                 "CD8_Tpex_pre_tumor" ,  "CD8_Tpex_post_tumor" ,  
                 "CD8_Tex_pre_tumor", "CD8_Tex_post_tumor" , 
                 "CD8_ISG+_pre_tumor", "CD8_ISG+_post_tumor",
                 "Proliferating_CD8_pre_tumor", "Proliferating_CD8_post_tumor")
  
  mat <- mat[, col_order]
  keep    <- !duplicated(marker)    
  marker  <- marker[keep]
  r_split <- r_split[keep]
  
  r_split <- r_split[marker %in% rownames(mat)]
  marker  <- marker[marker %in% rownames(mat)]
  mat     <- mat[marker, ]           
  mat <- t(mat)
  
  split <- factor(rep(c("CD8_EM_GZMK+", "CD8_Tpex", "CD8_Tex", "CD8_ISG+", "Proliferating_CD8"), each = 2),
                  levels = c("CD8_EM_GZMK+", "CD8_Tpex", "CD8_Tex", "CD8_ISG+", "Proliferating_CD8"))
  
  rowHA <- HeatmapAnnotation(which = "column", show_annotation_name = F, show_legend = T, annotation_name_gp = gpar(fontsize = 10), 
                             gp = gpar(col = "white"),
                             Celltype = split,
                             Tissue = factor(rep(c("pre_tumor", "post_tumor"), 5),
                                             levels = c("pre_tumor", "post_tumor")),
                             col = list(Celltype = c(
                               "CD8_EM_GZMK+" = "#F9DC5C",
                               "CD8_Tpex" = "#D4C95A",
                               "CD8_Tex" = "#B6D261",
                               "CD8_ISG+" = "#DD9F58",
                               "Proliferating_CD8" = "#C1BBDD"),
                               Tissue = c("pre_tumor" = "grey", "post_tumor" = "#00C496"))
  )
  p1 <- Heatmap(t(mat), name = "Expression",
                col = colorRamp2(c(-2, -1, 0, 1, 2), c("#4B0055FF", "#00588BFF", "#009B95FF", "#53CC67FF", "#FDE333FF" )),
                show_column_dend = F, show_row_dend = F, cluster_columns = F, cluster_row_slices = F, cluster_column_slices = F,
                row_order = marker, row_names_gp = gpar(fontsize = 10),
                row_split = r_split,
                row_title_rot = 90, row_title_gp = gpar(fontsize = 10), row_gap = unit(3, "mm"),
                column_order = col_order,
                column_title_side = "top",
                column_names_side = "top", column_names_rot = 90, column_names_centered = F, column_title_gp = gpar(fontsize = 7),
                row_names_side = "left", show_heatmap_legend = T,
                top_annotation = rowHA,
                heatmap_height = unit(30, "cm"), heatmap_width = unit(10, "cm"),
                rect_gp = gpar(col = "black"), show_column_names = F, show_row_names = T) 
  print(p1)
}


## I. intra-TLS Tex/TR, Exhuation association ----
## tls - Tex/TR ratio correlation 
df <- lunit %>% dplyr::select(patient_name, ERLR, tls_area_ca_ratio_post, tls_cluster_cnt_ca_ratio_post, IS_post)

tmp <- CD8 %>% subset(., subset = ident_tissue != "post_LN")
df1 <- tmp@meta.data %>% filter(celltype_robust == "TR") %>% 
  group_by(patient_name, ident_tissue, celltype_T) %>% dplyr::summarise(n = n()) %>% 
  group_by(patient_name, ident_tissue) %>% dplyr::mutate(sum = sum(n), ratio = n/sum) %>% 
  filter(!patient_name %in% c("P07", "P08")) %>%
  filter(celltype_T == "CD8_Tex") %>% filter(ident_tissue == "post_tumor")

test <- merge(df, df1)
print(lm_p(test, x = "tls_area_ca_ratio_post", y = "ratio"))

ggplot(test, aes(tls_area_ca_ratio_post, ratio)) + geom_point(aes(col = ERLR)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = values_response_3) +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") +
  ylab("Tex/TR ratio in post tumor") + NoLegend()

## TLS - exhaustion score correlation 
tmp <- subset(CD8, ident_tissue == "post_tumor") %>% subset(., celltype_T == "CD8_Tex")
score <- AddModuleScore(tmp, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Tumor_reactive1" = TR_gene1, 
                                        "Tpex_score" = pex_gene,
                                        "Exhaustion" = exh_gene,
                                        "cytotoxicity" = cyto_gene2
                        ),
                        name = c("Stemness_Memory", "Tumor_reactive1", "Tpex_score", "Exhaustion", "cytotoxicity"))

df1 <- score@meta.data %>% 
  filter(!(ident_tissue == "post_tumor" & patient_name %in% c("P07", "P08"))) %>%
  group_by(patient_name, celltype_T) %>% 
  dplyr::summarise(Stemness_Memory = mean(Stemness_Memory1),
                   Tumor_reactive1 = mean(Tumor_reactive12),
                   Tpex_score = mean(Tpex_score3),
                   Exhaustion = mean(Exhaustion4),
                   cytotoxicity = mean(cytotoxicity5)) 
lun <- merge(df, df1)

print(lm_p(lun, x = "tls_area_ca_ratio_post", y = "Exhaustion"))

ggplot(lun, aes(tls_area_ca_ratio_post, Exhaustion)) + geom_point(aes(col = ERLR)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = values_response_3) +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") + NoLegend() 

# Supple Fig11. ========================================================================================================================
## A. T/NK feature plot ----
list.rna <- rev(c("S1PR1", "CCR7", "TCF7","BACH2", "SELL",
                  "CD28",  "FOSB","FAS",
                  "GZMK", "EOMES",  "CD27", "PDCD1", "GZMM",
                  "XCL1", "STAT5A", "TNFRSF9",
                  "CXCL13", "TOX",
                  "ZNF683",  "HOPX", "KLRD1", "KIR2DL4", 
                  "IFIT1", "OAS1",  "RSAD2",
                  "MKI67",  "CDC25A"))

DefaultAssay(CD8) <- "RNA"
p1 <- FeaturePlot(CD8, features = rev(list.rna), min.cutoff = "q10", max.cutoff = "q90",
                  reduction = "standardintegrationharmony_umap", ncol = 8, raster = T) & NoAxes()
rasterize(p1, dpi = 500)

## B. T/NK UMAP ----
rasterise(DimPlot(CD8, reduction = "standardintegrationharmony_umap",  raster = F, group.by = "celltype_T", pt.size = 0.3, 
                  order = rev(levels_CD8), cols = rev(values_CD8)), dpi = 500) + 
  labs(title = "") & NoAxes()

## C. CD8 dot plot ----
CD8 <- subset(CD8, subset = celltype_T %in% c("CD8_Tex","CD8_CM","CD8_EM_GZMK+","CD8_ISG+","CD8_naive","CD8_Tpex", "Proliferating_CD8"))

p1 <- DotPlot(CD8, features = list.rna, group.by = "celltype_T", assay = "RNA", dot.scale = 5, col.max = 2) +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = levels_CD8_abb) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 5),
        axis.title = element_blank(),
        legend.position = 'top',
        panel.grid = element_line(size = 0.1, linetype = 2)) + NoLegend() 
p1

## D. TR CD8 signature scores, pre vs post per cell type ----
tmp2 <- CD8 %>% subset(., subset = ident_tissue != "post_LN") %>% subset(., celltype_robust == 'TR')

score <- AddModuleScore(tmp2, assay = "RNA",
                        features = list("Stemness_Memory" = Tscm.cm_gene,
                                        "Anti_PD1_reactive" = TR_gene3,
                                        "Tumor_reactive1" = TR_gene1, # oliveira
                                        "Tpex_score" = pex_gene,
                                        "Exhaustion" = exh_gene,
                                        "cytotoxicity" = cyto_gene2),
                        name = c("Stemness_Memory", "Anti_PD1_reactive", "Tumor_reactive1", "Tpex_score",
                                 "Exhaustion", "cytotoxicity"))

df_final <- score@meta.data %>% group_by(patient_ID, patient_name, ident.clinical_response_3, ident_tissue, celltype_T) %>% 
  dplyr::summarise(Stemness_Memory1 = mean(Stemness_Memory1), 
                   Tumor_reactive13 = mean(Tumor_reactive13), 
                   Exhaustion5 = mean(Exhaustion5), 
                   cytotoxicity6 = mean(cytotoxicity6))

df2 <- df_final %>% pivot_longer(., names_to = "func", values_to = "score", 
                                 cols = c("Stemness_Memory1", "Tumor_reactive13", "Exhaustion5", "cytotoxicity6")) %>% 
  mutate(func2 = case_when(func == "Stemness_Memory1" ~ "Stemness and memory\n(Galletti,Nat.Immunol,2020)", 
                           func == "Tumor_reactive13"~ "Tumor reactivity\n(Oliveira,Nat,2021)", 
                           func ==  "Exhaustion5" ~ "T cell Exhaustion\n(Beltra,Immunity,2020)", 
                           func == "cytotoxicity6" ~ "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)"))
df2$func2 <- factor(df2$func2, levels = c("Stemness and memory\n(Galletti,Nat.Immunol,2020)",
                                          "Cytotoxic_Tem\n(Galletti,Nat.Immunol,2020)", 
                                          "Tumor reactivity\n(Oliveira,Nat,2021)",
                                          "T cell Exhaustion\n(Beltra,Immunity,2020)"))
df2$celltype_T <- factor(df2$celltype_T, levels = rev(levels_CD8))

df3 <- df2 %>%
  filter(!patient_name %in% c("P07", "P08")) %>%
  dplyr::mutate(sample = paste0(celltype_T, "_", ident_tissue)) %>%
  group_by(func2, patient_name, celltype_T) %>%
  filter(n_distinct(ident_tissue) == 2, n() == 2) %>%
  ungroup()

lv_sample <- c("CD8_EM_GZMK+_pre_tumor",     "CD8_EM_GZMK+_post_tumor",
               "CD8_Tpex_pre_tumor",         "CD8_Tpex_post_tumor",
               "CD8_Tex_pre_tumor",          "CD8_Tex_post_tumor",
               "CD8_ISG+_pre_tumor",         "CD8_ISG+_post_tumor",
               "Proliferating_CD8_pre_tumor","Proliferating_CD8_post_tumor")

df3$sample       <- factor(df3$sample, levels = lv_sample)
df3$ident_tissue <- factor(df3$ident_tissue, levels = c("pre_tumor", "post_tumor"))
df3              <- df3 %>% filter(!is.na(sample))

stat_paired <- df3 %>%
  dplyr::select(func2, celltype_T, patient_name, ident_tissue, score) %>%
  tidyr::pivot_wider(names_from = ident_tissue, values_from = score) %>%
  tidyr::drop_na(pre_tumor, post_tumor) %>%
  group_by(func2, celltype_T) %>%
  dplyr::summarise(
    n_pair = dplyr::n(),
    p      = tryCatch(
      stats::wilcox.test(pre_tumor, post_tumor, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_),
    .groups = "drop") %>%
  group_by(func2) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH"), n_test = dplyr::n()) %>%
  ungroup() %>%
  dplyr::mutate(
    group1 = paste0(celltype_T, "_pre_tumor"),
    group2 = paste0(celltype_T, "_post_tumor"),
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

ggplot(df3, aes(sample, score, fill = ident_tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 1) +
  facet_wrap(func2 ~ ., scales = "free", nrow = 2) +
  geom_line(aes(group = interaction(patient_name, celltype_T)),
            col = "darkgrey", linewidth = 0.3) +
  theme_classic() +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 7)) +
  labs(x = "", y = "Scaled Average Expression") +
  ggpubr::stat_pvalue_manual(stat_paired, label = "label",
                             tip.length = 0.01, size = 2.5, bracket.size = 0.3) +
  stat_summary(aes(group = ident_tissue), fun = mean, geom = "point",
               shape = 18, size = 2, color = "darkred", fill = "darkred",
               position = position_dodge(.75)) +
  NoLegend()

## E. ISG score ----
tmp2 <- CD8 %>% subset(., subset = ident_tissue != "post_LN")
score <- AddModuleScore(tmp2, assay = "RNA", slot = "data",
                        features = list("good" = IFNG.GS,
                                        "bad"  = ISG.RS),
                        name = c("good", "bad"))

md <- score@meta.data %>%
  filter(celltype_robust == "TR") %>%        
  mutate(good_z = as.numeric(scale(good1)),  
         bad_z  = as.numeric(scale(bad2)),
         delta  = good_z - bad_z)

df_final <- md %>%
  group_by(patient_name, patient_ID, ident_tissue, celltype_T, ident.clinical_response_3) %>%
  dplyr::summarise(good_z = mean(good_z),  
                   bad_z  = mean(bad_z),
                   ratio  = mean(delta),
                   n_cell = n(), .groups = "drop")

inter <- df_final$ratio %>%  mean()
stat_11E1 <- adj_refall(df_final, y = "ratio", x = "celltype_T", padj = "BH") %>% add_p_lab()
print(stat_11E1 %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

p1 <- ggplot(df_final, aes(celltype_T, ratio, fill = celltype_T)) +
  geom_boxplot() + ylab("delta good/bad IFNG effect") + xlab("") +
  theme_classic() + geom_point(position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = inter) +
  scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  ggpubr::stat_pvalue_manual(stat_11E1, x = "group2", label = "p.adj.lab", size = 2.2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred") +
  NoLegend()


## 
df_final$ident_tissue <- factor(as.character(df_final$ident_tissue), levels = c("pre_tumor", "post_tumor"))
stat_11E2 <- adj_two(as.data.frame(df_final), y = "ratio", x = "ident_tissue",
                     facet_by = "celltype_T", method = "wilcox.test", padj = "BH") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_11E2 %>% dplyr::select(celltype_T, n, p, p.adj, n_test, p.adj.sig))

p2 <- ggplot(df_final, aes(ident_tissue, ratio, fill = ident_tissue)) +
  geom_boxplot() + geom_point(position = position_dodge(width = 0.75)) +
  facet_wrap(~ celltype_T, nrow = 1) +
  theme_classic() + theme(legend.position = "top") +
  ggpubr::stat_pvalue_manual(stat_11E2, label = "p.lab", size = 2.2, tip.length = 0.01) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred") +
  scale_fill_manual(breaks = levels_ident_tissue, values = values_ident_tissue) + xlab("")


p1 / p2



# Fig 5. ============================================================================================================
## A. TCR counts per tissues ----
plot_df <- do.call("rbind", tcr_final$data) %>%  filter(celltype %in% levels_CD8_abb) %>%
  group_by(tissue, celltype) %>%
  summarise(TCR_count = n(), .groups = "drop")

plot_df$celltype <- factor(plot_df$celltype, levels = levels_CD8_abb)

ggplot(plot_df, aes(x = tissue, y = TCR_count, fill = celltype)) +
  geom_col(position = "fill") +                 
  labs(x = "tissue", y = "Proportion of TCR count", fill = "celltype") +
  theme_classic() + scale_fill_manual(breaks = levels_CD8_abb, values = values_CD8_abb) 

## B. CD8 TR prop change ----
tmp <- subset(CD8, ident_tissue != "post_LN") %>% subset(., celltype_robust == "TR")
df <- tmp@meta.data %>% group_by(patient_name, ident.clinical_response_3, ident_tissue, celltype_T) %>% 
  dplyr::summarise(n = n()) %>% ungroup() %>% group_by(patient_name, ident.clinical_response_3, ident_tissue) %>% 
  dplyr::mutate(sum = sum(n), ratio = n/sum) %>% ungroup() %>% dplyr::select(-c(n, sum)) %>% 
  pivot_wider(names_from = "ident_tissue", values_from = "ratio") %>% as.data.frame() %>% 
  filter(!patient_name %in% c("P11", "P07", "P08")) %>% replace(., is.na(.), 0) %>% dplyr::mutate(ratio = post_tumor/pre_tumor)
df$celltype_T <- factor(df$celltype_T, levels = rev(levels_CD8))

ggplot(df, aes(celltype_T, sqrt(ratio), fill = celltype_T)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = F, label.y.npc = 0.85) +
  coord_flip() +
  theme_classic() + scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  labs(x = "", fill = "celltype (n=8)")

ggplot(df, aes(celltype_T, sqrt(ratio), fill = celltype_T)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = F, label.y.npc = 0.85) +
  coord_flip() +
  facet_grid(.~ident.clinical_response_3) +
  theme_classic() + scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  labs(x = "", fill = "celltype (n=8)")


## C. TCR expanded per timepoint and response ----
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(tissue != "post_LN") %>% 
  filter(celltype %in% levels_CD8_abb)

df <- tcr_df %>%
  filter(!CDR3.aa %in% c("None_None", "None_NA", "NA_None", "NA_NA")) %>%
  group_by(clinical_3, tissue, patient, CDR3.aa, celltype) %>%
  dplyr::summarise(cell = n(), .groups = "drop") %>%    
  group_by(clinical_3, tissue, patient, CDR3.aa) %>%
  dplyr::mutate(clone_size = sum(cell),                    
                clone_size_max2 = if_else(clone_size == 1, "1", ">=2"),
                clone_size_max3 = case_when(clone_size == 1 ~ "1",
                                            clone_size == 2 ~ "2",
                                            TRUE ~ ">=3")) %>%
  ungroup() %>%
  dplyr::mutate(tissue_res = paste(tissue, clinical_3, sep = "-"))

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
    
    clones <- df2$CDR3.aa
    top5_clone  <- clones[1:5]
    top10_clone <- clones[6:10]                        
    
    top5_cell <- tcr_df %>%
      filter(tissue == i, clinical_3 == j, CDR3.aa %in% top5_clone) %>%
      dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(Barcode, "-"), `[`, 1),
                                            "_", patient, "_", tissue)) %>% pull(barcode4merge2)
    
    top10_cell <- tcr_df %>%
      filter(tissue == i, clinical_3 == j, CDR3.aa %in% top10_clone) %>%
      dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(Barcode, "-"), `[`, 1),
                                            "_", patient, "_", tissue)) %>% pull(barcode4merge2)
    
    top5  <- c(top5,  top5_cell)
    top10 <- c(top10, top10_cell)
  }
}

CD8@meta.data <- CD8@meta.data %>% 
  dplyr::mutate(barcode4merge2 = paste0(sapply(strsplit(barcode, "-"), `[`, 1), "_", patient_name, "_", ident_tissue),
                top_clone = case_when(barcode4merge2 %in% top5 ~ "Top5 expanded",
                                      barcode4merge2 %in% top10 ~ "Top10 expanded",
                                      T ~ "Rest"))

for (i in c("pre_tumor", "post_tumor")) {
  for (j in c("LR", "ER")) {
    tmp <- subset(CD8, subset = ident_tissue != "post_LN") %>% subset(ident_tissue == i) %>% subset(., ident.clinical_response_3 == j)
    tmp@meta.data <- tmp@meta.data %>% dplyr::mutate(res_tissue = paste0(ident.clinical_response_3, "_", ident_tissue))
    rasterize(DimPlot(tmp, reduction = "standardintegrationharmony_umap",  raster = F,  pt.size = 0.5, group.by = "top_clone", 
                      order = c("Top5 expanded", "Top10 expanded", "Rest"),
                      cols = c("Top5 expanded" = "darkred", "Top10 expanded" = "darkblue", "Rest" = "grey")), dpi = 500) +
      labs(title = "", color = "", x = "umap_1", y = 'umap_2') + NoAxes() + NoLegend()
  }
}

tissues <- c("pre_tumor", "post_tumor")
groups  <- c("LR", "ER")

top_clones <- df %>%
  filter(tissue %in% tissues, clinical_3 %in% groups) %>%
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
  dplyr::mutate(celltype = factor(celltype, levels = levels_CD8_abb))

ggplot(plot_df_top10, aes(tissue_res, cell, fill = celltype)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(breaks = levels_CD8_abb, values = values_CD8_abb) +
  labs(x = "", y = "Fraction of T cells", fill = "celltype",
       title = "Top 10 expanded clones") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## D.pre GZMK_EM GSEA by response ----
tmp <- subset(CD8, subset = ident_tissue == "pre_tumor") %>% subset(., celltype_T == "CD8_EM_GZMK+")
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- "ident.clinical_response_3"

deg <- FindMarkers(tmp, ident.1 = "ER", assay = "RNA")
df <- deg %>% 
  mutate(gene = rownames(.),
         sig = case_when(p_val_adj <= 0.01 & avg_log2FC > log2(1.3) ~ "ER",
                         p_val_adj <= 0.01 & avg_log2FC < -log2(1.3) ~ "LR",
                         T ~ "NS"),
         ident1 = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "NS", T ~ sig),
         gene_label1 = case_when(grepl("^MT-|RPL|RPS|KRT", gene) ~ "", sig == "NS" ~ "", T ~ gene))

test <- df[order(df$avg_log2FC, decreasing = T), ] %>% filter(p_val_adj <= 0.05)
df2 <- df %>% filter(gene_label1 != "")

ggplot(df2, aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(aes(col = ident1)) + geom_text_repel(aes(label = gene_label1), max.overlaps = 3000, size = 3) + theme_classic() + 
  scale_color_manual(values = c("pre_tumor" = "black", "post_tumor" = "#5F967C", NS = "#ECECEC")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") + 
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)), linetype = "dashed", color = "darkgrey") +
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 8)) +  NoLegend() + xlim(c(-1.8, 2.7)) + ylim(c(0, 100))

## E. GZMK+ CD8 related association ----
df <- lunit %>% dplyr::select(patient_name, ERLR, tls_area_ca_ratio_post, tls_cluster_cnt_ca_ratio_post)

## (1) intra TLS cnt ~ GZMK+/TR ratio
tmp <- CD8 %>% subset(., subset = ident_tissue == "post_tumor") %>% subset(., celltype_robust == "TR")
df1 <- tmp@meta.data %>%
  group_by(patient_name, ident_tissue, celltype_T) %>% dplyr::summarise(n = n()) %>%
  group_by(patient_name, ident_tissue) %>% dplyr::mutate(sum = sum(n), ratio = n/sum) %>%
  filter(!patient_name %in% c("P07", "P08")) %>%
  filter(celltype_T == "CD8_EM_GZMK+")
test <- merge(df, df1)

## (2) CXCR5 ~ intra TLS area
tmp <- subset(CD8, ident_tissue == "post_tumor") %>% subset(., celltype_T == "CD8_EM_GZMK+")
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
  group_by(patient_name, celltype_T) %>%
  dplyr::summarise(CXCR5 = mean(CXCR56))

lun <- merge(df, df1) %>% filter(!patient_name %in% c("P07", "P08"))
fit_ratio <- lm(ratio ~ tls_cluster_cnt_ca_ratio_post, data = test)
fit_cxcr5 <- lm(tls_area_ca_ratio_post ~ CXCR5,        data = lun)
.pf <- function(m) { f <- summary(m)$fstatistic; unname(pf(f[1], f[2], f[3], lower.tail = FALSE)) }
p_5E   <- c(ratio = .pf(fit_ratio), CXCR5 = .pf(fit_cxcr5))
print(p_5E)

lab_5E <- function(m, key) sprintf("p = %s\nR2 = %.3f,  n = %d",
                                   format.pval(p_5E[[key]], digits = 2),
                                   summary(m)$r.squared, nobs(m))

## plot 1 : intra TLS cnt ~ GZMK+/TR ratio
ggplot(test, aes(tls_cluster_cnt_ca_ratio_post, ratio)) +
  geom_point(aes(col = ERLR)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = values_response_3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.2, color = "black",
           label = lab_5E(fit_ratio, "ratio")) +
  ylab("Tem_GZMK/TR ratio in post tumor") + NoLegend()

## plot 2 : CXCR5 ~ intra TLS area
ggplot(lun, aes(CXCR5, tls_area_ca_ratio_post)) +
  geom_point(aes(col = ERLR)) + theme_pubr() +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = values_response_3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 3.2, color = "black",
           label = lab_5E(fit_cxcr5, "CXCR5")) +
  NoLegend()



## F. Tumor-LN shared TCR cnt / circos ----
## bar plot
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(celltype %in% levels_CD8_abb)  %>% filter(tissue %in% c("post_LN", "post_tumor")) %>% 
  distinct( .keep_all = TRUE)

overlap <- tcr_df %>% group_by(CDR3.aa, tissue) %>% dplyr::summarise(cell_count = n(), clone = sum(Clones)) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% dplyr::select(-c(cell_count)) %>% 
  pivot_wider(names_from = "tissue", values_from = "clone") %>% 
  replace(., is.na(.), 0) %>% filter(post_LN > 0 & post_tumor > 0) %>% 
  as.data.frame() %>% pull(CDR3.aa)

df <- tcr_df %>% group_by(CDR3.aa, tissue, celltype) %>% dplyr::summarise(cell_count = n(), total_clone = sum(Clones)) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% arrange(CDR3.aa) %>% 
  dplyr::select(-c(cell_count))  %>% filter(CDR3.aa %in% overlap) %>%  
  dplyr::mutate(total_clone2 = case_when(tissue == "post_LN" ~ -total_clone, T ~ total_clone),
                unique = case_when(tissue == "post_LN" ~ -1, T ~ 1))


df3 <-df %>% group_by(tissue, celltype) %>% dplyr::summarise(uniq = length(unique(CDR3.aa)))
df3$celltype <- factor(df3$celltype, levels = rev(levels_CD8))
ggplot(df3, aes(tissue, uniq, fill = celltype)) + geom_col(position = "stack") +
  geom_text(position = position_stack(vjust = 0.5), aes(group = celltype, label = uniq)) + theme_linedraw() +
  scale_fill_manual(breaks = levels_CD8, values = values_CD8) +
  ggtitle("unique overlapped clone")

## circos plot
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(celltype %in% levels_CD8_abb)  %>% filter(tissue %in% c("post_LN", "post_tumor")) %>% 
  distinct( .keep_all = TRUE)

overlap <- tcr_df %>% group_by(CDR3.aa, tissue) %>% dplyr::summarise(cell_count = n(), clone = sum(Clones)) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% dplyr::select(-c(cell_count)) %>% 
  pivot_wider(names_from = "tissue", values_from = "clone") %>% 
  replace(., is.na(.), 0) %>% filter(post_LN > 0 & post_tumor > 0) %>% 
  as.data.frame() %>% pull(CDR3.aa)

df <- tcr_df %>% group_by(CDR3.aa, tissue, celltype) %>% dplyr::summarise(cell_count = n(), total_clone = sum(Clones)) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% arrange(CDR3.aa) %>% 
  dplyr::select(-c(cell_count))  %>% filter(CDR3.aa %in% overlap) %>% 
  dplyr::mutate(unique_clone = 1)

df2 <- df %>% dplyr::select(CDR3.aa, tissue, celltype, total_clone) %>% arrange(CDR3.aa) %>% 
  group_by(CDR3.aa) %>% dplyr::summarise(tissue = list(tissue)) %>% rowwise() %>% mutate(tissue = list(unique(tissue)))

overlap <- df2 %>% group_by(CDR3.aa) %>% dplyr::summarise(across(everything(), ~ list(unique(unlist(.))))) %>%
  dplyr::mutate(tissue_length = map_int(tissue, length)) %>% filter(tissue_length == 2) %>% pull(CDR3.aa)

df2 <- df %>% dplyr::select(c(CDR3.aa, tissue, celltype, total_clone)) %>% 
  filter(CDR3.aa %in% overlap) %>% pivot_wider(., values_from = total_clone, names_from = celltype)
df2[is.na(df2)] <- 0


matrix_list <- list()
cdr3_values <- unique(df$CDR3.aa)

for (cdr3 in cdr3_values) {
  sub_df <- subset(df, CDR3.aa == cdr3)
  tmp <- sub_df %>% mutate(sample = paste0(celltype, "_", tissue)) %>% dplyr::select(sample, unique_clone)
  
  result <- t(tmp)
  colnames(result) <- result[1, ]
  result <- result %>% as.data.frame() %>% subset(rownames(result) == "unique_clone")
  matrix_list[[cdr3]] <- result}

dff <- rbindlist(matrix_list, fill = T) %>%
  replace(., is.na(.), 0) %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate_all(as.numeric) %>%
  dplyr::mutate(sum = rowSums(across(where(is.numeric))))  %>%
  mutate(across(starts_with("V"), ~ . / sum)) %>%
  dplyr::select(-c(sum))

for (i in c(1)) {
  print(i)
  patterns <- apply(dff, 2, function(col) paste(as.numeric(col != 0), collapse = ""))
  
  unique_patterns <- unique(patterns)
  grouped_cols <- lapply(unique_patterns, function(pattern) {
    col_indices <- which(patterns == pattern)
    rowSums(dff[, col_indices, drop = FALSE])  
  })
  
  names(grouped_cols) <- paste0("Group_", seq_along(grouped_cols))  
  new_dff <- as.data.frame(grouped_cols)
  
  pattern_order <- c("CM", "EM_GZMK", "Tex", "ISG+", "Proliferating")
  reorder_by_nonzero_pattern <- function(df, pattern_order) {
    col_names <- colnames(df)
    patterns <- sapply(col_names, function(col) {
      paste(as.numeric(df[[col]] != 0), collapse = "")})
    
    sorted_patterns <- patterns[order(patterns, decreasing = TRUE)]
    sorted_col_names <- names(sorted_patterns)
    ordered_cols <- unlist(lapply(pattern_order, function(pattern) {
      grep(pattern, sorted_col_names, value = TRUE)}))
    
    remaining_cols <- setdiff(sorted_col_names, ordered_cols)
    final_order <- c(ordered_cols, remaining_cols)
    df[, final_order, drop = FALSE]}
  
  dff <- reorder_by_nonzero_pattern(new_dff, pattern_order)
  sectors <- c("CD8_naive_post_LN", "CD8_CM_post_LN", "CD8_EM_GZMK+_post_LN", "CD8_Tex_post_LN", "CD8_ISG+_post_LN" ,
               paste0(levels_CD8_abb, "_post_tumor"))
  dff <- dff[sectors, ]
  
  pattern <- apply(as.data.frame(t(dff)), 1, function(x) paste(as.integer(x != 0), collapse = ""))
  df_summarized <- as.data.frame(t(dff)) %>% 
    mutate(pattern = pattern) %>%
    group_by(pattern) %>%
    summarise(across(everything(), sum)) %>%
    dplyr::select(-pattern) %>% as.data.frame() %>% t() %>%
    as.data.frame()
  
  dff <- df_summarized
  
  df$celltype <- factor(df$celltype, levels_CD8_abb)
  sector_widths <- df %>% group_by(celltype, tissue) %>% dplyr::summarise(unique_clone = sum(unique_clone)) %>% ungroup() %>% 
    dplyr::mutate(sum = sum(unique_clone), ratio = unique_clone/sum) %>% 
    arrange(tissue, celltype) %>% pull(ratio)
  
  colors <- c("#ECC067", "#FDFD96", "#F9DC5C", "#B6D261", "#DD9F58",
              "#ECC067", "#FDFD96", "#F9DC5C", "#D4C95A", "#B6D261", "#DD9F58","#C1BBDD")
  
  gap.after <- c(1, 1, 1, 1, 10, 1, 1, 1, 1, 1, 1, 10)
  
  circos.par(cell.padding = c(0.02, 0, 0.02, 0), gap.after = gap.after)
  circos.initialize(factors = sectors, xlim = c(0, 1), sector.width = sector_widths / sum(sector_widths))
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    sector_index <- CELL_META$sector.index 
    xlim <- CELL_META$xlim  
    ylim <- CELL_META$ylim  
    
    circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], border = "black", col = NA)
    
    circos.text(mean(xlim), mean(ylim), labels = sector_index, facing = "clockwise", niceFacing = TRUE, cex = 0.6)
  }, bg.border = NA)
  
  for (z in 1:19) {
    for (i in 1:5) {
      for (j in 6:12) {
        if (i == j) {
          next}
        
        if (dff[i, z] == 0 || dff[j, z] == 0) {
          next}
        
        if (z == 1) {
          circos.link(
            sector.index1 = sectors[i],
            point1 = c(0, dff[i, z]),
            sector.index2 = sectors[j],
            point2 = c(0, dff[j, z]),
            col = colors[j],
            border = "black")
        } else {
          circos.link(
            sector.index1 = sectors[i],
            point1 = c(sum(dff[i, 1:z-1]), sum(dff[i, 1:z])),
            sector.index2 = sectors[j],
            point2 = c(sum(dff[j, 1:z-1]), sum(dff[j, 1:z])),
            col = colors[j],
            border = "black")}}}}
  
  circos.clear()
}

## G. intra tumor TCR overlap ---
### network graph 
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(celltype %in% levels_CD8_abb)  %>% 
  filter(tissue %in% c("post_tumor")) %>% filter(clinical_3 == "LR")

overlap <- tcr_df %>% group_by(CDR3.aa, celltype) %>% dplyr::summarise(cell_count = n(), clone = sum(Clones)) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% dplyr::select(-c(cell_count)) %>% 
  pivot_wider(names_from = "celltype", values_from = "clone")

overlap <- overlap[rowSums(!is.na(overlap)) >= 3, ] %>% 
  replace(., is.na(.), 0) %>% as.data.frame() %>% pull(CDR3.aa)

df <- tcr_df %>% group_by(CDR3.aa, celltype) %>% dplyr::summarise(cell_count = n(), total_clone = sum(Clones), unique_clone = 1) %>% ungroup() %>% 
  filter(!CDR3.aa %in% c("NA_NA", "NA_None", "None_NA", "None_None")) %>% arrange(CDR3.aa) %>% 
  dplyr::select(-c(cell_count))  %>% filter(CDR3.aa %in% overlap)

matrix_list <- list()
cdr3_values <- unique(df$CDR3.aa)

for (cdr3 in cdr3_values) {
  sub_df <- subset(df, CDR3.aa == cdr3) 
  tmp <- sub_df %>% dplyr::select(celltype, unique_clone)
  
  result <- t(tmp)
  colnames(result) <- result[1, ]
  result <- result %>% as.data.frame() %>% filter(rownames(.) == "unique_clone")
  matrix_list[[cdr3]] <- result}

dff <- rbindlist(matrix_list, fill = T) %>%
  replace(., is.na(.), 0) %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate_all(as.numeric) %>% 
  dplyr::mutate(sum = rowSums(across(where(is.numeric))))  %>% 
  mutate(across(starts_with("V"), ~ . / sum)) %>%
  dplyr::select(-c(sum))

test <- rbindlist(matrix_list, fill = T) %>%
  replace(., is.na(.), 0) %>% as.data.frame() %>% mutate_all(as.numeric) 

vertex <- df %>% group_by(celltype) %>% dplyr::summarise(sum = n())

df2 <- matrix(nrow = 7, ncol = 7)
rownames(df2) <- levels_CD8_abb
colnames(df2) <- levels_CD8_abb

for (i in 1:length(levels_CD8_abb)) {
  for (j in 1:length(levels_CD8_abb)) {
    if (i == j) {
      next} else {
        tmp <- test %>% 
          dplyr::select(all_of(c(levels_CD8_abb[i], levels_CD8_abb[j]))) %>% 
          filter(!!sym(levels_CD8_abb[i]) == 1 & !!sym(levels_CD8_abb[j]) == 1)
        print(tmp)
        print(nrow(tmp))
        df2[levels_CD8_abb[i], levels_CD8_abb[j]] <- nrow(tmp)
      }}}


dff <- data.frame(from = character(), to = character(), importance = numeric())
for ( i in 1:7) {
  for (j in 1:7) {
    dff[nrow(dff) + 1,] <- c(colnames(df2)[i], colnames(df2)[j], df2[i, j])}}

dff_unique <- dff %>%
  mutate(
    pair_min = pmin(from, to),  
    pair_max = pmax(from, to)  
  ) %>%
  distinct(pair_min, pair_max, importance, .keep_all = TRUE) %>% 
  dplyr::select(-pair_min, -pair_max) %>% filter(!is.na(importance))


g <- graph_from_data_frame(dff_unique, directed = F)
E(g)$width <- dff_unique$importance
V(g)$size <- vertex$sum[match(V(g)$name, vertex$celltype)]
V(g)$color <- values_CD8_abb[V(g)$name]

plot(g, edge.width = as.numeric(E(g)$width)*0.5,
     vertex.label.cex = 0.8, 
     vertex.size = V(g)$size*0.3,    
     layout = layout_saved) 


## H. Tumor-LN shared phenotype score ----
## LN/T only vs shared calc
label_clone_groups <- function(
    tcr_obj,
    meta,
    target_celltype = "CD8_EM_GZMK+",
    clone_col    = "CDR3.aa",
    celltype_col = "celltype",
    tissue_col   = "tissue",
    tumor_label  = "post_tumor",
    ln_label     = "post_LN"
) {
  if (length(tcr_obj$data) != nrow(meta)) {
    stop("stopped.")
  }
  
  meta2 <- meta %>%
    mutate(sample_id = as.character(seq_len(n())))
  
  df_tcr <- bind_rows(tcr_obj$data, .id = "Sample")
  
  df_join <- df_tcr %>%
    dplyr::select(-any_of(setdiff(intersect(colnames(df_tcr), colnames(meta2)), "Sample"))) %>%
    left_join(meta2, by = "Sample")
  
  df_sub <- df_join %>%
    filter(.data[[celltype_col]] %in% target_celltype)
  
  clone_tissue_summary <- df_sub %>%
    group_by(.data[[clone_col]]) %>%
    summarise(
      has_tumor = any(.data[[tissue_col]] == tumor_label),
      has_ln    = any(.data[[tissue_col]] == ln_label),
      .groups   = "drop"
    ) %>%
    mutate(
      clone_group = dplyr::case_when(
        has_tumor & !has_ln ~ "Tumor_only",
        has_tumor &  has_ln ~ "Tumor_LN_share",
        !has_tumor & has_ln ~ "LN_only",
        TRUE                ~ "Other_or_missing"
      )
    )
  
  return(clone_tissue_summary)
}

clone_groups <- label_clone_groups(
  tcr_obj        = tcr_final,
  meta           = tcr_final$meta,
  target_celltype = levels_CD8_abb,
  clone_col      = "CDR3.aa",
  celltype_col   = "celltype",
  tissue_col     = "tissue",
  tumor_label    = "post_tumor",
  ln_label       = "post_LN"
) %>% filter(!CDR3.aa %in% c("None_NA", "None_None", "NA_None"))

clone_groups$clone_group %>% table()

tcr_df <- do.call("rbind", tcr_final$data) %>%
  filter(celltype %in% levels_CD8_abb) %>%
  filter(tissue %in% c("post_tumor", "post_LN")) %>%
  distinct(.keep_all = TRUE)

tcr_with_group <- tcr_df %>%
  left_join(clone_groups %>% dplyr::select(CDR3.aa, clone_group), by = "CDR3.aa") %>% 
  mutate(
    barcode4merge2 = paste0(
      sapply(strsplit(Barcode, "-"), `[`, 1),
      "_", patient, "_", tissue)
  )

tcr_clone_map <- tcr_with_group %>%
  distinct(barcode4merge2, .keep_all = TRUE) %>%
  dplyr::select(barcode4merge2, clone_group)

CD8@meta.data <- CD8@meta.data %>%
  dplyr::select(-any_of("clone_group")) %>%   
  mutate(b = rownames(.)) %>%
  mutate(
    barcode4merge2 = paste0(
      sapply(strsplit(barcode, "-"), `[`, 1),
      "_", patient_name, "_", ident_tissue)
  ) %>%
  left_join(tcr_clone_map, by = "barcode4merge2") %>%
  mutate(
    clone_group = tidyr::replace_na(clone_group, "NA")) %>%
  column_to_rownames(var = "b")


tmp <- subset(CD8, subset = ident_tissue != "pre_tumor") %>% subset(., subset = clone_group != "NA") %>% 
  subset(., subset = celltype_T == "CD8_EM_GZMK+")
tmp <- AddModuleScore(tmp, assay = "RNA",
                      features = list("Stemness_Memory" = Tscm.cm_gene,
                                      "Exhaustion" = exh_gene,
                                      "cytotoxicity" = cyto_gene2,
                                      "Tumor_reactive1" = TR_gene1),
                      name = c("Stemness_Memory", "Exhaustion", "cytotoxicity", "Tumor_reactive1"))

tmp$clone_group %>% table()
tmp$clone_group %>% is.na() %>% table()
tmp$celltype_T %>% table()
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
tmp <- subset(CD8, subset = ident_tissue == "post_tumor")
DefaultAssay(tmp) <- "RNA"
Idents(tmp) <- "celltype_T"

deg <- FindMarkers(tmp, ident.1 = "CD8_EM_GZMK+", assay = "RNA")

a <- deg %>% mutate(gene = rownames(.))
a <- a[order(a$avg_log2FC, decreasing = T), ] # top = post

ranks <- a %>% filter(!grepl("^MT-|RPS|RPL|KRT|IGL|IGK|IGH", gene)) %>% 
  filter(p_val < 0.01) %>%
  as.data.frame() %>%
  dplyr::select(gene, avg_log2FC)  %>% deframe()

#GOBP gdbr_list
set.seed(42)
fgseaRes_G <- fgsea(pathways = gdbr_list , 
                    stats    = ranks,
                    eps = 0.0, nPermSimple = 10000,
                    minSize  = 3, maxSize  = 500)
fgseaResTidy_G <- fgseaRes_G %>% as_tibble() %>% arrange(desc(NES))
fgseaResTidy_G <- fgseaResTidy_G %>% filter(!is.na(.$NES))

print(fgseaResTidy_G %>% dplyr::filter(grepl(T_pat, pathway)) %>%
        dplyr::summarise(n_pathway = dplyr::n(),
                         n_pval_0.01 = sum(pval <= 0.01),
                         n_padj_0.05 = sum(padj <= 0.05)))  

df <- fgseaResTidy_G %>% filter(pval <= 0.01) %>% filter(grepl(T_pat, pathway)) %>%
  filter(!grepl("_THYMIC|REELIN|MONOCYTE|SYNAPSE|B_CELL|NEURO|NECRO", pathway)) %>%
  filter(NES > 0)

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score\n(Rest <-> CD8_EM_GZMK+)",
       title="GOBP pathways NES from GSEA \nin Whole tissue", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")

## B. Tumor-LN TCR similarity ----
tcr_df <- do.call("rbind", tcr_final$data) %>% filter(celltype %in% c(levels_CD8_abb))  %>% filter(tissue %in% c("post_LN", "post_tumor"))

Type <- c("CD8_naive", "CD8_CM", "CD8_EM_GZMK+", "CD8_Tex", "CD8_Tpex", "CD8_ISG+", "Proliferating_CD8")
Tissue <- c("post_LN", "post_tumor")
response <- c("LR", "ER")

tmp <- tcr_final
if (!is.null(Type)) {tmp <- repFilter(tmp, .method = "by.meta", .query = list(celltype = include(Type)))}
if (!is.null(Tissue)) {tmp <- repFilter(tmp, .query = list(tissue = include(Tissue)), .method = "by.meta")}
if (!is.null(response)) {tmp <- repFilter(tmp, .query = list(clinical_3 = include(response)), .method = "by.meta")}

df <- repOverlap(tmp$data, .method = "tversky", .verbose = F)
ms <- tmp$meta %>% dplyr::select(Sample, patient, tissue, celltype)

df2 <- df %>% as.data.frame() %>% rownames_to_column(var = "source") %>%
  pivot_longer(cols = -source, names_to = "target", values_to = "Similarity") %>%
  filter(!is.na(Similarity)) %>%
  left_join(ms %>% rename_with(~ paste0("source_", .x), -Sample), by = c("source" = "Sample")) %>%
  left_join(ms %>% rename_with(~ paste0("target_", .x), -Sample), by = c("target" = "Sample")) %>%
  filter(source_patient == target_patient, source_tissue != target_tissue) %>%
  dplyr::mutate(cells = paste(pmin(source_celltype, target_celltype),
                              pmax(source_celltype, target_celltype), sep = "|"))

a <- df2 %>% group_by(cells) %>%  dplyr::summarise(Mean = median(Similarity)) %>% arrange(-Mean) %>% pull(cells)
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
