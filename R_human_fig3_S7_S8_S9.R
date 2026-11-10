# Fig3. ================================================================================================================================
## B. GSEA of whole immune cells ----
tmp <- subset(adata, Tissue != "post_LN") %>% subset(., Immune == "immune")
Idents(tmp) <- "Tissue"
deg <- FindMarkers(tmp, ident.1 = "post_tumor", ident.2 = "pre_tumor", assay = "RNA",
                   logfc.threshold = 0.2, min.pct = 0.1)

a <- deg %>% mutate(gene = rownames(.))
a <- a[order(a$avg_log2FC, decreasing = TRUE), ]

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

df <- fgseaResTidy_G %>% filter(pval <= 0.05) %>% filter(grepl(immune_pat, pathway)) %>%
  filter(row_number() %in% c(1:15, (n() - 15):n()))

ggplot(df, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = pval)) +
  coord_flip() +
  labs(x = "", y = "Normalized Enrichment Score\n(Pre <-> Post)",
       title = "GOBP pathways NES from GSEA \nin Whole tissue", fill = "pval") +
  theme_classic() + scale_fill_gradient(low = "red", high = "lightgrey", limits = c(0, 0.1)) +
  theme(axis.text = element_text(size = 8, face = "bold"), axis.title = element_text(size = 10),
        plot.title = element_text(size = 10), legend.position = "right")


## C. B cell annotation, monocle3 ----
rasterize(DimPlot(B, reduction = "umap", group.by = "level2",
                  order = rev(levels(B$level2)), cols = B@misc$cols$level2), dpi = 500) +
  NoAxes() + labs(title = "") 

p1 <- plot_cells(cds, show_trajectory_graph = TRUE, trajectory_graph_segment_size = 1,
                 label_principal_points = FALSE, label_branch_points = FALSE, labels_per_group = TRUE,
                 label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_roots = TRUE,
                 rasterize = FALSE, alpha = 0.5,
                 color_cells_by = "pathway", trajectory_graph_color = "black",
                 norm_method = "log", cell_size = 0.5, graph_label_size = 3) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15), legend.position = "top") +
  scale_color_manual(breaks = names(metadata(cds)$cols$pathway), values = metadata(cds)$cols$pathway) +
  labs(x = "", y = "") + NoAxes() 

rasterise(p1, dpi = 500)


## D. EF/GC score feature plot ----
tmp <- subset(B, subset = Tissue != "post_LN")

GC <- EFGC_ASC %>% filter(Group == "GC") %>% arrange(LogFC) %>%
  filter(LogFC <= -1) %>% filter(P_adj <= 0.05) %>% pull(Gene)
EF <- EFGC_ASC %>% filter(Group == "EF") %>% arrange(-LogFC) %>%
  filter(LogFC >= 1) %>% filter(P_adj <= 0.05) %>% pull(Gene)

score <- AddModuleScore(tmp, assay = "RNA",
                        features = list("GC_related" = GC,
                                        "EF_related" = EF,
                                        "rafi_TIL_specific" = TILB_gene),
                        name = c("GC_related", "EF_related", "rafi_TIL_specific"))

Minval <- min(c(score$GC_related1, score$EF_related2, score$rafi_TIL_specific3))

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GC_related2        = GC_related1 - Minval,
                EF_related3        = EF_related2 - Minval,
                rafi_TIL_specific4 = rafi_TIL_specific3 - Minval,
                GCEF_ratio         = GC_related2 / (EF_related3 + 0.001))

tmp2 <- subset(score, level2 %in% c("B_naive", "B_PreGC", "B_Act_EGR1", "B_Act_NR4A2",
                                    "B_SwStromal", "B_AtM", "B_GCLZ", "B_GCDZ"))
rasterize(FeaturePlot(tmp2, features = c("GC_related1", "EF_related2"), reduction = "umap",
                      max.cutoff = "q75", min.cutoff = "q25"), dpi = 500) & NoAxes()

## E. GC/EF ratio by path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))
df$pathway   <- droplevels(factor(df$pathway))
df$log2_GCEF <- log2(df$GCEF_ratio)

stat_3E <- adj_two(df, y = "log2_GCEF", x = "pathway", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_3E %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(pathway, log2_GCEF, fill = pathway), col = "black") +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_3E, label = "p.lab", size = 3.2, tip.length = 0.01) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(values = B@misc$cols$pathway[c("Path1", "Path2")]) +
  labs(y = "log2(GC/EF ratio)")


## F. HPV specificity along path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))
df$pathway <- droplevels(factor(df$pathway))

df_3F <- df %>% dplyr::mutate(log2_HPV = log2(rafi_TIL_specific4)) %>%
  dplyr::filter(is.finite(log2_HPV), is.finite(pseudotime))
stat_3F <- dplyr::bind_rows(
  lm_p(df_3F %>% dplyr::filter(pathway == "Path1"), x = "pseudotime", y = "log2_HPV") %>%
    dplyr::mutate(pathway = "Path1"),
  lm_p(df_3F %>% dplyr::filter(pathway == "Path2"), x = "pseudotime", y = "log2_HPV") %>%
    dplyr::mutate(pathway = "Path2"))
print(stat_3F %>% dplyr::select(pathway, n, p, r2))

ggplot(df, aes(pseudotime, log2(rafi_TIL_specific4), col = pathway)) +
  rasterize(geom_point(alpha = 0.2, size = 0.3), dpi = 500) + geom_smooth(method = "lm") +
  scale_color_manual(values = B@misc$cols$pathway[c("Path1", "Path2")]) +
  stat_poly_eq(aes(label = paste("p =", format.pval(after_stat(p.value), digits = 2), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, label.y = c(0.98, 0.88), size = 3) +
  theme_classic()


## G. Cell proportion along pseudotime ----
tmp <- subset(B, subset = Tissue != "post_LN")

df <- tmp@meta.data %>% filter(pathway %in% c("Path1")) %>%
  filter(!level2 %in% c("Plasma", "B_GCDZ"))
df$level2 <- factor(as.character(df$level2), levels = rev(levels(B$level2)))

ggplot(df, aes(pseudotime, group = level2, fill = level2, col = level2)) +
  geom_density(position = "fill", alpha = 0.7) +
  facet_wrap(pathway ~ ., scales = "free_x", ncol = 1) +
  theme_minimal() +
  scale_color_manual(breaks = levels(B$level2), values = B@misc$cols$level2) +
  scale_fill_manual(breaks = levels(B$level2), values = B@misc$cols$level2) 


## H. HPV specificity by cell type ----
tmp <- subset(score, subset = pathway == "Path1")

df <- tmp@meta.data %>% filter(level2 != "B_GCDZ") %>%
  mutate(level2 = as.character(level2)) %>%
  group_by(Patient, Tissue, level2) %>% dplyr::summarise(HPV = mean(rafi_TIL_specific4), .groups = "drop")

order <- df %>% group_by(level2) %>% dplyr::summarise(mean = mean(HPV)) %>% arrange(-mean) %>% pull(level2)
df$level2 <- factor(df$level2, levels = order)

stat_3H <- adj_refall(df, y = "HPV", x = "level2", padj = "BH") %>% add_p_lab()
print(stat_3H %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(level2, HPV, fill = level2)) +
  geom_hline(yintercept = mean(df$HPV), lty = "dashed") +
  geom_boxplot() + geom_point() +
  scale_fill_manual(breaks = order, values = B@misc$cols$level2) + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_3H, x = "group2", label = "p.adj.lab", size = 2.5)


## I. SHM of B-AtM ----
df <- vdj %>% filter(Tissue != "post_LN") %>% filter(pathway == "Path1") %>% filter(level2 == "B_AtM") %>%
  mutate(Tissue         = factor(Tissue, levels = c("pre_tumor", "post_tumor")),
         isotype_switch = factor(isotype_switch, levels = c("UnSw", "Sw")))

stat_3I <- adj_two(df, y = "mu_frerm(tlsq", x = "Tissue", facet_by = "isotype_switch",
                   method = "t.test", padj = "BH") %>% add_p_lab(prefix = "t.test p = ")
print(stat_3I %>% dplyr::select(isotype_switch, group1, group2, n, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(Tissue, mu_freq, fill = Tissue)) + geom_boxplot() +
  facet_grid(~isotype_switch) + theme_bw() +
  scale_fill_manual(values = B@misc$cols$Tissue) +
  ggpubr::stat_pvalue_manual(stat_3I, label = "p.lab", size = 3, tip.length = 0.01)

## K. intra-TLS score ----
df2 <- lunit %>%
  dplyr::select(Patient, Response2, tls_area_ca_ratio_pre, tls_area_ca_ratio_post) %>%
  pivot_longer(cols = c("tls_area_ca_ratio_post", "tls_area_ca_ratio_pre"),
               names_to = "tissue", values_to = "tls_score") %>%
  dplyr::mutate(tissue = case_when(tissue == "tls_area_ca_ratio_post" ~ "post", TRUE ~ "pre")) %>%
  filter(!is.na(tls_score)) %>%
  group_by(Patient) %>% filter(n() == 2) %>% ungroup()     

df2$tissue <- factor(df2$tissue, levels = c("pre", "post"))

stat_3K <- adj_two(df2, y = "tls_score", x = "tissue",
                   method = "wilcox.test", paired = TRUE, id = "Patient") %>%
  add_p_lab(prefix = "wilcox p = ")
print(stat_3K %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df2, aes(tissue, tls_score)) + geom_boxplot() +
  ggpubr::stat_pvalue_manual(stat_3K, label = "p.lab", size = 3.2, tip.length = 0.01) +
  theme_classic() + geom_point()


## L. Association between intra-TLS score and pathologic tumor regression ----
df2 <- lunit %>%
  dplyr::select(Patient, Response2, tls_area_ca_ratio_post, tumor_regression) %>%
  filter(!Patient %in% c("P010", "P007"))    

stat_3L <- lm_p(df2, x = "tls_area_ca_ratio_post", y = "tumor_regression")
print(stat_3L)

ggplot(df2, aes(tls_area_ca_ratio_post, tumor_regression)) + geom_point(aes(col = Response2)) +
  geom_smooth(method = "lm", col = "darkred") +
  scale_color_manual(values = adata@misc$cols$Response2) +
  theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 8), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black") + NoLegend()


# Supple Fig7.================================================================================================================================
## A. whole cell annotation ----
DimPlot(adata, reduction = "umap", raster = FALSE, group.by = "level1", pt.size = 0.2, label = FALSE,
        order = rev(levels(adata$level1)), cols = adata@misc$cols$level1) + NoAxes() + labs(title = "")


## B. level1 dot plot ----
list.rna <- c("CPA3",  "KIT", "FCGR3A", "CD14", "CD68", "APOE", "C1QC", "LAMP3", "CD1C",
              "TRDC", "NCAM1", "CD3E", "CD8A", "PDCD1", "CD4", "FOXP3", "IL2RA", "CD19", "MS4A1",
              "MZB1", "IGHG1", "IGHA1", "EPCAM", "MUC1", "PECAM1", "VWF", "COL1A1", "MKI67", "TOP2A")

dot_order <- c(setdiff(levels(adata$level1), "Proliferating"), "Proliferating")

DotPlot(adata, features = list.rna, group.by = "level1", dot.scale = 7, assay = "RNA") +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = rev(dot_order)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_line(linewidth = 0.1, linetype = 2)) + NoLegend()

## C. level1 feature plot ----
rasterize(FeaturePlot(adata, reduction = "umap", raster = TRUE,
                      min.cutoff = "q5", max.cutoff = "q95",
                      features = list.rna, ncol = 7), dpi = 500) & NoAxes() +
  theme(text = element_text(size = 5),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.1, "cm"))

## D. level1 proportion change ----
tmp <- subset(adata, Tissue != "post_LN") %>% subset(., Immune == "immune")
df <- tmp@meta.data %>%
  mutate(across(c(Patient, Response2, Tissue, level1), as.character)) %>%
  group_by(Patient, Response2, Tissue, level1) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(Patient, Response2, Tissue) %>%
  dplyr::mutate(sum = sum(n), ratio = n / sum) %>% ungroup() %>% dplyr::select(-c(n, sum)) %>%
  pivot_wider(names_from = "Tissue", values_from = "ratio") %>% as.data.frame() %>%
  filter(!Patient %in% c("P004", "P010", "P007")) %>%   
  replace(., is.na(.), 0) %>% dplyr::mutate(ratio = post_tumor / pre_tumor) %>%
  filter(Response2 %in% c("ER", "LR"))
df$level1 <- factor(df$level1, levels = rev(levels(adata$level1)))

ggplot(df, aes(level1, sqrt(ratio), fill = level1)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = FALSE, label.y.npc = 0.85) +
  coord_flip() +
  theme_classic() + scale_fill_manual(breaks = levels(adata$level1), values = adata@misc$cols$level1) +
  labs(x = "", fill = "celltype (n=8)")

ggplot(df, aes(level1, sqrt(ratio), fill = level1)) +
  geom_hline(yintercept = 1, lty = "dashed") +
  geom_boxplot() + geom_point() +
  stat_compare_means(method = "wilcox", ref.group = ".all.", hide.ns = FALSE, label.y.npc = 0.85) +
  coord_flip() +
  facet_grid(. ~ Response2) +
  theme_classic() + scale_fill_manual(breaks = levels(adata$level1), values = adata@misc$cols$level1) +
  labs(x = "", fill = "celltype (n=8)")

# Supple Fig8. ================================================================================================================================
## A. B feature plot ----
list.rna <- c("IGHD", "FCER2", "TCL1A", "YBX3",  "ZNF331", "CCR7", "KLF2",
              "HSPA1A", "EGR1", "NR4A2",
              "CD69", "NFKB1", "VPS37B",
              "CRIP1", "ITGB1",
              "FCRL5", "IFI44L", "TBX21", "ITGAX", "TLR7", "STAT1", "ZEB2",
              "DUSP4", "SYK", "HHEX", "MYC", "FOXO1", "SOX4",
              "CXCL13", "CXCR5", "CXCR3", "CCL4", "CCL5", "CD40", "IL21R", "AICDA", "MIR155HG", "EBI3", "TRAF1",
              "BCL6", "IRF4", "NEIL1", "SUGCT", "MEF2B", "STMN1", "MKI67",
              "CD27", "XBP1", "PRDM1",
              "IGHG1", "IGHA1")

rasterize(FeaturePlot(B, reduction = "umap", raster = TRUE, min.cutoff = "q5", max.cutoff = "q95", pt.size = 1.5,
                      features = list.rna, ncol = 6), dpi = 500) & NoAxes() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width = unit(0.1, "cm"))

## B. B dot plot ----
DotPlot(B, features = list.rna, group.by = "level2", dot.scale = 5) +
  scale_color_viridis(oob = scales::squish) +
  scale_y_discrete(limits = rev(levels(B$level2))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        panel.grid = element_line(linewidth = 0.1, linetype = 2)) + NoLegend()


# Supple Fig9. ================================================================================================================================
## A. HPV specificity per path ----
df <- score@meta.data %>% filter(pathway %in% c("Path1", "Path2"))
df$pathway  <- droplevels(factor(df$pathway))
df$log2_HPV <- log2(df$rafi_TIL_specific4)

stat_9A <- adj_two(df, y = "log2_HPV", x = "pathway", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9A %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(pathway, log2_HPV, fill = pathway), col = "black") +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_9A, label = "p.lab", size = 3.2, tip.length = 0.01) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(values = B@misc$cols$pathway[c("Path1", "Path2")]) +
  labs(y = "log2(TIL-B specificity score)")


## B. association of GC - HPV or SHM ----
tmp <- subset(B, subset = Tissue != "post_LN") %>% subset(., pathway == "Path1")

GC <- EFGC_ASC %>% filter(Group == "GC") %>% arrange(LogFC) %>%
  filter(LogFC <= -1) %>% filter(P_adj <= 0.05) %>% pull(Gene)
EF <- EFGC_ASC %>% filter(Group == "EF") %>% arrange(-LogFC) %>%
  filter(LogFC >= 1) %>% filter(P_adj <= 0.05) %>% pull(Gene)

score <- AddModuleScore(tmp, assay = "RNA",
                        features = list("GC_related" = GC,
                                        "EF_related" = EF,
                                        "rafi_TIL_specific" = TILB_gene),
                        name = c("GC_related", "EF_related", "rafi_TIL_specific"))

Minval <- min(c(score$GC_related1, score$EF_related2, score$rafi_TIL_specific3))

score@meta.data <- score@meta.data %>%
  dplyr::mutate(GC_related2        = GC_related1 - Minval,
                EF_related3        = EF_related2 - Minval,
                rafi_TIL_specific4 = rafi_TIL_specific3 - Minval,
                GCEF_ratio         = GC_related2 / (EF_related3 + 0.001))

meta <- score@meta.data %>%
  dplyr::mutate(cell_id = rownames(.)) %>%
  dplyr::select(cell_id, GC_related2, EF_related3, rafi_TIL_specific4, GCEF_ratio)


## GC - HPV
meta2 <- vdj %>% filter(Tissue != "post_LN") %>%
  dplyr::select(cell_id, Patient, level2, Tissue, Response2, mu_freq) %>%
  left_join(meta, by = "cell_id")

df <- meta2 %>% filter(!is.na(GC_related2)) %>% filter(level2 != "B_GCDZ")
df <- df %>% group_by(Patient, level2) %>%
  dplyr::summarise(GC = mean(GC_related2), Mu = mean(mu_freq), HPV = mean(rafi_TIL_specific4),
                   GCEF = mean(GCEF_ratio), .groups = "drop")

df_9B1   <- df
stat_9B1 <- lm_p(df_9B1, x = "GC", y = "HPV")
print(stat_9B1)

ggplot(df, aes(GC, HPV)) + geom_point() + geom_smooth(method = "lm", color = "darkred") + theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black")

## GC - Mu
meta2 <- vdj %>% filter(Tissue != "post_LN") %>%
  dplyr::select(cell_id, Patient, level2, Tissue, Response2, mu_freq, isotype_switch) %>%
  left_join(meta, by = "cell_id")

df <- meta2 %>% filter(!is.na(GC_related2)) %>% filter(level2 != "B_GCDZ")
df <- df %>% group_by(Patient, level2, isotype_switch) %>%
  dplyr::summarise(GC = mean(GC_related2), Mu = mean(mu_freq), HPV = mean(rafi_TIL_specific4), .groups = "drop")

df_9B2   <- df
stat_9B2 <- lm_p(df_9B2, x = "GC", y = "Mu")
print(stat_9B2)

ggplot(df, aes(GC, Mu)) + geom_point() + geom_smooth(method = "lm", color = "darkred") + theme_pubr() +
  stat_poly_eq(aes(label = paste("p =", round(after_stat(p.value), digits = 5), ",\n", after_stat(rr.label))),
               formula = y ~ x, parse = FALSE, color = "black")


## C. GC/EF ratio by response ----
df <- score@meta.data %>% filter(Tissue != "post_LN") %>% filter(pathway == "Path1") %>% filter(level2 == "B_AtM")
df$Tissue <- factor(as.character(df$Tissue), levels = c("pre_tumor", "post_tumor"))

stat_9C <- adj_two(df, y = "GCEF_ratio", x = "Tissue",
                   facet_by = "Response2", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9C %>% dplyr::select(Response2, n, p, p.adj, n_test, p.adj.sig))

ggplot(df, aes(Tissue, GCEF_ratio, fill = Tissue)) + geom_violin() + geom_boxplot() +
  facet_grid(. ~ Response2) + theme_classic() +
  stat_summary(fun = "mean", geom = "point", shape = 18, color = "darkred", size = 3) +
  ggpubr::stat_pvalue_manual(stat_9C, label = "p.lab", size = 3, tip.length = 0.01) +
  scale_fill_manual(values = B@misc$cols$Tissue)


## D. Cell proportion along pseudotime in Path2 ----
tmp <- subset(B, subset = Tissue != "post_LN")

df <- tmp@meta.data %>% filter(pathway %in% c("Path2")) %>%
  filter(!level2 %in% c("Plasma", "B_GCDZ"))
df$level2 <- factor(as.character(df$level2), levels = rev(levels(B$level2)))

ggplot(df, aes(pseudotime, group = level2, fill = level2, col = level2)) +
  geom_density(position = "fill", alpha = 0.7) +
  facet_wrap(pathway ~ ., scales = "free_x", ncol = 1) +
  theme_minimal() +
  scale_color_manual(breaks = levels(B$level2), values = B@misc$cols$level2) +
  scale_fill_manual(breaks = levels(B$level2), values = B@misc$cols$level2) 

## E. GC score by cell type in Path1 ----
df <- score@meta.data %>% filter(pathway %in% c("Path1")) %>% filter(level2 != "B_GCDZ") %>%
  mutate(level2 = as.character(level2))
df1 <- df %>% group_by(Patient, Tissue, level2) %>% dplyr::summarise(Mean = mean(GC_related2), .groups = "drop")
order <- df1 %>% group_by(level2) %>% dplyr::summarise(mean = mean(Mean)) %>% arrange(-mean) %>% pull(level2)

df1$level2 <- factor(df1$level2, levels = order)
df1$level2 <- factor(df1$level2, levels = intersect(levels(B$level2), unique(df1$level2)))

stat_9E <- adj_refall(df1, y = "Mean", x = "level2", padj = "BH") %>% add_p_lab()
print(stat_9E %>% dplyr::select(group2, p, p.adj, n_test, p.adj.sig))

ggplot(df1, aes(level2, Mean, fill = level2), col = "black") +
  geom_boxplot() + geom_point() +
  ggpubr::stat_pvalue_manual(stat_9E, x = "group2", label = "p.adj.lab", size = 2.5) +
  theme_classic() +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2) +
  scale_fill_manual(values = B@misc$cols$level2) + NoLegend()


## F. BCR diversity by cell type in Path1 ----
row_name <- c("UnSw_pre_tumor", "UnSw_post_tumor", "Sw_pre_tumor", "Sw_post_tumor")
col_name <- levels(B$level2)

tmp2 <- vdj %>% filter(Tissue != "post_LN") %>%
  mutate(iso_tissue = paste0(isotype_switch, "_", Tissue),
         sample     = paste0(Patient, "_", Tissue)) %>%        
  filter(pathway == "Path1") %>% filter(level2 != "B_GCDZ")

mat_list <- list()
for (k in unique(tmp2$sample)) {
  tmp <- tmp2 %>% filter(sample == k)
  mat <- matrix(0, nrow = length(row_name), ncol = length(col_name), dimnames = list(row_name, col_name))
  
  for (i in seq_along(row_name)) {
    for (j in seq_along(col_name)) {
      clone_counts <- tmp %>% filter(iso_tissue == row_name[i], level2 == col_name[j]) %>%
        group_by(clone_id) %>% summarise(n = n()) %>% pull(n)
      mat[i, j] <- diversity(clone_counts, index = "invsimpson")
    }
  }
  
  mat_list[[k]] <- as.data.frame(mat) %>%
    dplyr::mutate(sample = k, iso_tissue = rownames(.)) %>%
    pivot_longer(cols = all_of(col_name), names_to = "celltype", values_to = "div") %>% as.data.frame()
}

sample_info <- tmp2 %>% distinct(sample, Patient, Response2)

div <- bind_rows(mat_list) %>% filter(div != Inf) %>%
  left_join(sample_info, by = "sample") %>%
  dplyr::mutate(response       = factor(as.character(Response2), levels = levels(B$Response2)),
                isotype_switch = str_extract(iso_tissue, "UnSw|Sw"),
                Tissue         = factor(str_extract(iso_tissue, "pre_tumor|post_tumor"),
                                        levels = c("pre_tumor", "post_tumor")))

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
  scale_fill_manual(values = B@misc$cols$level2) + NoLegend() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


## G. GC score of B-AtM per timepoint ----
df <- score@meta.data %>% filter(pathway == "Path1") %>% filter(level2 == "B_AtM")
df$Tissue <- factor(as.character(df$Tissue), levels = c("pre_tumor", "post_tumor"))

stat_9G <- adj_two(df, y = "GC_related2", x = "Tissue", method = "t.test", padj = "BH") %>%
  add_p_lab(prefix = "t.test p = ")
print(stat_9G %>% dplyr::select(group1, group2, n, p, p.adj, n_test))

ggplot(df, aes(Tissue, GC_related2, fill = Tissue)) +
  geom_violin() + geom_boxplot() + theme_classic() +
  ggpubr::stat_pvalue_manual(stat_9G, label = "p.lab", size = 3.2, tip.length = 0.01) +
  scale_fill_manual(values = B@misc$cols$Tissue) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "darkred", fill = "darkred", position = position_dodge(.75))

