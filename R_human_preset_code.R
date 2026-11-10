# Load packages ------------------------------------------------------------------------------------------------------------------------------------
# R version 4.1.3 
library(openxlsx)        # v4.2.8.1
library(readxl)          # v1.4.5
library(readr)           # v2.1.6
library(Seurat)          # v4.4.0 (objects: SeuratObject v4.1.3 format, Assay v3)
library(dplyr)           # v1.1.4
library(ggplot2)         # v4.0.1
library(circlize)        # v0.4.17
library(tidyr)           # v1.3.1
library(tibble)          # v3.3.0
library(ggpmisc)         # v0.6.3
library(ggpubr)          # v0.6.2
library(stringr)         # v1.6.0
library(viridis)         # v0.6.5
library(ggrepel)         # v0.9.6
library(ComplexHeatmap)  # v2.10.0
library(scales)          # v1.4.0
library(fgsea)           # v1.20.0
library(msigdbr)         # v7.5.1
library(data.table)      # v1.17.8
library(immunarch)       # v0.9.1
library(monocle3)        # v1.3.4
library(igraph)          # v2.2.1
library(ggrastr)         # v1.0.2
library(vegan)           # v2.7-2
library(patchwork)       # v1.3.2

# Define project path ------------------------------------------------------------------------------------------------------------------------------
proj_dir <- "/home/yeju/1.YMC_HNC/HNSC_paper/1.public_relaease/0.R_public"   
data_dir <- file.path(proj_dir, "0.data")
ref_dir  <- file.path(proj_dir, "2.reference")

# Load data ------------------------------------------------------------------------------------------------------------------------------------
## Seurat objects
adata <- readRDS(file.path(data_dir, "HNSCC_whole_adata_public.rds"))
B     <- readRDS(file.path(data_dir, "HNSCC_B_subset_public.rds"))
CD8   <- readRDS(file.path(data_dir, "HNSCC_CD8_subset_public.rds"))

## B cell monocle 3 object
cds <- readRDS(file.path(data_dir, "HNSCC_B_monocle3_cds_public.rds"))

## BCR data
vdj <- read_tsv(file.path(data_dir, "HNSCC_BCR_contigs_public.tsv"))

## TCR data
tcr <- readRDS(file.path(data_dir, "HNSCC_CD8_TCR_immunarch_public.rds"))

## Lunit TLS data
lunit <- read_tsv(file.path(data_dir, "HNSCC_TLS_lunit_public.tsv"))

## Network layout for shared CD8 TCRs (Fig5G)
load(file.path(data_dir, "Fig5G_network_layout.RData"))


# GSEA database ------------------------------------------------------------------------------------------------------------------------------------
m_gobp <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP') # package version: 7.5.1

T_pat = '_T_CELL|IMMUNE|INNATE|ADAPTIVE|INFLAM|INTERLEUKIN|INTERFERON|NECROSIS|APOPTOSIS|SENESCENCE|NATURAL|LYMPHOCYTE|LEUKOCYTE|TRANSFORMING|CHEMO|CYTOTO|CYTOKINE|ANTIGEN|SIGNALING|EXHAUSTION|CHEMOKINE'
immune_pat <- "_T_CELLIMMUNE|IMMUNO|INNATE|ADAPTIVE|INFLAM|INTERLEUKIN|INTERFERON|NECROSIS|APOPTOSIS|SENESCENCE|NATURAL|LYMPHOCYTE|LEUKOCYTE|TRANSFORMING|CHEMO|CYTOTO|CYTOKINE|ANTIGEN|SIGNALING|EXHAUSTION|KINASE|KAPPA|ACUTE|HUMORAL|ANTIBODY|ANTIGEN|_MHC_|AUTOPHA|_B_CELL|_HOST|_CD40|_CELL_CYCLE|VIRUS|VIRAL|SENESCENCE|COMPLEMENT|TUMOR_CELL|DENDRITIC_CELL|EOSINOPHIL|ERK1|_FC_|ISOTYPE|LYMPH|MACROPHAGE|MAST_CELL|MEMORY|MITOTIC|VIRAL|MONOCYTE|MYELOID_CELL|NATURAL_KILLER_CELL|PHAGO|TOLL_LIKE|_TUMOR_|MAPK|PROGRAMMED|STAT|APOPTOTIC|GRANULOCYTE|_T_HELPER|TRANSFORMING_GROWTH_FACTOR|ORGANISM|CELL_KILLING|GERMINAL"

gdbr_list <- split(x = m_gobp$gene_symbol, f = m_gobp$gs_name)


# Define functions ------------------------------------------------------------------------------------------------------------------------------------
gene2func <- function(gene, celltype = "T") {
  
  marker_T <- c("PRF1", "NKG7", "GZMB", "GZMK",  "GZMA" , "GZMM", "GZMH",  "GZMA", "KLRC1",  "KLRB1", "TNF", "IFNG", "GNLY", "CLEC2D", "KLRD1",
                "LY6E", "LY6D", "GATA3", "HOPX", "CXCR3", "IL2RB", "IL2RA",  "STAT3", "CD74",  "JUN", "CD38",  "HLA-DQA1", "HLA-DQA2", "ZAP70", "APOBEC3G", "BHLHE40",
                "TNFRSF9",  "HLA-DRB1", "HLA-DPA1", "HLA-DRA", "CD69", "NFKB1A", "STAT4", "DUSP1", "FOS", "FOSB", "HLA-DMA",
                "HLA-DPB1" , "HLA-DRB5",  "HLA-DQB1", "TNFSF9",  "STAT1", "TAP1", "JUNB", "FOSL2", "TNFSF4", "TNFSF14",
                "CXCR4", "IL7R", "BCL2L11", "IRF4", "SLAMF6", "BCL3", "TXNIP", "SLAMF7", "IL6ST", "PDCD1", "LAG3", "TIGIT", "PRDM1", "GBP1" , "GBP5", "IL10RA",  "HAVCR2", "TOX", "TOX2",
                "ENTPD1", "CTLA4", "CXCL13",  "SOX4", "CD101", "NR4A2", "KLF13", "DUSP2", 
                "ISG20",  "IFITM3", "IFI6", "MX1", "OAS1", "APOBEC3C", "IFI44L", "IFITM1", "IFI44", "MX2", "IFIT3", "OASL", "ISG15", "OAS2", "OAS3", "CCL20", "CCL4",   
                "CCL5", "CCL4L2",  "CCL3L3" , "CXCL8", "CCL3", "IRF1", "XCL1", "IKZF2", "BATF","ZNF683", "RUNX3", "VCAM1", 
                "ZFP36L2", "MYO1E", "ZBTB38", "KLF9", "IL32", "S100A11", "LGALS1", "CD9", "TNFSF10", "ITGB7", "CD63", "IKZF1", "CD82", "ITGB2",
                "ITK", "RELB", "HSP90AB1", "HSP90AA1", "HSPH1", "HSPD1", "HSPA1A", "HSPA1B")
  
  case_when(grepl("IGH|IGK|IGL|TRAV|TRBV", gene) ~ "non_functional",
            grepl("CXCL|CXCR|CCL|CCR|KLF|^IL|DUSP|FLRK|^ITGA|^CD|^ISG|SLAMF|HLA-|^IFI|^IRF|TNFSF|TNFRSF|ZBTB|^GZM|IKZF|^KLR|^SOX|^TLR", gene) ~ "functional",
            gene %in% marker_T ~ "functional",
            TRUE ~ "non_functional")
}

## statistical function
# adj_refall() : ref.group = ".all." 
adj_refall <- function(data, y, x, facet_by = NULL,
                       method = "wilcox.test", padj = "BH",
                       y.frac = 1.06, digits = 2) {
  data <- as.data.frame(dplyr::ungroup(data))
  data <- data[!is.na(data[[y]]) & !is.na(data[[x]]), , drop = FALSE]
  data <- droplevels(data)
  fml  <- stats::as.formula(paste0("`", y, "` ~ `", x, "`"))
  
  st <- if (is.null(facet_by)) {
    ggpubr::compare_means(fml, data = data, method = method, ref.group = ".all.")
  } else {
    ggpubr::compare_means(fml, data = data, group.by = facet_by,
                          method = method, ref.group = ".all.")
  }
  
  st <- st %>% dplyr::select(-dplyr::any_of(c("p.adj", "p.format", "p.signif")))
  st <- if (is.null(facet_by)) {
    st %>% dplyr::mutate(p.adj = stats::p.adjust(p, padj), n_test = dplyr::n())
  } else {
    st %>% dplyr::group_by(.data[[facet_by]]) %>%
      dplyr::mutate(p.adj = stats::p.adjust(p, padj), n_test = dplyr::n()) %>%
      dplyr::ungroup()
  }
  st <- dplyr::mutate(st, p.adj.lab = paste0("p.adj =\n", format.pval(p.adj, digits = digits)),
                      p.adj.sig = .sig(p.adj))
  
  ypos <- .ypos(data, y, facet_by, y.frac)
  if (is.null(facet_by)) dplyr::mutate(st, y.position = ypos$y.position[1])
  else dplyr::left_join(st, ypos, by = facet_by)
}


# adj_two() : unpaired/paired
adj_two <- function(data, y, x, facet_by = NULL, method = c("t.test", "wilcox.test"),
                    paired = FALSE, id = NULL, padj = "BH",
                    y.frac = 1.06, digits = 2) {
  
  method <- match.arg(method)
  data   <- as.data.frame(dplyr::ungroup(data))        
  data   <- data[!is.na(data[[y]]) & !is.na(data[[x]]), , drop = FALSE]
  lv     <- levels(factor(data[[x]]))
  stopifnot(length(lv) == 2)
  if (paired && is.null(id)) stop("paired = TRUE")
  
  data$.facet <- if (is.null(facet_by)) "all" else as.character(data[[facet_by]])
  
  out <- lapply(split(data, data$.facet), function(dd) {
    if (paired) {
      w <- dd[, c(id, x, y)] %>%
        tidyr::pivot_wider(names_from = dplyr::all_of(x), values_from = dplyr::all_of(y),
                           values_fn = mean) %>%
        tidyr::drop_na(dplyr::all_of(lv))
      n <- nrow(w)
      p <- if (n < 2) NA_real_ else tryCatch(
        if (method == "t.test") stats::t.test(w[[lv[1]]], w[[lv[2]]], paired = TRUE)$p.value
        else stats::wilcox.test(w[[lv[1]]], w[[lv[2]]], paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_)
    } else {
      a <- dd[[y]][dd[[x]] == lv[1]]; b <- dd[[y]][dd[[x]] == lv[2]]
      n <- min(length(a), length(b))
      p <- if (n < 2) NA_real_ else tryCatch(
        if (method == "t.test") stats::t.test(a, b)$p.value
        else stats::wilcox.test(a, b, exact = FALSE)$p.value,
        error = function(e) NA_real_)
    }
    data.frame(.facet = dd$.facet[1], group1 = lv[1], group2 = lv[2],
               n = n, p = p, stringsAsFactors = FALSE)
  })
  st <- dplyr::bind_rows(out)
  st$p.adj     <- stats::p.adjust(st$p, padj)         
  st$n_test    <- nrow(st)
  st$p.adj.lab <- paste0(if (method == "t.test") "t.test p.adj = " else "wilcox p.adj = ",
                         format.pval(st$p.adj, digits = digits))
  st$p.adj.sig <- .sig(st$p.adj)
  
  ypos <- .ypos(data, y, facet_by, y.frac)
  if (is.null(facet_by)) { st$y.position <- ypos$y.position[1]; st$.facet <- NULL }
  else { names(st)[names(st) == ".facet"] <- facet_by
  st[[facet_by]] <- as.character(st[[facet_by]])
  ypos[[facet_by]] <- as.character(ypos[[facet_by]])
  st <- dplyr::left_join(st, ypos, by = facet_by) }
  st
}


# adj_list() 
adj_list <- function(data, y, x, comparisons, facet_by = NULL,
                     method = c("t.test", "wilcox.test"), padj = "BH",
                     y.frac = 1.06, step = 0.09, digits = 2) {
  
  method <- match.arg(method)
  data   <- as.data.frame(dplyr::ungroup(data))        
  data   <- data[!is.na(data[[y]]) & !is.na(data[[x]]), , drop = FALSE]
  data$.facet <- if (is.null(facet_by)) "all" else as.character(data[[facet_by]])
  
  st <- dplyr::bind_rows(lapply(split(data, data$.facet), function(dd) {
    dplyr::bind_rows(lapply(comparisons, function(cp) {
      a <- dd[[y]][dd[[x]] == cp[1]]; b <- dd[[y]][dd[[x]] == cp[2]]
      p <- if (min(length(a), length(b)) < 2) NA_real_ else tryCatch(
        if (method == "t.test") stats::t.test(a, b)$p.value
        else stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA_real_)
      data.frame(.facet = dd$.facet[1], group1 = cp[1], group2 = cp[2],
                 n1 = length(a), n2 = length(b), p = p, stringsAsFactors = FALSE)
    }))
  }))
  st$p.adj     <- stats::p.adjust(st$p, padj)
  st$n_test    <- nrow(st)
  st$p.adj.lab <- paste0("p.adj = ", format.pval(st$p.adj, digits = digits))
  st$p.adj.sig <- .sig(st$p.adj)
  
  rng <- .ypos(data, y, facet_by, 1)
  st <- st %>% dplyr::group_by(.facet) %>%
    dplyr::mutate(.k = dplyr::row_number()) %>% dplyr::ungroup()
  rng$.facet <- if (is.null(facet_by)) "all" else as.character(rng[[facet_by]])
  st <- dplyr::left_join(st, rng[, c(".facet", "y.position", ".rng")], by = ".facet")
  st$y.position <- st$y.position + st$.rng * ((y.frac - 1) + step * (st$.k - 1))
  st$.k <- NULL; st$.rng <- NULL
  if (is.null(facet_by)) st$.facet <- NULL else names(st)[names(st) == ".facet"] <- facet_by
  st
}


# lm_p() : F-test p + R^2
lm_p <- function(data, x, y) {
  d <- data.frame(.x = as.numeric(data[[x]]), .y = as.numeric(data[[y]]))
  d <- d[stats::complete.cases(d), ]
  m <- stats::lm(.y ~ .x, data = d)
  f <- summary(m)$fstatistic
  data.frame(x = x, y = y, n = nrow(d),
             p  = unname(stats::pf(f[1], f[2], f[3], lower.tail = FALSE)),
             r2 = summary(m)$r.squared, stringsAsFactors = FALSE)
}



# internal utility
.sig <- function(p) dplyr::case_when(is.na(p) ~ "NA", p <= 1e-4 ~ "****", p <= 1e-3 ~ "***",
                                     p <= 1e-2 ~ "**",  p <= 0.05 ~ "*",   TRUE ~ "ns")

.ypos <- function(data, y, facet_by, y.frac) {
  f <- function(v) data.frame(y.position = max(v, na.rm = TRUE) +
                                diff(range(v, na.rm = TRUE)) * (y.frac - 1),
                              .rng = diff(range(v, na.rm = TRUE)))
  if (is.null(facet_by)) f(data[[y]])
  else data %>% dplyr::group_by(.data[[facet_by]]) %>%
    dplyr::group_modify(~ f(.x[[y]])) %>% dplyr::ungroup()
}

add_p_lab <- function(st, digits = 2, prefix = "p = ") {
  st$p.lab <- paste0(prefix, format.pval(st$p, digits = digits))
  st$p.sig <- .sig(st$p)
  st
}



# AddModuleScore gene sets --------------------------------------------------------------------------------------------------------------------------------
# EF/GC ASC DEG list
# Ma J, et al. A blueprint for tumor-infiltrating B cells across human cancers. Science 2024;384:eadj4857. (Table: "5. DEGs of EF and GC ASCs")
EFGC_ASC <- read.xlsx(file.path(ref_dir, "sci2024_GCEF_DEG_list.xlsx"),
                      sheet = "5. DEGs of EF and GC ASCs", colNames = TRUE, startRow = 2)


# IFNG.GS (IFNG hallmark signature) / ISG.RS (ISG resistance signature)
# Benci JL, et al. Opposing functions of interferon coordinate adaptive and innate immune responses to cancer immune checkpoint blockade. Cell 2019;178:933-948.
# (ISG.RS was originally defined in Benci JL, et al. Cell 2016;167:1540-1554.)
IFNG.GS <- c("TNFSF10", "NLRC5", "CCL5", "PFKP", "GZMA", "IRF9", "NMI", "RNF31", "CD38", "IFNAR2", "PSMB9", "IDO1", "SOD2", "ZBP1", "CD74", "EPSTI1", "PSMB10", "PSMA3", "BANK1", "RAPGEF6", "PARP12", "CXCL11", "RNF213", "TOR1B", "CASP4", "TRIM25", "ITGB7", "PELI1", "RBCK1", "FAS", "LAP3", "SAMHD1", "CFB", "PDE4B", "OGFR", "CASP7", "CMPK2", "CD86", "MVP", "ARL4A", "UPP1", "SAMD9L", "TXNIP", "IL7", "LYSMD2", "B2M", "RTP4", "HLA-DQA1", "BPGM", "CSF2RB", "IRF4", "PTPN2", "GCH1", "FTSJ2", "ST3GAL5", "SRI", "PARP14", "PNP", "AUTS2", "C1R", "NFKBIA", "TNFAIP2", "CCL7", "RIPK2", "CASP3", "IFIT2", "IFITM2", "PTPN6", "CD69", "CMKLR1", "OAS2", "SOCS1", "SPPL2A", "MYD88", "METTL7B", "TAP1", "CASP1", "IL4R", "PSMA2", "ST8SIA4", "EIF2AK2", "ICAM1", "PNPT1", "PIM1", "XCL1", "RSAD2", "WARS", "DHX58", "NOD1", "IL2RB", "IRF1", "PSME1", "BTG1", "CFH", "VAMP5", "XAF1", "ISG20", "CASP8", "TAPBP", "IL18BP", "SP110", "IRF2", "IFI30", "SLC25A28", "ZNFX1", "PSMB8", "FCGR1A", "CCL2", "PTPN1", "ARID5B", "IFITM3", "MARCH1", "FGL2", "TNFAIP3", "APOL6", "GBP4", "SOCS3", "SECTM1", "SSPN", "STAT4", "IRF8", "JAK2", "IL15RA", "NUP93", "PML", "HLA-DMA", "CD40", "MTHFD2", "IFIH1", "TNFAIP6", "TRAFD1", "CDKN1A", "UBE2L6", "TRIM26", "HLA-DRB1", "NFKB1", "ADAR", "VCAM1", "GBP6", "BATF2", "STAT2", "CD274", "LCP2", "LATS2", "CXCL9", "CIITA", "MT2A", "IRF5", "IL10RA", "NAMPT", "RIPK1", "SLAMF7", "PLA2G4A", "SELP", "KLRK1", "ISOC1", "TRIM21", "GPR18", "PSMB2", "P2RY14", "PTGS2", "FPR1", "TDRD7", "STAT3", "C1S", "PRIC285", "HIF1A", "NCOA3", "DDX58", "PSME2", "EIF4E3", "HLA-A", "IL15", "SERPING1", "VAMP8", "IL6")

ISG.RS <- c("IFI27", "IRF7", "USP18", "BST2", "CXCL10", "DDX60", "HERC6", "HLA-B", "HLA-G", "IFI35", "IFI44", "IFI44L", "IFIT1", "IFIT3", "ISG15", "LGALS3BP", "LY6E", "MX1", "MX2", "OAS3", "OASL", "PLSCR1", "STAT1", "TRIM14", "HSD17B1", "OAS1", "CA2", "CCNA1", "CXCL1", "GALC", "IFI6", "IFITM1", "LAMP3", "MCL1", "ROBO1", "SLC6A15", "THBS1", "TIMP3")


# TIL-B signature (naive vs TIL, Flu vs TIL intersection)
# Wieland A, et al. Defining HPV-specific B cell responses in patients with head and neck cancer. Nature 2021;597:274-278. (online 2020)
TILB_gene <- c("CYS1", "OLR1", "SSTR2", "RBFOX2", "RN7SL627P", "HLA-DRB6", "LAG3", "NLRP4", "CXCL13", "RGS13", "NLRP11", "OSR2", "CCL22", "MUC16", "SHANK2", "TEKT4P2", "LIMA1", "RGS1", "TG", "RASSF6", "CAMP", "CTTN", "CCDC88A", "TNFSF9", "HMSD", "PDGFD", "BCAR1", "PTPRS", "SMPDL3A", "LINC00877", "ATP1B1", "ARHGAP44", "ID2", "RTN2", "SYNE2", "SLCO5A1", "RGS2", "ASPH", "PTPN7", "GRAMD2B", "POLR1E", "TPM4", "KNTC1", "ANGPTL6")


# Galletti G, et al. Two subsets of stem-like CD8+ memory T cell progenitors with distinct fate commitments in humans. Nat Immunol 2020;21:1552-1562.
# https://www.nature.com/articles/s41590-020-0791-5#Sec33
# stem-cell memory T (TSCM) cells self-renew and generate long-lived central memory T (TCM) cells
Tscm.cm_gene <- c('KLRC2','RTKN2','SOX4','LEF1','ZNF683','ACTN1','SPINT2','CD7','XCL1','TXNIP','NUCB2','ARID5B','SELL','HSPB1','MAL','NGFRAP1','KLRD1','C1orf162','PLAC8','LDLRAP1','JUNB','MAP3K1','PABPC1','RP11-291B21.2','TRGC2','TRABD2A','TCF7','FUT7','ELF1','TXK','AREG','IFITM3','CXCR3','FAM65B','LYAR','BEX4','EVL','CD27','ITM2C','PIK3IP1','TRG-AS1','RAB37','TRDC','FCRL3','IFITM1','FCRL6','NOSIP','ITM2A','TRBC1','MBP','CD300A','CCDC69','SYPL1','TMEM123','PPP1R2','NCF1','BEX2','NELL2','HMGN3','DUSP2','PRMT2','SMAP2','PTGDR', # C10; conventional memory T cell intermediate memory and effector
                  'PRKCA','GRASP','DPP4','IRS2','MAP3K1','TRABD2A','GNAQ','SERINC5','NOSIP','DENND5A','SATB1','MAL','VIPR1','MYC','IL6ST','MGAT5','CISH','PIM1','TMEM123','THEM4','CD55','TOB1','CASP10','RUNX2','FAM153A','NT5E','PTGER2','TESPA1','SLC12A6','KLHL6','CHD7','TNFRSF25','ACTN1','RFLNB','CHMP7','ATF7IP2','ACVR1','KLF7','CD226','TXK','SPART','IL11RA','GCSAM','PDE8A','NDFIP1','RAB3GAP1','RNF157','FNIP2','ARHGAP5','RCAN3','IL6R','AQP3','CR1','XBP1','NFKB1','LEF1','PLAG1','ESYT2' # Tscm.cm (vs Tpex)
)

# TEM-like CCR7-GZMKhi (C0); Galletti G, et al. Nat Immunol 2020 (same as above)
#  CCR7- GZMKhi cells expressed cytolytic molecules and generally lacked the memory markers LEF1, CD27, CD28 and CD127
cyto_gene2 <- c('CMC1','RGS1','CD74','GZMK','HLA-DRB1','HLA-DPB1','TIGIT','CD8B','CCL5','MALAT1','CST7','HLA-DPA1','HLA-DQA1','PTPRC','CCL4','HLA-DRA','JUN','JUNB','SH2D1A','HLA-DQB1','GZMM','HLA-DMA','SRSF7','ARID5B','HLA-DRB5','IL10RA','CD27','FCRL6','DNAJB1','DUSP2','LYST','LINC00152','PIK3R1','CD84','IFNG','NR4A2','FYN','PATL2','COTL1','TUBA4A','IFRD1','RHOH','ZNF331','ITM2A','TSPYL2','TRBC1')


# Beltra JC, et al. Developmental relationships of four exhausted CD8+ T cell subsets reveals underlying transcriptional and epigenetic landscape control mechanisms. Immunity 2020;52:825-841. (Table S2)
# Tterm (7) + inhibitory receptors
exh_gene <- c('Pdgfb', 'Aak1', 'Abcb1a', 'Abcb9', 'Abi3', 'Adam8', 'Ap5b1', 'Aplp1', 'Apol10b', 'Apol7e', 'Appl2', 'Arhgap10', 'AW011738', 'BC037034', 'Btbd11', 'Ccrl2', 'Cd160', 'Cd200r1', 'Cd200r4', 'Cd38', 'Cd7', 'Cela1', 'Chst2', 'Cipc', 'Ckb', 'Coro2a', 'Crim1', 'Cysltr2', 'Dtx4', 'Dusp6', 'Enpp5', 'Etv3', 'Fam102a', 'Fam78a', 'Gadd45b', 'Gdf11', 'Glcci1', 'Glrx', 'Gm12253', 'Hic1', 'Il12rb1', 'Il21r', 'Itih5', 'Kcna3', 'Klk8', 'L3mbtl3', 'Litaf', 'Lpar6', 'Nedd9', 'Padi2', 'Pdcd1', 'Pde4b', 'Pde7a', 'Pelo', 'Pglyrp2', 'Pik3ap1', 'Ppp1r16b', 'Ptger4', 'Ret', 'Rgs1', 'Rgs3', 'Rundc3b', 'Sepn1', 'Sh2d2a', 'Slc16a6', 'Slc22a15', 'Slc39a6', 'Snai3', 'Soat2', 'Srbd1', 'St14', 'Syne2', 'Tcrg-C1', 'Tmem135', 'Tmem140', 'Tnfsf10', 'Wbp1', 'Zbtb38', "LAG3", "HAVCR2", "PDCD1", "ENTPD1", "CTLA4", "TOX", "TOX2", "TIGIT", "CD244") %>% toupper(.) %>% unique() %>% as.vector()

# Tpex; Beltra JC, et al. Immunity 2020 (same as above)
pex_gene <- c('Crtam', 'Cxcr5', 'Ifi16', 'Mnda', 'Tnfsf4', 'Tcf7') %>% toupper()


# Tumor-reactive CD8 T cell signature (melanoma)
# Oliveira G, et al. Phenotype, specificity and avidity of antitumour CD8+ T cells in melanoma. Nature 2021;596:119-125.
TR_gene1 <- c("GNLY", "KRT86", "CXCR6", "PRF1", "CHST12", "SLA2", "ACP5", "CD7", "ISG15", "ISG20", "OASL", "ANXA1", "GZMB", "MX1", "RPS26", "PBXIP1", "LEPROTL1", "RSAD2", "SELPLG", "ARHGAP9", "ATP1B3", "CCL5", "RARRES3", "S100A11", "ANXA2", "IFIT1", "IFIT3", "CREM", "CLIC3", "CRIP1", "TAGLN2", "SRGN", "SLC7A5", "GALNT2", "RALA", "HMOX1", "APOBEC3G", "PLSCR1", "LITAF", "FTH1", "LAYN", "TGFB1", "TNFSF10", "GGA2", "IFITM3","HSPA6", "HSPA1B", "DNAJB1", "CCL4", "CCL4L2", "BAG3", "HSPA1A", "DUSP1", "FOS", "SERPINH1", "NR4A2", "IFNG", "NFKBIA", "ZFAND2A", "GADD45B", "CCL3", "JUNB", "NR4A1", "DNAJA1", "PPP1R15A", "HSPB1", "HSPH1", "KLF6", "RGS1", "ATF3", "DEDD2", "RGS2", "XCL1", "HSP90AA1", "HSPD1", "CRTAM", "HSP90AB1", "ZFP36L1", "CMC1", "GZMK", "HSPE1", "FOSB", "CACYBP", "XCL2", "JUN", "CD69", "MALAT1", "TAGAP", "CHORDC1", "PMAIP1", "GLA", "MRPL18", "ID2", "HSPA8", "RSRP1","CCR7", "CXCL13", "XCL1", "XCL2", "NMB", "CRTAM", "GNG4", "LANCL2", "HLA-DRA", "DUT", "REL", "TYMS", "HLA-DPA1", "GSTP1", "HLA-DQA1", "PKM",  "CD74", "NME2", "RPLP0", "TCF7", "HLA-DRB1", "TUBB", "STMN1", "HIST1H4C", "TUBA1B", "MKI67", "HMGN2", "TYMS", "UBE2C", "HIST1H1B", "CENPF", "HIST1H1E", "HMGB1", "HIST1H1D", "CKS1B","H2AFZ", "KIF22", "HMGB2", "PKMYT1", "H2AFV", "CDT1", "PCNA", "DUT", "NUSAP1", "CKS2", "ASPM", "TPX2", "BIRC5", "TOP2A", "CENPW", "ZWINT","ANP32B", "CDCA5", "CDK1", "CKAP2L", "SMC4", "TK1", "PTTG1", "ASF1B", "CENPM", "CDKN3", "MAD2L1", "CDC20", "DNAJC9", "KIFC1", "CCNA2", "TMEM106C", "CENPN", "CENPU", "HIST2H2AC", "ACTB", "IL7R", "GPR183", "LMNA", "ANXA1", "VIM", "FTH1", "MBP", "CD55", "TGFB1", "CCR7", "S100A10", "KLRB1", "MARCKSL1", "ANXA2", "GNLY", "TCF7", "KLF3", "FOS", "RPS12", "PRNP", "KLF2", "PABPC1", "ABLIM1", "RPS16", "KLRD1", "GZMM", "CDKN1A", "RPL3", "SLC2A3", "GLUL", "RPL13", "RPL32", "RPS29", "CD44", "RPS27", "TPT1") %>% unique() %>% as.vector()


# Anti-PD1 neoantigen-reactive signature (metastatic melanoma)
# Lowery FJ, et al. Molecular signatures of antitumor neoantigen-reactive T cells from metastatic human cancers. Science 2022;375:877-884.
TR_gene3 <- c("CXCL13", "ADGRG1", "PDE7B", "PDCD1", "ETV7", "HMOX1", "CCDC50", "AHI1", "CD4", "TSHZ2", "ENTPD1", "TOX", "GADD45G", "IFNG", "GNG4", "MAF", "KIAA0319L", "LAG3", "MYL6B", "NMB", "ZNRF1", "DUSP4", "JAK3", "YPEL2", "ELMO1", "TNFRSF18", "ORMDL3", "IGFLR1", "F2R", "TIGIT", "MIS18BP1", "CCDC167", "CCL4", "ITGAL", "NIN", "LIMS1", "PIM2", "SECISBP2L", "TMEM173", "GZMK", "PTMS", "ID2", "LYST", "RILPL2", "ZKSCAN1", "COTL1", "BHLHE40", "TLK1", "GALM", "ITGB2", "PRDM1", "NAP1L4", "LBH", "PAG1", "C16orf87", "CPT1A", "NUTM2B-AS1", "ITGA4", "CXCR6", "DDA1", "MAL", "VOPP1", "ZBTB38", "SUPT3H", "CYTH1", "RNF213", "CNST", "STT3B", "CTSB", "NME3", "STAT1", "LUC7L3", "BATF", "CARHSP1", "CCNL1", "PDE4D", "CD82", "AKAP13", "HNRNPR", "RPS4Y1", "HNRNPLL", "SFXN1", "ADI1", "SYNGR2", "ATRX", "ILKAP", "PDE3B", "PPM1G", "PPP1CC", "CORO1B",  "ITM2A", "SPOCK2", "VPS25", "PHYKPL", "NR3C1", "ARGLU1") %>% unique() %>% as.vector()

