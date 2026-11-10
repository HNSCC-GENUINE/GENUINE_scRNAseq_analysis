# GENUINE_scRNAseq_analysis
Analysis codes for single-cell RNA sequencing and immunological profiling


# 1. Human Data Analysis Using R

## Data availability

The processed data objects (`0.data/`) are not included in this repository.
During peer review, raw data are accessible via GSE353363 (private reviewer token required, given in manuscript), and processed data will be provided upon reviewer request. 


## Data files 

| File | Description |
|---|---|
| `HNSCC_whole_adata_public.rds` | Seurat object, all cells (106,122 cells) |
| `HNSCC_B_subset_public.rds` | Seurat object, B/plasma cells (23,669 cells) |
| `HNSCC_B_monocle3_cds_public.rds` | Monocle 3 `cell_data_set` of the B subset (same cells as above), with trajectory |
| `HNSCC_CD8_subset_public.rds` | Seurat object, CD8 T / innate-like T / NK cells (22,375 cells) |
| `HNSCC_BCR_contigs_public.tsv` | BCR contigs (AIRR-style columns, SHM, clone IDs) for cells in the B subset |
| `HNSCC_CD8_TCR_immunarch_public.rds` | TCR data as an *immunarch* object (`$data`, `$meta`); 7 CD8 T cell subsets only |
| `HNSCC_CD8_TCR_cells_public.tsv` | Same TCR data as a single cell-level table |
| `HNSCC_TLS_lunit_public.tsv` | Patient-level tertiary lymphoid structure (TLS) metrics from Lunit SCOPE IO image analysis |
| `Fig5G_network_layout.RData` | Saved node coordinates (`layout_saved`) for the Fig 5G network graph |

All Seurat objects contain the `RNA` assay only (counts and log-normalized data) and the `pca_harmony` and `umap` reductions. They were saved in SeuratObject v4.1.3 format (Assay v3). They load directly in Seurat v4. In Seurat v5 they can be updated with `UpdateSeuratObject()`.

## Metadata

| Column | Values | Description |
|---|---|---|
| `Patient` | `P001`–`P015` (11 patients) | De-identified patient ID |
| `Tissue` | `pre_tumor`, `post_tumor`, `post_LN` | Pre-treatment tumor biopsy, post-treatment tumor, post-treatment lymph node |
| `Response` | `SD`, `mPR`, `pCR` | Pathologic response |
| `Response2` | `LR`, `ER` | Limited responder (SD) / effective responder (mPR, pCR) |
| `level1` | 13 major cell types | Major annotation (whole object) |
| `level2` | B: 10 subsets; CD8: 10 subsets | Subset annotation (B and CD8 objects) |
| `level2_robust` | e.g. `TR` | CD8 object only. `TR` = tumor-reactive-like subsets (CD8_EM_GZMK+, CD8_Tpex, CD8_Tex, CD8_ISG+, Proliferating_CD8) |
| `Immune` | `immune`, `non_immune` | Whole object only |
| `pseudotime`, `pathway` | numeric; `Path1`, `Path2`, `Mixed`, `None` | B object and cds only; Monocle 3 trajectory |

Factor levels define plotting order. Colors are stored with the objects:
`obj@misc$cols$<column>` (Seurat) and `metadata(cds)$cols$<column>` (Monocle 3).

## Linking cells across files

- All files are linked by the **Seurat cell name** (`colnames(obj)`).
- In the BCR and TCR tables, `cell_id` is the cell name of the matching cell in the B or CD8 object.
- Samples are identified by `Patient` + `Tissue`. The numeric suffix of the cell names (e.g. `_1_1`) comes from object merging and is **not** a sample ID.
- For each patient, the post-treatment tumor and lymph node samples were multiplexed with cell hashing (CMO) and sequenced in one library.

## Code execution notes

Run `1.rscript/R_human_preset_code.R`, then each figure script **from top to bottom**.

Order dependencies within the scripts:
- `R_human_fig3_S7_S8_S9.R`: Fig 3D creates the `score` object used by Fig 3E–H and Supplementary Fig 9A. Supplementary Fig 9B recomputes `score` on Path1 cells, which Supplementary Fig 9C, 9E and 9G then use.
- `R_human_fig4_S11.R` / `R_human_fig5_S12.R`: main panels use `CD8_7`, the 7 CD8 T cell subsets created at the top of the script. NK, MAIT and rdT cells appear only in Supplementary Fig 11A–B. Fig 4F uses the DEG result from Fig 4E.

## Software

R 4.1.3. Main packages: Seurat 4.4.0, monocle3 1.3.4, fgsea 1.20.0, msigdbr 7.5.1, immunarch 0.9.1, ComplexHeatmap 2.10.0, circlize 0.4.17, igraph 2.2.1, vegan 2.7-2, ggplot2 4.0.1. Full versions are in `renv.lock` and `261111_sessionInfo.txt`.

## Gene sets

- GO Biological Process gene sets from **msigdbr v7.5.1**. Because the MSigDB release differs from the one used for the manuscript figures, GSEA NES values may differ slightly; the direction of the main pathways is unchanged. Pathways shown in the manuscript were selected from the GSEA output.
- EF/GC antibody-secreting cell DEGs: Ma J, et al. *Science* 2024;384:eadj4857 (Table "5. DEGs of EF and GC ASCs").
- Other signatures (TIL-B, stem-like/memory, cytotoxic, exhaustion, Tpex, tumor-reactive, IFNG/ISG) are listed in `R_human_preset_code.R` with their sources.

## Notes on patient exclusion

- **TLS association analyses** (Fig 3L, Fig 4I, Fig 5E): patients without evaluable paired slides were excluded. The number of patients is given in each figure legend.
- **Proportion change analyses** (Supplementary Fig 7D, Fig 5B): patients without paired pre- and post-treatment samples for the compartment analyzed were excluded.
- **Paired signature scores** (Fig 4G, Supplementary Fig 11D): n = 9; only patient × subset pairs with both time points are shown.
