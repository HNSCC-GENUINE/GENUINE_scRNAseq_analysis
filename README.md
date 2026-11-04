# GENUINE_scRNAseq_analysis
Analysis codes for single-cell RNA sequencing and immunological profiling


## External reference data

`261103_human_fig_preset.Rmd` loads one external file, used to build the GC/EF
antibody-secreting cell (ASC) gene signatures (Fig. 3D–E, Supplementary Fig. 9B):

| File | Source | Sheet |
|---|---|---|
| `sci2024_GCEF_DEG_list.xlsx` | Supplementary table of Ma J, et al. *A blueprint for tumor-infiltrating B cells across human cancers.* Science 2024;384:eadj4857. https://doi.org/10.1126/science.adj4857 | `5. DEGs of EF and GC ASCs` |

Download the supplementary Excel file from the publisher's website, rename it to
`sci2024_GCEF_DEG_list.xlsx`, and place it in the working directory.
The sheet is read from row 2 (`startRow = 2`), and the columns `Gene`, `Group`,
`LogFC` and `P_adj` are used. Genes are selected as follows:
- GC signature: `Group == "GC"`, `LogFC <= -1`, `P_adj <= 0.05`
- EF signature: `Group == "EF"`, `LogFC >= 1`, `P_adj <= 0.05`
