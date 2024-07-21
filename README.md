# Xenium_Line61_spatial

Code used to analyze data and generate figures for Horan-Portelance et al., 2024, *Single-cell spatial transcriptomics reveals molecular patterns of selective neuronal vulnerability to α-synuclein pathology in a transgenic mouse model of PD/DLB*. DOI: 

We performed Xenium, an imaging-based spatial transcriptomics assay, in Line 61 mice, a mouse model of α-synucleinopathy. There are 8 total samples, 4 non-transgenic (non-tg or nTg) and 4 transgenic (α-syn-tg or Tg), with 2 males and 2 females in each genotype. We ran Xenium with the base mouse brain panel and a 100-plex custom panel of probes, all of which can be found in Supplementary Table 1 of the paper. 

We also performed post-Xenium immunofluorescence in the same sections used in the assay, staining for phosphorylated α-synuclein. These images were overlaid with Xenium cells and cells were assigned as pSyn+ or pSyn- for downstream analysis. These images can be found in our Zenodo repository. 

Raw Xenium data from this experiment can be found on GEO.

We also generated a Zenodo repository which holds the images as well as processed Seurat objects.

The code in this repository is everything which can analyze the data as well as generate all the figures for the paper, both main and supplementary. A full list of figures and where you can find the code for them is below: 

Folder | Script | Figures
--- | --- | ---
2_CTX_analysis | 4_ctx_basic_plots | 1a, 1b, 1c, 1d, 1e, S4a, S4b
2_CTX_analysis | 5_ctx_pSyn_visualization | 2a, 2b
2_CTX_analysis | 6_ctx_gene_exp_analysis | 2c, 2d, 2e
2_CTX_analysis | 7_ctx_gene_expression_linregs | 2f
2_CTX_analysis | 10_ctx_DE_GLM_visualization | 6b, 6c, 6d, 6e
2_CTX_analysis | 11_ctx_pSyn_percentiles_analysis | 5a, 5b, S6a
2_CTX_analysis | 13_ctx_pSyn_percentiles_visualization | 5c
2_CTX_analysis | 15_ABC_Xenium_comparison | S3a, S3b, S3c
2_CTX_analysis | 16_ctx_DE_euler_plots | S7a
2_CTX_analysis | 17_ctx_QC_metrics | S2a, S2b, S2c, S2d
3_HIP_analysis | 4_hip_basic_plots | 3a, 3b, 3c, 3d
3_HIP_analysis | 5_hip_pSyn_visualization | 4a, 4b
3_HIP_analysis | 6_hip_gene_expression_analysis | 4c, 4d
3_HIP_analysis | 7_hip_QC_metrics | S5a, S5b, S5c, S5d
