# Xenium_Line61_spatial

Code used to analyze data and generate figures for Horan-Portelance et al., 2024, "Single-cell spatial transcriptomics reveals molecular patterns of selective neuronal vulnerability to α-synuclein pathology in a transgenic mouse model of PD/DLB". DOI: 

We performed Xenium, an imaging-based spatial transcriptomics assay, in Line 61 mice, a mouse model of α-synucleinopathy. There are 8 total samples, 4 non-transgenic (non-tg or nTg) and 4 transgenic (α-syn-tg or Tg), with 2 males and 2 females in each genotype. We ran Xenium with the base mouse brain panel and a 100-plex custom panel of probes, all of which can be found in Supplementary Table 1 of the paper. 

We also performed post-Xenium immunofluorescence in the same sections used in the assay, staining for phosphorylated α-synuclein. These images were overlaid with Xenium cells and cells were assigned as pSyn+ or pSyn- for downstream analysis. These images can be found in our Zenodo repository. 

Raw Xenium data from this experiment can be found on GEO: URL

We also generated a Zenodo repository which holds the images as well as processed Seurat objects: URL

The code in this repository is everything which can analyze the data as well as generate all the figures for the paper, both main and supplementary. A full list of figures and where you can find the code for them is below: 

Figure | Folder | Script
--- | --- | ---
Fig. 1a | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 1b | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 1c | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 1d | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2a | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2b | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2c | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2d | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2e | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 2f | 2_CTX_analysis | 4_ctx_basic_plots
Fig. 3a | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 3b | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 3c | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 3d | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 4a | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 4b | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 4c | 3_HIP_analysis | 4_ctx_basic_plots
Fig. 4d | 3_HIP_analysis | 4_ctx_basic_plots
