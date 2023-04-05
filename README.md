# TCGADEPMAP
This repo contains code for gene expression alignment presented in TCGADEPMAP – Mapping Translational Dependencies and Synthetic Lethalities within The Cancer Genome Atlas. Part of the codes are incorportated from [Celligner](https://github.com/broadinstitute/Celligner_ms)

## Manuscript
The preprint version of the manuscript is available at https://www.biorxiv.org/content/10.1101/2022.03.24.485544v2.

## Manuscript Supplemental Data
The predicted gene essentiality scores for TCGADEPMAP(Table S5), GTEXDEPMAP(Table S11) and PDXEDEPMAP (Table S10) are available at [figshare](https://figshare.com/projects/TCGADEPMAP_Mapping_Translational_Dependencies_and_Synthetic_Lethalities_within_The_Cancer_Genome_Atlas/130193)

## Data

The TCGA expression data is available from the [treehouse dataset](https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

The DEPMAP cell line expression data is downloaded from [DEPMAP portal](https://depmap.org/portal/).

The metadata for the analysis is downloaded from https://figshare.com/articles/Celligner_data/11965269.

### Dependencies:

Here is the list of the dependencies in R packages which are available in Bioconductor.
'here', 'tidyverse', 'reshape2', 'plyr', 'data.table', 'Seurat', 'pheatmap', 'pdist', 'gridExtra', 'ggpubr', 'grDevices', 'RColorBrewer', 'FNN', 'ggrepel', 'ggridges', 'irlba', 'viridis', 'limma', 'edgeR', 'batchelor', 'BiocParallel', 'sva', 'GSEABase', 'piano', 'fgsea', 'preprocessCore'
