# MicroArrayPipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3770853.svg)](https://doi.org/10.5281/zenodo.3770853)


R Shiny app to process Affymetrix raw data

Main workflow:
- Cel file or GEO input 
- Reading raw data and metadata
- Plotting raw data QC plots
- RMA Normalization 
- Plotting Post-normalization QC plots
- PCA & clustering
- Limma call for Differentially Expressed Genes DEG 
- Pathway enrichment analysis 
