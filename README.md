# MicroArrayPipeline

R Shiny app to process Affymetrix raw data

Main workflow:
- Reading raw data and metadata
- Plotting raw data QC plots
- RMA Normalization 
- Plotting Post-normalization QC plots
- PCA & clustering
- Limma call for Differentially Expressed Genes DEG 
- KEGG/GO Enrichment analysis 
