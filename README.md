# MAAPster

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3770853.svg)](https://doi.org/10.5281/zenodo.3770853) [![GitHub releases](https://img.shields.io/github/release/CCBR/MicroArrayPipeline)](https://github.com/CCBR/MicroArrayPipeline/releases) [![GitHub issues](https://img.shields.io/github/issues/CCBR/MicroArrayPipeline)](https://github.com/CCBR/MicroArrayPipeline/issues) [![GitHub license](https://img.shields.io/github/license/CCBR/MicroArrayPipeline)](https://github.com/CCBR/MicroArrayPipeline/blob/master/LICENSE)

MicroArray Analysis Pipeline, also known as `MAAPster`, is a comprehensive Shiny application and R package that performs transcriptome analysis of human or mouse Affymetrix gene expression data. Please see the [Platform Support](#Platform-Support) section for a list of all currently supported Affymetrix chips.

Samples may be uploaded locally or accessed from published data by entering a GEO Series identification number. MAAPster can analyze data from a single experiment that includes multiple samples, and the user may investigate multiple contrasts between groups of samples. Raw CEL files are analyzed using several Bioconductor packages in R, including limma and oligo. Output includes array probe quality control plots, sample quality control plots such as 3D PCA, differential gene expression analysis, gene expression heatmaps, pathway analysis and single sample GSEA analysis.



The main workflow consists of the following steps:
- [Reading in raw data (local or GEO) and metadata](#Input-Data)
- [Generating raw data QC plots](#Array-Probe-Quality-Control)
- [RMA or Cyclic Loess Normalization](#Normalization)
- [Generating Post-normalization QC plots](#Post-normalization-QC-plots)
- [PCA and clustering](#Sample-Quality-Control)
- [Limma call for identifying DE genes](#Differential-Gene-Expression) 
- [Pathway enrichment analysis](#Pathway-Analysis)

## Input Data

There are two input options. Raw data can be uploaded locally or from the Gene Expression Omnibus (GEO): 

**Local File Upload:** Analyze expression data by uploading CEL files from a local computer. Data must be expression data in CEL file format, and all data must be from the same type of Affymetrix chip.

**GEO:** Analyze expression data from the Gene Expression Omnibus ([GEO](https://www.ncbi.nlm.nih.gov/geo/)) by providing the relevant series ID from the experiment, as noted below. Enter the entire ID, starting with GSE. Samples from different chips cannot be analyzed together. If the series ID contains samples from multiple chips, choose the type of chip from the drop-down menu. Samples analyzed on each chip will populate as the chip is chosen

Once the files are uploaded, groups are assigned for each sample. At least 2 groups must be created. The control group should be entered last. Each group must have at least two replicates.

***Optional:*** If batch effects are present in the dataset, assign batch names as shown. Batch effects will be corrected for differential gene expression analysis by incorporating batches into the linear model. For plots based on normalized gene expression and before ssGSEA analysis, normalized gene expression is corrected for batch effects with the ComBat function from the sva package *(Leek, Johnson, Parker, Jaffe, & Storey, 2012)*. 

## Array Probe Quality Control

Two plots, NUSE and RLE, display quality control metrics for the microarray probes. Plots are created with the oligo package:

**NUSE:** The Normalized Unscaled Standard Error (NUSE) plot shows standard error estimates for each sample. These estimates are normalized to 1 across all samples, so a sample with a significantly higher value may be of lower quality *(Carvalho & Irizarry, 2010)*.

**RLE:** The Relative Log Expression (RLE) plot compares the expression level of one probeset to the median expression of that probeset across all samples.

## Normalization

**RMA:** Each sample is normalized to correct for technical artifacts using the Robust Multichip Averaging (RMA) method described in *(Irizarry et al., 2003)*. 

**Cyclic loess:** This additional round of normalization is also available and can be performed if RMA is insufficient *(Ballman, Grill, Oberg, & Therneau, 2004; Bolstad, Irizarry, Astrand, & Speed, 2003)*. Only selected samples are used for normalization. 

## Post-normalization QC plots

Pre-normalization and post-normalization histograms, MA plots and box plots display relative expression (represented by probe intensity) before and after normalization. Data for the plots are generated with the oligo package *(Carvalho & Irizarry, 2010)*.

**Histogram:** The histogram displays the density of probes at each log-intensity. Each sample is represented by a curved line. Before normalization, samples should be distributed randomly. An outlier sample may be due to technical problems, and separate groupings of samples may indicate batch effects. Either of these potential issues should be investigated further before downstream analyses. After normalization, all samples lines should follow the same curve. Again, if samples are clustered into groups, batch effects may be present.

**MA plot:** Log-ratios (M) are plotted against the average log-intensities (A). Plots are generated that compare each sample to the median of all samples. After RMA, the red sample line should be relatively equal to the blue median line.

**Boxplot:** Displays the distribution of log-intensities for each sample. After RMA, distributions across all samples should be approximately the same.

## Sample Quality Control

The sample similarity heatmap and 3D PCA provide information about the quality of replicates in the groups.

**PCA:** In general, samples in the same group should cluster together, and groups of samples should cluster separately from other groups. If batch effects are present, re-run the analysis and assign samples to respective batches for correction. For example, the dataset below has batch effects. The first principal component, PC1, separates replicates 1 and 2 (rep1 and rep2) from replicates 3 and 4 (rep3 and rep4):

## Differential Gene Expression

The DEG-Enrichments Results tab displays differentially expressed genes between groups. MAAPster runs analysis in the background with the limma package. Linear modeling is performed using limma’s lmFit function, and differential gene expression is determined using the contrasts.fit and eBayes functions *(Ritchie et al, 2015)*. Using the toptable function, false discovery rates are calculated to adjust p-values for multiple testing *(Ritchie et al., 2015)*.

## Pathway Analysis

**l2P:** The top 500 significantly upregulated and top 500 significantly downregulated genes are extracted for pathway analysis (significance is determined as unadjusted p-value < 0.05), Pathway analysis is performed with CCBR’s [l2p R package](https://github.com/CCBR/l2p).

**Single-sample GSEA:** ssGSEA is performed using the gsva package as previously described *(Hanzelmann et al., 2013)*. Pathway enrichment scores are then analyzed to determine fold changes and p-values between groups of samples, similar to the differential gene expression analysis described above. Differential pathway enrichment is ranked by p-value, and the top 50 pathways are displayed in a heatmap. Human gene set modules were downloaded from the BROAD Institute’s MSigDB, and mouse gene set modules were downloaded from the Gene Set Knowledgebase.

## Platform Support

The pipeline supports the following Affymetrix chips, more maybe added upon request:

**Human:** `Human Genome U133 Plus 2.0 Array`, `GeneChip™ Human Genome U133A Array`, `GeneChip™ Human Genome U133A 2.0 Array`, `GeneChip™ Human Genome U133B Array`, `GeneChip™ Human Gene 1.0 ST Array (v1 and v2)`, `GeneChip™ Human Gene 1.1 ST Array Version 1`, `GeneChip™ Human Gene 2.0 ST Array`, `GeneChip™ Human Gene 2.1 ST`, `Clariom™ S Assay HT, human`, `Clariom™ S Assay, human`, `Human Genome U219 Array`, `GeneChip Human Genome U95 Version 2`, `GeneChip™ Human Transcriptome Array 2.0`, `Human Exon 1.0 ST Array`, `HT Human Genome U133 Array`, `Clariom™ D Assay, human`

**Mouse:** `GeneChip™ Mouse Gene 1.0 ST Array`, `GeneChip™ Mouse Gene 1.1 ST Array`, `GeneChip™ Mouse Gene 2.0 ST Array`, `Clariom™ S Assay HT, mouse`, `Clariom™ S Assay, mouse`, `GeneChip™ Mouse Genome 430 2.0 Array`, `GeneChip™ Mouse Genome 430A 2.0 Array`, `GeneChip® Mouse Expression Set 430`, `GeneChip® Murine Genome U74v2 Set`, `Mouse Exon 1.0 ST Array`

## References

<sup>Ballman, K. V., Grill, D. E., Oberg, A. L., & Therneau, T. M. (2004). Faster cyclic loess: normalizing RNA arrays via linear models. Bioinformatics, 20(16), 2778-2786. doi:10.1093/bioinformatics/bth327</sup>

<sup>Bolstad, B. M., Irizarry, R. A., Astrand, M., & Speed, T. P. (2003). A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics, 19(2), 185-193.</sup>

<sup>Carvalho, B. S., & Irizarry, R. A. (2010). A framework for oligonucleotide microarray preprocessing. Bioinformatics, 26(19), 2363-2367. doi:10.1093/bioinformatics/btq431</sup>

<sup>Hanzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics, 14, 7. doi:10.1186/1471-2105-14-7</sup>

<sup>Irizarry, R. A., Hobbs, B., Collin, F., Beazer-Barclay, Y. D., Antonellis, K. J., Scherf, U., & Speed, T. P. (2003). Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Biostatistics, 4(2), 249-264. doi:10.1093/biostatistics/4.2.249</sup>

<sup>Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E., & Storey, J. D. (2012). The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics, 28(6), 882-883. doi:10.1093/bioinformatics/bts034</sup>

<sup>Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43(7), e47. doi:10.1093/nar/gkv007</sup>

<hr>

<p align="center">
	<a href="#maapster">Back to Top</a>
</p>
