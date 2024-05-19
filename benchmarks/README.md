### ibaqpy - IBAQ Python benchmarks and tests 

This repository contains the benchmarks and tests for the IBAQ Python package (ibaqpy). The aim of this document is to provide a detailed description and documentation of different benchmarks on different datasets. An original benchmark of the package was published in [Proteomics in 2023](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.202300188), but it was only aiming to demonstrate if the package and quantms could be used to analyze large-scale TMT and LFQ datasets. 

In these series of benchmarks, we aim to test the complete packages including the different methods for protein quantification, feature and peptide selection and normalization. Also, we aim to benchmark different metrics and parameters to remove low-quality features and peptides; and finally, we aim to test the performance of the package in different datasets. These are the following questions we aim to answer with these benchmarks:

- How does the package perform in different datasets?
- What method brings the less Coefficient of Variation (CV) in the quantification results across samples? 
- What method brings the best correlation between the ibaq values in TMT vs. LFQ?
- What method brings less missing values across samples?
- What method brings less technical variability across samples? 

last created: 2024-05-20 

### Dataset PXD007683

This dataset from Gygi's lab was originally published in JPR as [Proteome-Wide Evaluation of Two Common Protein Quantification Methods
](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00016). The dataset is available at [PXD007683](https://www.ebi.ac.uk/pride/archive/projects/PXD007683).
The authors test the ability of two common methods, a tandem mass tagging (TMT) method and a label-free quantitation method (LFQ), to achieve comprehensive quantitative coverage by benchmarking their capacity to measure 3 different levels of change (3-, 2-, and 1.5-fold) across an entire data set. The authors reported Both methods achieved comparably accurate estimates for all 3-fold-changes. 

#### Coefficient of Variation (CV)

#### Correlation between TMT and LFQ samples

#### Variability of specific proteins across samples. 

#### Missing values across samples

#### Fold-change detection 3-, 2-, and 1.5-fold. 

