# ibaqpy

[![Python application](https://github.com/bigbio/ibaqpy/actions/workflows/python-app.yml/badge.svg)](https://github.com/bigbio/ibaqpy/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/ibaqpy/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/ibaqpy/actions/workflows/python-publish.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/6a1961c7d57c4225b4891f73d58cac6b)](https://app.codacy.com/gh/bigbio/ibaqpy/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![PyPI version](https://badge.fury.io/py/ibaqpy.svg)](https://badge.fury.io/py/ibaqpy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ibaqpy)

iBAQ (Intensity-Based Absolute Quantification) determines the abundance of a protein by dividing the total precursor intensities by the number of theoretically observable peptides of the protein [manuscript here](https://pubmed.ncbi.nlm.nih.gov/16219938/). ibaqpy is a Python package that computes iBAQ values starting from a feature parquet from [quantmsio](https://github.com/bigbio/quantms.io) and a [SDRF](https://github.com/bigbio/proteomics-sample-metadata) file. In addition, the package computes other ibaq values including rIBAQ, log2, and ppb.

ibaqpy also allows computing the TPA value (Total Protein Approach), protein copy number, and protein concentration. TPA is determined by summing peptide intensities of each protein and then dividing by the molecular mass to determine the relative concentration of each protein. By using [ProteomicRuler](https://www.sciencedirect.com/science/article/pii/S1535947620337749), it is possible to calculate the protein copy number and absolute concentration. The OpenMS tool was used to calculate the theoretical molecular mass of each protein. Similar to the calculation of IBAQ, the TPA value of protein-group was the sum of its intensity divided by the sum of the theoretical molecular mass.

The protein copy calculation follows the following formula:

```
protein copies per cell = protein MS-signal *  (avogadro / molecular mass) * (DNA mass / histone MS-signal)
```

For cellular protein copy number calculation, the uniprot accession of histones was obtained from species first, and the molecular mass of DNA was calculated. Then the dataframe was grouped according to different conditions, and the copy number, molar number and mass of proteins were calculated. In the calculation of protein concentration, the volume is calculated according to the cell protein concentration first, and then the protein mass is divided by the volume to calculate the intracellular protein concentration.

### Overview of ibaq-base values computation

As mentioned before, ibaq values are calculated by dividing the total precursor intensities by the number of theoretically observable peptides of the protein. We use the following steps to calculate the iBAQ values:

- _Observable peptides_, the protein sequence is digested in silico using a specific enzyme. The current version of this tool uses OpenMS method to load fasta file, and use [ProteaseDigestion](https://openms.de/current_doxygen/html/classOpenMS_1_1ProteaseDigestion.html) to enzyme digestion of protein sequences, and finally get the theoretical peptide number of each protein.

- _Total precursor intensities_, the total intensity of a protein is calculated by summing the intensity of all peptides that belong to the protein. The intensity values are obtained from the feature parquet file in [quantms.io](https://github.com/bigbio/quantms.io).

> Note: If protein-group exists in the peptide intensity dataframe, the intensity of all proteins in the protein-group is summed based on the above steps, and then divided by the number of proteins in the protein-group.

### Other values  

- `IbaqNorm` - normalize the ibaq values using the total ibaq of the sample `ibaq / sum(ibaq)`, the sum is applied for proteins in the same _sample + condition_.

- `IbaqLog`  - The ibaq log is calculated as `10 + log10(IbaqNorm)`. This normalized ibaq value was developed [by ProteomicsDB Team](https://academic.oup.com/nar/article/46/D1/D1271/4584631).

- `IbaqPpb` - The resulted IbaqNorm is multiplied by 100M `IbaqNorm * 100'000'000`. This method was developed originally [by PRIDE Team](https://www.nature.com/articles/s41597-021-00890-2).

- `TPA` - TPA value is calculated as `NormIntensity / MolecularWeight`

- `CopyNumber` - Protein copy number is calculated by a proteomic ruler approach.

- `Concentration[nM]` - Protein concentration is calculated using the total weight and a provided concentration per cell (cpc).

### From quantms to Ibaq values

![Ibaq](./data/ibaqpy.drawio.png "IBAQ")

The output of quantms is converted into quantms.io feature file. quantms.io provides a unified format for processing report files, including peptide intensity information. In quantms.io, you can use the `convert-ibaq` command, providing a **feature file** and an **SDRF file**, to inject experimental information into the feature file, generating an ibaqpy use case.

```asciidoc
>$ quantmsioc convert-feature --sdrf_file PXD004452-Hella-trypsin.sdrf.tsv --msstats_file PXD004452-Hella-trypsin.sdrf_openms_design_msstats_in.csv --mztab_file PXD004452-Hella-trypsin.sdrf_openms_design_openms.mzTab --file_num 30 --output_folder res --duckdb_max_memory 64GB --output_prefix_file PXD004452
>$ quantmsioc convert-ibaq --feature_file res/PXD004452-6c224f5a-7c1f-46f9-9dae-1541baeef8fe.feature.parquet --sdrf_file PXD004452-Hella-trypsin.sdrf.tsv --output_folder ibaq --output_prefix_file PXD004452
```

A feature in quantms.io is the combination of the following columns: 

- `ProteinName`: Protein name
- `Peptidoform`: Peptide sequence including post-translation modifications `(e.g. .(Acetyl)ASPDWGYDDKN(Deamidated)GPEQWSK)`
- `PEPTIDE_CANONICAL`: Canonical peptide sequence
- `PrecursorCharge`: Precursor charge
- `Channel`: Lable channel
- `Condition`: Condition label `(e.g. heart)`
- `BioReplicate`: Biological replicate index `(e.g. 1)`
- `Run`: Run index `(e.g. 1)`
- `Fraction`: Fraction index `(e.g. 1)`
- `Intensity`: Peptide intensity
- `Reference`: reference file
- `SampleID`: Sample ID `(e.g. PXD003947-Sample-3)`


In summary, each feature is the unique combination of a peptide sequence including modifications (peptidoform), precursor charge state, condition, biological replicate, run, fraction, reference_file_name, sample_accession, and a given intensity. In order to go from these features into protein ibaq values, the package does the following: 

#### Data preprocessing

In this section`features2peptides`, ibaqpy will do: 
- Parse the identifier of proteins and retain only unique peptides.
- Remove lines where intensity or study condition is empty: This could happen in the following cases: 
  - The DIA pipeline sometimes for some features releases intensities with value 0.
  - The quantms.io do not contain feature information for some conditions. This extreme case could happen when not ID/Quant was found for a given condition during the analysis.
- Filter peptides with less amino acids than min_aa.
- Low-confidence proteins were removed according to the threshold of unique peptides: We use a thershold of 2 unique peptides to consider a protein as a high-confidence protein. This parameter is applied if not specified by the user, and the default value is 2. If users want to change this threshold, they can use the `--min_unique` parameter.
- Filter decoy, contaminants, entrapment: Proteins with the following prefix are removed by default: `DECOY, CONTAMINANT, ENTRAPMENT` could be removed, by default, the filter is not applied. If users want to remove these proteins, they can use the `--remove_decoy_contaminants` parameter.
- Filter user-specified proteins: The user can provide a list of protein identifiers to remove from the analysis using the `--remove_ids` parameter. The remove ids parameters will remove proteins from the analysis that could be potential to influence the intensity normalization. For example, ALBU_HUMAN could be over expressed in human tissues, and that is why we may want to remove it when analyzing tissue data.
- Normalize at feature level between ms runs (technical repetitions):
  - When `MS runs > 1` in the sample, the `mean` of all average(`mean`, `median` or `iqr`) in each MS run is calculated(SampleMean)
  - The ratio between SampleMean and the average MS run is used as a reference to scale the original intensity
- Merge peptidoforms across fractions and technical repetitions: Combine technical replicates and fragments from the same sample.
- Normalize the data at the sample level:
  - `globalMedian`: A global median that adjusts the median of all samples.
  - `conditionMedian`: All samples under the same conditions were adjusted to the median value under the current conditions.
- Remove peptides with low frequency if `sample number > 1`: This parameter is applied always unless the user specifies the `--remove_low_frequency_peptides` parameter. The default value is 20% of the samples. If users want to change this threshold, they can use the `--remove_low_frequency_peptides` parameter.
- Assembly peptidoforms to peptides:
A peptidoform is a combination of a `PeptideSequence(Modifications) + Charge + BioReplicate + Fraction` (among other features), and a peptide is a combination of a `PeptideSequence(Canonical) + BioReplicate`. ibaqpy will do:
  - Select peptidoforms with the highest intensity across different modifications, fractions, and technical replicates
  - Merge peptidoforms across different charges and combined into peptides. In order to merge peptidoforms, the package will applied the `sum` of the intensity values of the peptidoforms.
- Intensity transformation to log: The user can specify the `--log2` parameter to transform the peptide intensity values to log2 before normalization.

> Note: At the moment, ibaqpy computes the ibaq values only based on unique peptides. Shared peptides are discarded. However, if a group of proteins share the same unique peptides (e.g., Pep1 -> Prot1;Prot2 and Pep2 -> Prot1;Prot2), the intensity of the proteins is summed and divided by the number of proteins in the group.

#### Calculate the IBAQ Value
First, peptide intensity dataframe was grouped according to protein name, sample name and condition. The protein intensity of each group was summed. Due to the experimental type, the same protein may exhibit missing peptides in different samples, resulting in variations in the number of peptides detected for the protein across different samples. To handle this difference, normalization within the same group can be achieved by using the formula `sum(peptides) / n`(n represents the number of detected peptide segments). Finally, the normalized intensity of the protein is divided by the number of theoretical peptides.See details in `peptides2proteins`.
 
> Note: In all scripts and result files, *uniprot accession* is used as the protein identifier.

### How to install ibaqpy

Ibaqpy is available in PyPI and can be installed using pip:

```asciidoc
pip install ibaqpy
```

You can install the package from code:

1. Clone the repository:

```asciidoc
>$ git clone https://github.com/bigbio/ibaqpy
>$ cd ibaqpy
```

2. Install conda environment:

```asciidoc
>$ mamba env create -f conda-environment.yaml
```

3. Install ibaqpy:

```asciidoc
>$ python setup.py install
```

### Collecting intensity files from quantms.org 

Absolute quantification files have been stored in the following url:

```
https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/quantms-data/
```

Inside each project reanalysis folder, the folder proteomicslfq contains the msstats input file with the structure `{Name of the project}.{Random uuid}.feature.parquet	`. 

E.g. http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/absolute-expression/quantms-data/MSV000079033.1/MSV000079033.1-bd44c7e3-654c-444d-9e21-0f701d6dac94.feature.parquet

### Major commands  


#### Features to peptides

```asciidoc
ibaqpy features2peptides -p tests/PXD003947/PXD003947-feature.parquet -s tests/PXD003947/PXD003947.sdrf.tsv --remove_ids data/contaminants_ids.tsv --remove_decoy_contaminants --remove_low_frequency_peptides --output tests/PXD003947/PXD003947-peptides-norm.csv
```
```asciidoc
Usage: features2peptides.py [OPTIONS]

Options:
  -p, --parquet TEXT              Parquet file import generated by quantms.io
  -s, --sdrf TEXT                 SDRF file import generated by quantms
  --min_aa INTEGER                Minimum number of amino acids to filter
                                  peptides
  --min_unique INTEGER            Minimum number of unique peptides to filter
                                  proteins
  --remove_ids TEXT               Remove specific protein ids from the
                                  analysis using a file with one id per line
  --remove_decoy_contaminants     Remove decoy and contaminants proteins from
                                  the analysis
  --remove_low_frequency_peptides
                                  Remove peptides that are present in less
                                  than 20% of the samples
  --output TEXT                   Peptide intensity file including other all
                                  properties for normalization
  --skip_normalization            Skip normalization step
  --nmethod TEXT                  Normalization method used to normalize
                                  feature intensities for tec
                                  (options: mean, median, iqr, none)
  --pnmethod TEXT                 Normalization method used to normalize
                                  peptides intensities for all samples
                                  (options: globalMedian,conditionMedian,none)
  --log2                          Transform to log2 the peptide intensity
                                  values before normalization
  --save_parquet                  Save normalized peptides to parquet
  --help                          Show this message and exit.
```


####  Compute IBAQ/TPA

```asciidoc
ibaqpy peptides2protein -f Homo-sapiens-uniprot-reviewed-contaminants-decoy-202210.fasta -p PXD017834-peptides.csv -e Trypsin -n -t -r --ploidy 2 --cpc 200 --organism human --output PXD003947.tsv --verbose
``` 

```asciidoc
Usage: peptides2protein [OPTIONS]

Options:
  -f, --fasta TEXT     Protein database to compute IBAQ values
  -p, --peptides TEXT  Peptide identifications with intensities following the
                       peptide intensity output
  -e, --enzyme TEXT    Enzyme used during the analysis of the dataset
                       (default: Trypsin)
  -n, --normalize      Normalize IBAQ values using by using the total IBAQ of
                       the experiment
  --min_aa INTEGER     Minimum number of amino acids to consider a peptide
  --max_aa INTEGER     Maximum number of amino acids to consider a peptide
  -t, --tpa            Whether calculate TPA (is_flag=True)
  -r, --ruler          Whether to use ProteomicRuler (is_flag=True)
  -i, --ploidy         Ploidy number (default=2)
  -m, --organism       Organism source of the data (default human)
  -c, --cpcCellular    protein concentration(g/L)(default 200)
  -o, --output TEXT    Output file with the proteins and ibaq values
  --verbose            Print addition information about the distributions of
                       the intensities, number of peptides remove after
                       normalization, etc.
  --qc_report TEXT     PDF file to store multiple QC images
  --help               Show this message and exit.
```

### How to cite ibaqpy

Wang H, Dai C, Pfeuffer J, Sachsenberg T, Sanchez A, Bai M, Perez-Riverol Y. Tissue-based absolute quantification using large-scale TMT and LFQ experiments. Proteomics. 2023 Oct;23(20):e2300188. doi: [10.1002/pmic.202300188](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/pmic.202300188). Epub 2023 Jul 24. PMID: 37488995.

### Credits 

- [Julianus Pfeuffer](@jpfeuffer)
- [Yasset Perez-Riverol](@ypriverol)
- [Hong Wang](@WangHong007)
- [Ping Zheng](@zprobot)
- [Joshua Klein](@mobiusklein) 
