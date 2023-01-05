# ibaqpy

ibaqpy compute IBAQ values for proteins starting from an msstats input file and a SDRF experimental design file. This package provides multiple tools: 

- `peptide_file_generation.py`: generate a peptide file from a msstats input file and a SDRF experimental design file. 

- `peptide_normalization.py`: Normalize the input output file from script `peptide_file_generation.py`. It includes multiple steps such as peptidoform normalization, peptidorom to peptide summarization, peptide intesity normalization, and imputation. 

- `compute_ibaq.py`: Compute IBAQ values from the output file from script `peptide_normalization.py`.

### Collecting intensity files 

Absolute quantification files has been store in the following url: 

```
http://ftp.pride.ebi.ac.uk/pride/data/proteomes/absolute-expression/
```

This FTP structure the information with the following form: 

- Project (Project reanalyzed)
  - Project Sample results (This can be one folder or multiple depending if the samples are analyzed together or splitted)

Collecting all the Sample intensities for all the projects can be moved to the folder all-data using the script: 

Moving triqler: 

```bash
find ./ -name "out_triqler.tsv" | while read line; do echo $line ; project="$(basename "$(dirname "$(dirname "$line")")")"; echo $project; cp -v $line all-data/${project}-triqler.tsv; done
```

Moving mzTab: 

```bash 
find ./ -name "out.mzTab" | while read line; do echo $line ; project="$(basename "$(dirname "$(dirname "$line")")")"; echo $project; cp -v $line all-data/${project}.mzTab; done
```

Moving msstats: 

```bash
find ./ -name "out_msstats.csv" | while read line; do echo $line ; project="$(basename "$(dirname "$(dirname "$line")")")"; echo $project; cp -v $line all-data/${project}-msstats.tsv; done
```

### Generate the peptide intensity files

```asciidoc
python  peptide_file_generation.py --mztab data/PXD004682-out.mzTab.gz --msstats data/PXD004682-out_msstats.csv.gz --triqler data/PXD004682-out_triqler.tsv.gz --sdrf data/PXD004682.sdrf.tsv.gz --output data/PXD004682-Peptide-Intensities.tsv --compression
```

### Credits 

Julianus and Yasset Perez-Riverol
