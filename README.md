# ibaqpy

IBAQPy compute IBAQ values for Triqler outputs of Sample specific datasets. Additional tasks:

- Remove contaminants, decoy and proteins with non-tryptic peptides
- Normalize IBAQ values for each protein.

### Generate the peptide intensity files

```asciidoc
python  peptide_file_generation.py --mztab data/PXD004682-out.mzTab.gz --msstats data/PXD004682-out_msstats.csv.gz --triqler data/PXD004682-out_triqler.tsv.gz --sdrf data/PXD004682.sdrf.tsv.gz --output data/PXD004682-Peptide-Intensities.tsv --compression
```


Credit: Julianus and Yasset Perez-Riverol
