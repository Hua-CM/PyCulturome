
# PyCulturome
![License](https://img.shields.io/badge/GPL-3.0-brightgreen) 
![Python](https://img.shields.io/badge/Python-3.8.0-brightgreen)

A one-click Python workflow for high-throughput cultivation and identification of bacteria from the microbiota (Culturome)

## Installation

~~~python
git clone 
cd PyCulturome
pip install -r reuirements.txt
python scripts/pipe.py
~~~

## Major usage

An example

~~~shell
python /path/to/PyCulturome/scripts/pipe.py \
 -1 /path/to/your/fastq1 \
 -2 /path/to/your/fastq2 \
 -f example/input/fwd_bar.tsv \
 -r example/input/rev_bar.tsv \
 -d /home/database/micro/SILVA138_RESCRIPt.fasta \
 -o example/output
~~~

 - **-1**: Read1 input file path
 - **-2**: Read2 input file path
 - **-f**: Forward barcode information in tab-seperated format with two columns: well id \t sequence
 - **-r**: Reverse barcode information in tab-seperated format with two columns: plate id \t sequence
 - **-d**: Taxonomic database for vsearch sequence annotation.
 - **-o**: Output directory

 > At present, only 799F and 1193R with forward 10bp barcode and reverse 6bp barcode are supported.

## Output

### sequence

- **\<name\>.rep.fa**:The fasta file contains representative sequence for ASVs.

### Tables
- **\<name\>.abun.tsv**: The OTU abundance table in wider table format. Generated by usearch directly

- **readscount_well.tsv**: The OTU abundance table in long table format.

- **purified_well.tsv**: The purified well 

- **recommend_well4asv**.tsv
Recommended well numbers for separating each ASV.

- **taxonomy_8.tsv**: Taxonomic information for each ASV

### Figures
- **Positive_Negative_Control.pdf**: The boxplot displays the number of reads in negative and positive control wells.

- **Reads_frequency.pdf**: The histogram plot displays the number of reads per well.

- **Purity_frequency.pdf**: The histogram plot displays the purity of wells (ranges from 0 ~ 100%)

- **Rarefication.pdf**: The rarefication curve of ASVs along with the number of wells

- **Tree.pdf**: The taxonomic tree of isolated strains.

## Other usage 

This script although provide some gadget to help researchers implement other steps
in NGS culturome

### Purification

A considerable amount of Sanger sequencing is required during the purification process. PyCulturome provides a tool for automatically assembling these Sanger sequencing reads and summarizing the results.

~~~shell
python src/sanger.py 
-i /path/to/input_directory
-d /path/to/16S_db
~~~

The NCBI 16S rRNA gene database (16S_ribosomal_RNA.tar.gz) could be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/.



## Dependency
1. [usearch](https://github.com/rcedgar/usearch_old_binaries) > v11
2. [vsearch](https://github.com/torognes/vsearch) >= v2.22
3. python >= 3.8
    - Biopython == 1.80
    - pandas > 2.0.0
    - numpy > 1.23.0
    - ete3 > 3.0.0
    - seaborn > 0.13.0
    - matplotlib > 3.6.0

## TODO

1. Compare and/or combine multiple libraris
2. Add parameter support for additional primers and barcodes.
3. Make the output table iterative development

## Reference
[Culturome-YongxinLiu](https://github.com/YongxinLiu/Culturome/tree/master)