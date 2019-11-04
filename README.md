# GTF to split annotations BED files

This tool splits a comprehensive GTF file into 8 separate BED files:

- *PREFIX*.exons.bed - All exons positions
- *PREFIX*.introns.bed - All introns positions
- *PREFIX*.utrs.bed - All UTRs positions
- *PREFIX*.intergenic.bed - All intergenic positions
- *PREFIX*.genes.bed - All genes
- *PREFIX*.miRNA.bed - All miRNA (pseudo)genes
- *PREFIX*.rRNA.bed - All rRNA (pseudo)genes
- *PREFIX*.lincRNA.bed - All lincRNA (pseudo)genes

## Requirements

There are three requirement for this script:

- Bash
- Python 3
- Bedtools 

## Usage

    ./main.sh -i REFERENCE.fasta.fai -g ANNOTATION.gtf -o PREFIX

### Example

    ./main.sh -i GRCh38.genome.fa.fai -g gencode.v28.annotation.gtf -o gencode.v28