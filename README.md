# GTF to split annotations BED files

This tool splits a comprehensive GTF file into 8 separate BED files:

- *PREFIX*.exons.bed - All exons positions
- *PREFIX*.introns.bed - All introns positions
- *PREFIX*.utrs.bed - All UTRs positions
- *PREFIX*.intergenic.bed - All intergenic positions
- *PREFIX*.genes.bed - All genes
- *PREFIX*.miRNA.bed - All miRNA genes
- *PREFIX*.rRNA.bed - All rRNA genes
- *PREFIX*.lincRNA.bed - All lincRNA genes

## Requirements

There are three requirement for this script:

- Bash
- Python 3
- Bedtools 

## Usage

    ./main.sh REFERENCE.fasta.fai ANNOTATION.gtf PREFIX

### Example

    ./main.sh GRCh38.genome.fa.fai gencode.v28.annotation.gtf gencode.v28