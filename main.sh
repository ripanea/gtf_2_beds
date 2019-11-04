#!/usr/bin/env bash

#Script Name: main.sh

#Description: This script takes a GENCODE version of GTF annotation file with reference genome FASTA index file to
# generate BED files for different gene biotype such as exon, utr, intron, etc.

#Required files and softwares:

##Files:
#1. GENCODE version of comprehensive gene annotation GTF file for all chromosomes and contigs
#2. GENCODE version of reference genome fasta index file for all chromosomes and contigs

##Software Requiements:
#1. BEDTOOLS

#Syntax:
#./main.sh -g <GTF_file> -i <index_fa_file> -o <output_prefix>

#Usage:
#./make_bed_files.sh -g gencode.v28.primary_assembly.annotation.gtf
#                    -i GRCh38.primary_assembly.genome.fa.fai gencode.v28.primary_assembly
#                    -o gencode.v28

# modify the behavior of the shell that will exit whenever a command exits with a non zero status (with some exceptions)
set -e
# makes the shell to treat undefined variables as errors
set -u
# changes the way command inside pipe are evaluated
set -o pipefail

#usage message to display on failure
usage() {
    echo "script usage: $(basename $0) [-g gtf file] [-i reference index file ] [-o output prefix]
                                       [-h|-? print usage]
         "
}

# print the usage message if no arguments to program provided
if [[ $# -eq 0 ]]
    then
    usage
    exit 1;
fi

# check the provided arguments
while getopts 'g:i:o:h' OPTION; do
    case "${OPTION}" in
        g)
            GTF=`realpath ${OPTARG}`
            ;;
        i)
            FASTA_IDX=`realpath ${OPTARG}`
            ;;
        o)
            PREFIX=`realpath ${OPTARG}`
            ;;
        h|?)
            usage
            exit 1
            ;;
    esac
done

# Create a temporary directory
mkdir tmp

# Generate the bed files
./make_beds.py "${GTF}" tmp/

# Sort and merge each bed file
for bed_file in $(find ./tmp -name "*.bed" -printf "%f\n");
do
    # Merge all files except "genes"
    if [[ ! ${bed_file} =~ ".genes." ]]
    then
        bedtools sort -faidx "${FASTA_IDX}" -i tmp/${bed_file} > tmp/sorted${bed_file}
        bedtools merge -i tmp/sorted${bed_file} > tmp/merged${bed_file}
        bedtools sort -faidx "${FASTA_IDX}" -i tmp/merged${bed_file} > "${PREFIX}"${bed_file}
    else
        bedtools sort -faidx "${FASTA_IDX}" -i tmp/${bed_file} > "${PREFIX}"${bed_file}
    fi
done

# Generate genome.txt file for bedtools complement
cut -f1,2 "${FASTA_IDX}" > tmp/genome.txt

# Generate the intergenic file
bedtools complement -i "${PREFIX}".genes.bed -g tmp/genome.txt > "${PREFIX}".intergenic.bed

# Remove temporary directory
rm -r tmp

exit;
