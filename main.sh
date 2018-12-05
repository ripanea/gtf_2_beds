#!/usr/bin/env bash

# Initialize main variables
FASTA_IDX=$1
GTF=$2
PREFIX=$3

# Create a temporary directory
mkdir tmp

# Generate the bed files
./make_beds.py ${GTF} tmp/

# Sort and merge each bed file
for bed_file in $(find ./tmp -name "*.bed" -printf "%f\n");
do
    # Merge all files except "genes"
    if [[ ! ${bed_file} =~ ".genes." ]]
    then
        bedtools sort -faidx ${FASTA_IDX} -i tmp/${bed_file} > tmp/sorted${bed_file}
        bedtools merge -i tmp/sorted${bed_file} > tmp/merged${bed_file}
        bedtools sort -faidx ${FASTA_IDX} -i tmp/merged${bed_file} > ${PREFIX}${bed_file}
    else
        bedtools sort -faidx ${FASTA_IDX} -i tmp/${bed_file} > ${PREFIX}${bed_file}
    fi
done

# Generate genome.txt file for bedtools complement
cut -f1,2 ${FASTA_IDX} > tmp/genome.txt

# Generate the intergenic file
bedtools complement -i ${PREFIX}.genes.bed -g tmp/genome.txt > ${PREFIX}.intergenic.bed

# Remove temporary directory
rm -r tmp