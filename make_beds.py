#!/usr/bin/env python3
import argparse
from collections import OrderedDict


class BEDRecord(object):

    def __init__(self, chrom, start, end, info=None):

        self.chrom = chrom
        self.start = start-1
        self.end = end

        self.info = info

    def __str__(self):
        if self.info is not None:
            return "{0}\t{1}\t{2}\t{3}\n".format(self.chrom, self.start, self.end, "\t".join(self.info))
        else:
            return "{0}\t{1}\t{2}\n".format(self.chrom, self.start, self.end)

    def __repr__(self):
        return self.__str__()


class Gene(object):

    def __init__(self, gene_id, chrom, start, end, gene_type):

        self.gene_id = gene_id

        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)

        self.gene_type = gene_type

        self.transcripts = OrderedDict()

    def add_transcript(self, transcript):
        self.transcripts[transcript.transcript_id] = transcript

    def get_exons(self):

        bed_records = []
        for transcript in self.transcripts.values():
            bed_records.extend(transcript.get_exons())

        return bed_records

    def get_introns(self):

        bed_records = []
        for transcript in self.transcripts.values():
            bed_records.extend(transcript.get_introns())

        return bed_records

    def get_utrs(self):

        bed_records = []
        for transcript in self.transcripts.values():
            bed_records.extend(transcript.get_utrs())

        return bed_records

    def as_bed(self):
        return BEDRecord(self.chrom, self.start, self.end)


class Transcript(object):

    def __init__(self, trans_id, chrom, start, end, strand):

        self.transcript_id = trans_id

        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)

        self.is_forward = strand == "+"

        self.exons = []
        self.utrs = []

    def add_exon(self, start, end):

        if self.is_forward:
            self.exons.append((int(start), int(end)))
        else:
            self.exons.insert(0, (int(start), int(end)))

    def get_exons(self):
        return [BEDRecord(self.chrom, start, end) for start, end in self.exons]

    def get_introns(self):

        introns = []

        for exon_id in range(1, len(self.exons)):

            prev_exon_end = self.exons[exon_id-1][1]
            cur_exon_start = self.exons[exon_id][0]

            introns.append(BEDRecord(self.chrom, prev_exon_end+1, cur_exon_start-1))

        return introns

    def add_utr(self, start, end):
        self.utrs.append((int(start), int(end)))

    def get_utrs(self):
        return [BEDRecord(self.chrom, start, end) for start, end in self.utrs]


def parse_gtf(gtf_file):

    # Initialize gene dataset
    genes = OrderedDict()

    with open(gtf_file) as inp:
        for line_count, line in enumerate(inp):

            # Skip header
            if line.startswith("#"):
                continue

            # Parse annotation
            annot_data = line.strip().split("\t")
            annot_info = {}
            for info in annot_data[-1].strip(";").split(";"):

                # Split info data
                info_type, info_value = info.strip().split(" ", 1)

                # Record info data
                annot_info[info_type] = info_value.strip('"')

            # Check if annotation is new gene
            if annot_data[2] == "gene":

                # Create new gene object
                new_gene = Gene(gene_id=annot_info["gene_id"],
                                chrom=annot_data[0],
                                start=annot_data[3],
                                end=annot_data[4],
                                gene_type=annot_info["gene_type"])

                # Record new gene transcript
                genes[annot_info["gene_id"]] = new_gene

            # Check if annotation is new transcript
            elif annot_data[2] == "transcript":

                # Create new transcript object
                new_trans = Transcript(trans_id=annot_info["transcript_id"],
                                       chrom=annot_data[0],
                                       start=annot_data[3],
                                       end=annot_data[4],
                                       strand=annot_data[6])

                # Add new transcript to correct gene
                gene_obj = genes[annot_info["gene_id"]]
                gene_obj.add_transcript(new_trans)

            # Check if annotation is new exon
            elif annot_data[2] == "exon":

                # Add new exon to correct transcript
                gene_obj = genes[annot_info["gene_id"]]
                trans_obj = gene_obj.transcripts[annot_info["transcript_id"]]
                trans_obj.add_exon(annot_data[3], annot_data[4])

            # Check if annotation is UTR
            elif annot_data[2] == "UTR":

                # Add new exon to correct transcript
                gene_obj = genes[annot_info["gene_id"]]
                trans_obj = gene_obj.transcripts[annot_info["transcript_id"]]
                trans_obj.add_utr(annot_data[3], annot_data[4])

            if line_count % 1000 == 0:
                print("  Processed {0} lines...\r".format(line_count), end="")

    return genes


def obtain_chr_lengths(fasta_idx_file):

    # Initialize the dictionary of chromosome lengths
    chr_lengths = {}

    # Obtain chromosome lengths
    with open(fasta_idx_file) as inp:
        for line in inp:

            # Parse chromosome info
            chr_name, length, _ = line.split("\t", 2)

            # Record the chromosome length
            chr_lengths[chr_name] = length

    return chr_lengths


def configure_argparser(argparser_obj):

    # Annotation GTF file
    argparser_obj.add_argument("annotation_gtf_file",
                               action="store",
                               help="Uncompressed annotation file (.gtf)")

    # Output prefix
    argparser_obj.add_argument("output_prefix",
                               action="store",
                               help="Path prefix of the output BED files")


def write_exons_bed(genes, prefix):

    with open("{0}.exons.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            for exon in gene_obj.get_exons():
                out.write(str(exon))


def write_introns_bed(genes, prefix):

    with open("{0}.introns.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            for intron in gene_obj.get_introns():
                out.write(str(intron))


def write_utrs_bed(genes, prefix):

    with open("{0}.utrs.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            for utr in gene_obj.get_utrs():
                out.write(str(utr))


def write_rrna_bed(genes, prefix):

    with open("{0}.rRNA.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            if gene_obj.gene_type.lower() == "rrna":
                out.write(str(gene_obj.as_bed()))


def write_mirna_bed(genes, prefix):

    with open("{0}.miRNA.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            if gene_obj.gene_type.lower() == "mirna":
                out.write(str(gene_obj.as_bed()))


def write_lincrna_bed(genes, prefix):

    with open("{0}.lincRNA.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            if gene_obj.gene_type.lower() == "lincrna":
                out.write(str(gene_obj.as_bed()))


def write_genes_bed(genes, prefix):

    with open("{0}.genes.bed".format(prefix), "w") as out:

        for gene_obj in genes.values():

            out.write(str(gene_obj.as_bed()))


def main():

    # Generate argument parser
    argparser_obj = argparse.ArgumentParser()

    # Configure argparser
    configure_argparser(argparser_obj)

    # Parse arguments
    args = argparser_obj.parse_args()

    # Parse annotation file
    print("Step 1: Parsing GTF annotation file")
    genes = parse_gtf(args.annotation_gtf_file)

    # Writing output BED files
    print("Step 2: Writing BED file with exon positions")
    write_exons_bed(genes, args.output_prefix)

    print("Step 3: Writing BED file with intron positions")
    write_introns_bed(genes, args.output_prefix)

    print("Step 4: Writing BED file with UTR positions")
    write_utrs_bed(genes, args.output_prefix)

    print("Step 5: Writing BED file with rRNA genes")
    write_rrna_bed(genes, args.output_prefix)

    print("Step 6: Writing BED file with miRNA genes")
    write_mirna_bed(genes, args.output_prefix)

    print("Step 7: Writing BED file with lincRNA genes")
    write_lincrna_bed(genes, args.output_prefix)

    print("Step 8: Writing BED file with all the genes")
    write_genes_bed(genes, args.output_prefix)


if __name__ == "__main__":
    main()
