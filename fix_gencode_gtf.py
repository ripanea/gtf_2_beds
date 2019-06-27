#!/usr/bin/env python3

import argparse
from collections import OrderedDict


def fix_gtf(gtf_file):

    with open(gtf_file) as inp:

        # Holds the header lines
        header_lines = list()

        # Holds the annotation lines
        annotation_lines = list()

        for line_count, line in enumerate(inp):

            # Skip header
            if line.startswith("#"):
                header_lines.append(line.strip())
                continue

            # Parse annotation
            annot_data = line.strip().split("\t")
            annot_info = OrderedDict()

            # Iterate through all the annotations
            for info in annot_data[-1].strip(";").split(";"):

                # Split info data
                info_type, info_value = info.strip().split(" ")

                # Record info data
                annot_info[info_type] = info_value.strip('"')

            annot_data[-1] = '; '.join("{0} \'{1}\'".format(k, v) for (k, v) in annot_info.items())

            annotation_lines.append('\t'.join(annot_data))

        return header_lines, annotation_lines


def write_gtf(header, annotations, prefix):

    with open("{0}.gtf".format(prefix), "w") as out:

        # Write the header lines
        for header_line in header:
            out.write("{0}\n".format(header_line))

        # Write the annotations
        for annotation_line in annotations:
            out.write("{0}\n".format(annotation_line))


def configure_argparser(argparser_obj):

    # Annotation GTF file
    argparser_obj.add_argument("annotation_gtf_file",
                               action="store",
                               help="Uncompressed annotation file (.gtf)")

    # Output prefix
    argparser_obj.add_argument("output_prefix",
                               action="store",
                               help="Path prefix of the output GTF files")


def main():

    # Generate argument parser
    argparser_obj = argparse.ArgumentParser()

    # Configure argparser
    configure_argparser(argparser_obj)

    # Parse arguments
    args = argparser_obj.parse_args()

    # Parse annotation file
    print("Step 1: Parsing GTF annotation file")
    header, annot_info = fix_gtf(args.annotation_gtf_file)

    # Writing output BED files
    print("Step 2: Writing GTF file with correct specifications")
    write_gtf(header, annot_info, args.output_prefix)


if __name__ == "__main__":
    main()
