#!/usr/bin/env python3

import argparse
from collections import OrderedDict


def fix_gtf(gtf_file, prefix):

    with open("{0}.gtf".format(prefix), "w") as out:

        with open(gtf_file) as inp:

            for line_count, line in enumerate(inp):

                # Identify the header lines
                if line.startswith("#"):

                    # Remove the new line character
                    line = line.strip()

                    # Write header lines
                    out.write("{0}\n".format(line))

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

                annot_data[-1] = '; '.join("{0} \"{1}\"".format(k, v) for (k, v) in annot_info.items())

                # Write annotation
                out.write("{0}\n".format('\t'.join(annot_data)))


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
    print("Parsing GTF annotation file and fixing it")
    fix_gtf(args.annotation_gtf_file, args.output_prefix)


if __name__ == "__main__":
    main()
