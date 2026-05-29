# -*- coding: utf-8 -*-
"""
@author:        Carl-Eric Wegner, PhD
@affiliation:   Biomics Group - Chair of Bioinorganic Chemistry
                Heinrich Heine University Düsseldorf

@contact:       carl-eric.wegner@hhu.de
                https://biomics.hhu.de
                https://github.com/wegnerce
                https://exploringmicrobes.science
"""

import sys
import argparse

def parse_gene_id(attributes):
    """Extract gene ID from GFF3 attributes column.

    Tries locus_tag, ID, and Name in order of preference.
    Falls back to the last key=value pair if none are found.
    In the original version, the gene ID was parsed from the 
    "locus_tag" attribute, but this is not guaranteed to be 
    present in all GFF3 files, so we try multiple keys.
    """
    attr_dict = {}
    for attr in attributes.strip().split(";"):
        if "=" in attr:
            key, _, value = attr.partition("=")
            attr_dict[key.strip()] = value.strip()

    for key in ("locus_tag", "ID", "Name"):
        if key in attr_dict:
            return attr_dict[key]

    # PLAN B: return the value of the last key=value pair, 
    # if any attributes are present. This is a bit of a hack, 
    # but it allows us to recover at least some gene IDs from 
    # GFF3 files that don't follow the expected scheme
    if attr_dict:
        return list(attr_dict.values())[-1]

    return None


def convert(gff3, saf):
    with open(gff3, "r") as infile, open(saf, "w") as outfile:
        outfile.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for lineno, raw in enumerate(infile, 1):
            if raw.startswith("#"):
                continue
            line = raw.rstrip("\n").split("\t")
            if len(line) < 9:
                continue
            if line[2] != "gene":
                continue

            chrom = line[0]
            start = line[3]
            end = line[4]
            strand = line[6]
            gene_id = parse_gene_id(line[8])

            if gene_id is None:
                print(f"Warning: no GeneID on line {lineno}, OMIT.",
                      file=sys.stderr)
                continue

            outfile.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert a .gff3 annotation file to .saf format for featureCounts."
                    "»» https://subread.sourceforge.net/featureCounts.html ««"
    )
    parser.add_argument("gff3", help="Input GFF3 file")
    parser.add_argument("saf", help="Output SAF file")
    args = parser.parse_args()

    try:
        convert(args.gff3, args.saf)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
