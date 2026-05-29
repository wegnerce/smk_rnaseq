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
    # Extract geneID from .gff3 attributes column.
    # 
    # Looks for locus_tag, ID, and Name. Used the last 
    # key=value pair if none are found (PLAN B).
    # In the original version, the geneID was parsed from the
    # "locus_tag" attribute, however, this is not always
    # present in .gff3 files.
    attr_dict = {}
    for attr in attributes.strip().split(";"):
        if "=" in attr:
            key, _, value = attr.partition("=")
            attr_dict[key.strip()] = value.strip()

    for key in ("locus_tag", "ID", "Name"):
        if key in attr_dict:
            return attr_dict[key]

    # PLAN B: 
    # Return the value of the last key=value pair, 
    # if any attributes are present. This is not ideal, 
    # but it allows us to recover at least some geneIDs from 
    # .gff3 files that don't follow the expected scheme
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
    parser.add_argument("gff3", help="Input .gff3 file")
    parser.add_argument("saf", help="Output .saf file")
    args = parser.parse_args()

    try:
        convert(args.gff3, args.saf)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
