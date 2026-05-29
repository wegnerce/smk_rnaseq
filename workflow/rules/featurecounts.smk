###############################################################################
# @author:      Carl-Eric Wegner, PhD
# @affiliation: Biomics Group - Chair of Bioinorganic Chemistry
#               Heinrich Heine University Düsseldorf
#
# @contact      carl-eric.wegner@hhu.de
#               https://biomics.hhu.de
#               https://github.com/wegner
#               https://exploringmicrobes.science
###############################################################################

# Note:
# So far the pipeline was only used for AL1 RNAseq data. To make it more 
# flexible, now that we also target other organisms, a little script was 
# added to generate a reference .saf file from a provided .gff3 file as 
# needed.
#
# Check the Snakefile to get an idea about how we check for the presence
# of .gff3 and .saf files.

if REF_GFF3:
    rule gff3_to_saf:
        """ generate .saf annotation file from .gff3 for featureCounts """
        input:
            gff3=REF_GFF3,
        output:
            saf=REF_SAF,
        shell:
            """
            python {workflow.basedir}/scripts/gff3_to_saf.py {input.gff3} {output.saf}
            """


rule featurecounts_readcounts:
    """ deduce readcounts per gene based on the bbmap output and the reference .saf file """
    input:
        bams=expand("results/03_MAPPED/{sample}_mapped_sorted.bam", sample=SAMPLES.index),
        saf=REF_SAF,
    output:
        readcounts="results/04_COUNTS/readcounts_featureCounts.txt",
    conda: "../envs/featurecounts.yaml"
    params:
        pe_flags="-p --countReadPairs" if IS_PAIRED else "",
    threads: 4
    shell:
        """
        featureCounts {params.pe_flags} -T {threads} -F SAF -a {input.saf} -o {output.readcounts} {input.bams}
        """
