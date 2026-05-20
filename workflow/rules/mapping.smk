###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################

rule bbmap_mapping:
    """ map the output from sortmerna onto the reference genome """
    input:
        read1="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
        read2="results/02_FILTERED/{sample}_mRNA_" + PAIRS[1] + ".fastq",
    output:
        mapped=temp("results/03_MAPPED/{sample}_mapped.bam"),
        mapped_sorted="results/03_MAPPED/{sample}_mapped_sorted.bam",
        mapping_stats="logs/bbmap/{sample}_stats_mapping.txt",
        flagstat="logs/bbmap/{sample}_flagstat.txt",
    conda: "../envs/bbmap.yaml"
    params:
        ref_genome=REF_GENOME,
    threads: 4
    shell:
        """
        bbmap.sh slow k=11 threads={threads} \
            in={input.read1} in2={input.read2} \
            ref={params.ref_genome} \
            out={output.mapped} statsfile={output.mapping_stats}
        samtools sort -@ {threads} {output.mapped} -o {output.mapped_sorted}
        samtools index {output.mapped_sorted}
        samtools flagstat {output.mapped_sorted} > {output.flagstat}
        """
