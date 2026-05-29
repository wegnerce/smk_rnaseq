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
# The "slow" setting is used for increased sensitivity, however, given that
# we deal with data from one organism, one could benchmark if this is really
# needed, e.g. "slow k=11" vs just "k13" - TO DO.

if IS_PAIRED:
    rule bbmap_mapping:
        """ map the output from sortmerna onto the reference genome (paired-end) """
        input:
            read1="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
            read2="results/02_FILTERED/{sample}_mRNA_" + PAIRS[1] + ".fastq",
        output:
            mapped=temp("results/03_MAPPED/{sample}_mapped.bam"),
            mapped_sorted="results/03_MAPPED/{sample}_mapped_sorted.bam",
            mapping_stats="logs/bbmap/{sample}_stats_mapping.txt",
            flagstat="logs/bbmap/{sample}_flagstat.txt",
            samtools_stats="logs/bbmap/{sample}_samtools_stats.txt",
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
            samtools stats {output.mapped_sorted} > {output.samtools_stats}
            """
else:
    rule bbmap_mapping:
        """ map the output from sortmerna onto the reference genome (single-end) """
        input:
            read="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
        output:
            mapped=temp("results/03_MAPPED/{sample}_mapped.bam"),
            mapped_sorted="results/03_MAPPED/{sample}_mapped_sorted.bam",
            mapping_stats="logs/bbmap/{sample}_stats_mapping.txt",
            flagstat="logs/bbmap/{sample}_flagstat.txt",
            samtools_stats="logs/bbmap/{sample}_samtools_stats.txt",
        conda: "../envs/bbmap.yaml"
        params:
            ref_genome=REF_GENOME,
        threads: 4
        shell:
            """
            bbmap.sh slow k=11 threads={threads} \
                in={input.read} \
                ref={params.ref_genome} \
                out={output.mapped} statsfile={output.mapping_stats}
            samtools sort -@ {threads} {output.mapped} -o {output.mapped_sorted}
            samtools index {output.mapped_sorted}
            samtools flagstat {output.mapped_sorted} > {output.flagstat}
            samtools stats {output.mapped_sorted} > {output.samtools_stats}
            """
