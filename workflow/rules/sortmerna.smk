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
# Long overdue update from SortMeRNA 2.1b to 4.x. 2.1b was not able to handle
# paired-end data natively, which made this step unnecessary slow.
#
# Using dedicated (temporary) working directories makes the processing
# cleaner, especially if things go wrong.
# 
# Keep in mind, SortMeRNA output .fq vs .fastq --> .f*q 
# --> took me ages to see this during debugging

if IS_PAIRED:
    rule sortmerna_filter_mRNA:
        """ filter rRNA-derived sequences from QCed paired-end RNAseq data """
        input:
            read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq",
            read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq",
        output:
            read1="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
            read2="results/02_FILTERED/{sample}_mRNA_" + PAIRS[1] + ".fastq",
        log:
            "logs/sortmerna/{sample}_stats_filtering.txt",
        conda: "../envs/sortmerna.yaml"
        resources:
            mem_mb=32000
        params:
            ssu_bac_db=SILVA_16S_BAC_DB,
            lsu_bac_db=SILVA_23S_BAC_DB,
            workdir="results/02_FILTERED/{sample}_sortmerna_wd",
            other_prefix="results/02_FILTERED/{sample}_mRNA_tmp",
            rRNA_prefix="results/02_FILTERED/{sample}_rRNA_tmp",
        threads: 4
        shell:
            """
            sortmerna \
                --ref {params.ssu_bac_db} \
                --ref {params.lsu_bac_db} \
                --reads {input.read1} \
                --reads {input.read2} \
                --aligned {params.rRNA_prefix} \
                --other {params.other_prefix} \
                --workdir {params.workdir} \
                --fastx --out2 --paired_out \
                --threads {threads} \
                -v 2> {log}
            mv {params.other_prefix}_fwd.f*q {output.read1}
            mv {params.other_prefix}_rev.f*q {output.read2}
            rm -rf {params.workdir}
            rm -f {params.rRNA_prefix}_fwd.f*q {params.rRNA_prefix}_rev.f*q
            """
else:
    rule sortmerna_filter_mRNA:
        """ filter rRNA-derived sequences from QCed single-end RNAseq data """
        input:
            read="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq",
        output:
            read="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
        log:
            "logs/sortmerna/{sample}_stats_filtering.txt",
        conda: "../envs/sortmerna.yaml"
        resources:
            mem_mb=32000
        params:
            ssu_bac_db=SILVA_16S_BAC_DB,
            lsu_bac_db=SILVA_23S_BAC_DB,
            workdir="results/02_FILTERED/{sample}_sortmerna_wd",
            other_prefix="results/02_FILTERED/{sample}_mRNA_tmp",
            rRNA_prefix="results/02_FILTERED/{sample}_rRNA_tmp",
        threads: 4
        shell:
            """
            sortmerna \
                --ref {params.ssu_bac_db} \
                --ref {params.lsu_bac_db} \
                --reads {input.read} \
                --aligned {params.rRNA_prefix} \
                --other {params.other_prefix} \
                --workdir {params.workdir} \
                --fastx \
                --threads {threads} \
                -v 2> {log}
            mv {params.other_prefix}.f*q {output.read}
            rm -rf {params.workdir}
            rm -f {params.rRNA_prefix}.f*q
            """
