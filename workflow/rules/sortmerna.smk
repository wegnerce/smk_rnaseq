###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################

# SortMeRNA 4.x handles paired-end reads natively — the merge-paired-reads.sh /
# unmerge-paired-reads.sh approach from 2.1 is no longer needed, and pre-built
# indexes are no longer required (built on-the-fly in --workdir).
#
# Each sample gets its own --workdir to prevent kvdb conflicts when samples
# run in parallel. The workdir is removed after completion.
#
# --out2 splits non-rRNA output into _fwd.fastq / _rev.fastq; these are
# renamed to match the R1/R2 naming convention used downstream.

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
