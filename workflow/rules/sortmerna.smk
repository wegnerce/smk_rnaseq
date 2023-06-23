###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################


# SortMeRNA 2.1 appends '.fastq' to specified output files, which leads to 
# messy file names '.fastq.fastq'.
# We prevent that by using parameters that specify output file names minus
# .fastq.
# 
# Needed reference DBs are should be stored under 'resources'.
#
# Here we only use the bacterial SSU and LSU DBs, dealing with an isolate this
# is sufficient.

rule sortmerna_filter_mRNA:
	""" filter mRNA- and rRNA-derived seqs from QCed RNAseq data """
	input:
		read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq",
		read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq",
	output:
		read1="results/02_FILTERED/{sample}_mRNA_" + PAIRS[0] + ".fastq",
		read2="results/02_FILTERED/{sample}_mRNA_" + PAIRS[1] + ".fastq",
	log: 
		"logs/sortmerna/{sample}_stats_filtering.txt",
	conda: "envs/sortmerna.yaml"
	params:
		merged="results/02_FILTERED/{sample}_merged.fastq",
		merged_mRNA="results/02_FILTERED/{sample}_merged_non_rRNA",
		merged_rRNA="results/02_FILTERED/{sample}_merged_rRNA",
		merged_mRNA_fq="results/02_FILTERED/{sample}_merged_non_rRNA.fastq",
		merged_rRNA_fq="results/02_FILTERED/{sample}_merged_rRNA.fastq",
		ssu_bac_db=SILVA_16S_BAC_DB,
		ssu_bac_idx=SILVA_16S_BAC_IDX,
		lsu_bac_db=SILVA_23S_BAC_DB,
		lsu_bac_idx=SILVA_23S_BAC_IDX,
	threads: 4
	shell:
		"""
		merge-paired-reads.sh {input.read1} {input.read2} {params.merged}
		sortmerna --ref {params.ssu_bac_db},{params.ssu_bac_idx}:{params.lsu_bac_db},{params.lsu_bac_idx} --reads {params.merged} --other {params.merged_mRNA} --aligned {params.merged_rRNA} --fastx --best 5 --paired_in -a 4 -v --blast 1 --log 2> {log}
		unmerge-paired-reads.sh {params.merged_mRNA_fq} {output.read1} {output.read2}
		rm {params.merged}
		rm {params.merged_rRNA_fq}
		rm {params.merged_mRNA_fq}
		"""

