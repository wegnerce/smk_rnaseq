###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################

rule featurecounts_readcounts:
	"""deduced readcounts per gene based on the bbmap output and the reference .saf file""" 
	input:
		bams=expand("results/03_MAPPED/{sample}_mapped_sorted.bam", sample=SAMPLES.index),
	output:
		readcounts="results/04_COUNTS/readcounts_featureCounts.txt",
	conda: "../envs/featurecounts.yaml"
	params:
		saf=REF_SAF
	threads: 4
	shell:
		"""
		featureCounts -p -a {params.saf} -o {output.readcounts} -F SAF -t {input.bams}
		"""

