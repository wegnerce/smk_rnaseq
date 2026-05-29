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
# The current set-up with fastqc and bbduk is inefficient, fastQC has limited
# multi-threading capabilities.
# Using fastp, would reduce the number of rules, as it combines the 
# functionality of fastQC and bbduk - TO DO.

rule fastqc_raw:
    """ generate QC reports for the raw data """
    input:
        read=raw_data_dir + "/{sample}_{pair}.fastq.gz",
    output:
        qual="logs/fastqc/raw/{sample}_{pair}_fastqc.html",
        zip ="logs/fastqc/raw/{sample}_{pair}_fastqc.zip",
    resources:
        mem_mb=2000
    conda: "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.read} -f fastq --outdir logs/fastqc/raw
        """

if IS_PAIRED:
    rule bbduk_trim:
        """ adapter removal and sequence trimming with bbduk (paired-end) """
        input:
            read1=raw_data_dir + "/{sample}_" + PAIRS[0] + ".fastq.gz",
            read2=raw_data_dir + "/{sample}_" + PAIRS[1] + ".fastq.gz",
        output:
            read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq",
            read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq",
            trim_stats="logs/bbduk/{sample}_stats_QC.txt",
        conda: "../envs/bbmap.yaml"
        params:
            adapter=ADAPTER,
            minlen=config["BBDUK"]["minlen"],
            trimq=config["BBDUK"]["trimq"],
            ftl=config["BBDUK"]["force_trim_left"],
        threads: 4
        shell:
            """
            bbduk.sh -Xmx8g in1={input.read1} in2={input.read2} \
                out1={output.read1} out2={output.read2} \
                stats={output.trim_stats} \
                minlen={params.minlen} qtrim=rl trimq={params.trimq} ktrim=r k=25 mink=11 \
                ref={params.adapter} hdist=1 ftl={params.ftl} threads={threads}
            """
else:
    rule bbduk_trim:
        """ adapter removal and sequence trimming with bbduk (single-end) """
        input:
            read=raw_data_dir + "/{sample}_" + PAIRS[0] + ".fastq.gz",
        output:
            read="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq",
            trim_stats="logs/bbduk/{sample}_stats_QC.txt",
        conda: "../envs/bbmap.yaml"
        params:
            adapter=ADAPTER,
            minlen=config["BBDUK"]["minlen"],
            trimq=config["BBDUK"]["trimq"],
            ftl=config["BBDUK"]["force_trim_left"],
        threads: 4
        shell:
            """
            bbduk.sh -Xmx8g in={input.read} \
                out={output.read} \
                stats={output.trim_stats} \
                minlen={params.minlen} qtrim=rl trimq={params.trimq} ktrim=r k=25 mink=11 \
                ref={params.adapter} hdist=1 ftl={params.ftl} threads={threads}
            """

rule fastqc_trim:
    """ generate QC reports for the trimmed data """
    input:
        read="results/01_TRIMMED/{sample}_trimmed_{pair}.fastq",
    output:
        qual="logs/fastqc/trimmed/{sample}_trimmed_{pair}_fastqc.html",
        zip ="logs/fastqc/trimmed/{sample}_trimmed_{pair}_fastqc.zip",
    resources:
        mem_mb=2000
    conda: "../envs/fastqc.yaml"
    shell:
        """
        fastqc {input.read} -f fastq --outdir logs/fastqc/trimmed
        """

rule read_stats:
    """ count reads and bases at each processing stage for a quick dataset overview """
    input:
        raw       = expand(raw_data_dir + "/{sample}_{pair}.fastq.gz",
                           sample=SAMPLES.index, pair=PAIRS),
        trimmed   = expand("results/01_TRIMMED/{sample}_trimmed_{pair}.fastq",
                           sample=SAMPLES.index, pair=PAIRS),
        filtered  = expand("results/02_FILTERED/{sample}_mRNA_{pair}.fastq",
                           sample=SAMPLES.index, pair=PAIRS),
        bam_stats = expand("logs/bbmap/{sample}_samtools_stats.txt",
                           sample=SAMPLES.index),
    output:
        tsv="results/00_QC/read_stats.tsv",
    run:
        import gzip

        def count_fastq(path):
            opener = gzip.open if path.endswith('.gz') else open
            reads = bases = 0
            with opener(path, 'rt') as fh:
                for i, line in enumerate(fh):
                    if i % 4 == 1:
                        reads += 1
                        bases += len(line.rstrip('\n'))
            return reads, bases

        def parse_samtools_stats(path):
            reads = bases = 0
            with open(path) as fh:
                for line in fh:
                    if not line.startswith('SN'):
                        continue
                    parts = line.split('\t')
                    if parts[1].strip() == 'reads mapped:':
                        reads = int(parts[2])
                    elif parts[1].strip() == 'bases mapped (cigar):':
                        bases = int(parts[2])
            return reads, bases

        samples = list(SAMPLES.index)
        n_pairs = len(PAIRS)
        bam_stats_list = list(input.bam_stats)

        with open(output.tsv, 'w') as out:
            out.write("sample\tstage\treads\ttotal_bp\tpct_reads\tpct_bp\n")
            for idx, sample in enumerate(samples):
                counts = {}
                for stage, file_list in [
                    ("raw",      list(input.raw)),
                    ("trimmed",  list(input.trimmed)),
                    ("filtered", list(input.filtered)),
                ]:
                    # expand() is sample-major: [s0_R1, s0_R2, s1_R1, s1_R2, ...]
                    files    = [file_list[idx * n_pairs + p] for p in range(n_pairs)]
                    reads, bases = count_fastq(files[0])
                    total_bp = bases + sum(count_fastq(f)[1] for f in files[1:])
                    counts[stage] = (reads, total_bp)

                counts["mapped"] = parse_samtools_stats(bam_stats_list[idx])

                raw_reads, raw_bp = counts["raw"]
                for stage in ("raw", "trimmed", "filtered", "mapped"):
                    reads, total_bp = counts[stage]
                    pct_reads = 100.0 * reads    / raw_reads if raw_reads else 0.0
                    pct_bp    = 100.0 * total_bp / raw_bp    if raw_bp    else 0.0
                    out.write(f"{sample}\t{stage}\t{reads}\t{total_bp}\t{pct_reads:.1f}\t{pct_bp:.1f}\n")


rule multiqc:
    """ aggregate all QC reports into a single MultiQC report """
    input:
        expand("logs/fastqc/raw/{sample}_{pair}_fastqc.zip", sample=SAMPLES.index, pair=PAIRS),
        expand("logs/fastqc/trimmed/{sample}_trimmed_{pair}_fastqc.zip", sample=SAMPLES.index, pair=PAIRS),
        expand("logs/bbduk/{sample}_stats_QC.txt", sample=SAMPLES.index),
        expand("logs/bbmap/{sample}_stats_mapping.txt", sample=SAMPLES.index),
        expand("logs/bbmap/{sample}_flagstat.txt", sample=SAMPLES.index),
        expand("logs/bbmap/{sample}_samtools_stats.txt", sample=SAMPLES.index),
        expand("logs/sortmerna/{sample}_stats_filtering.txt", sample=SAMPLES.index),
    output:
        report="results/00_QC/multiqc_report.html",
    conda: "../envs/multiqc.yaml"
    shell:
        """
        multiqc logs/ --outdir results/00_QC --force
        """
