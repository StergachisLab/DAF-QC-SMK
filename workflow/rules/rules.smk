# Basic rules

rule deduplicate:
    input:
        data= get_input_bam
    params:
        dup_end_len=DUP_END_LENGTH,
        min_id=MIN_ID_PERC
    conda:
        "../envs/cmd.yaml"
    output:
        dedup_bam=temp("temp/{sm}/align/{sm}.dupmark.bam"),
        dup_index=temp("temp/{sm}/align/{sm}.dupmark.bam.bai")
    log:
        "logs/{sm}/dedup/{sm}.dedup.log"
    threads: 8
    shell:
        """
        mkdir -p temp/{wildcards.sm}/align logs/{wildcards.sm}/dedup && \
        pbmarkdup -j {threads} --end-length {params.dup_end_len}  --min-id-perc {params.min_id} --log-file {log} {input} {output.dedup_bam} && \
        samtools index {output.dedup_bam}
        """
    


# TODO add option for alignment and check if bam is aligned
# Output this as a CRAM file instead
rule align:
    input:
        fa= REF,
        data= "temp/{sm}/align/{sm}.dupmark.bam"
    conda:
        "../envs/cmd.yaml"
    output:
        aligned_bam="results/{sm}/align/{sm}.mapped.bam",
        index="results/{sm}/align/{sm}.mapped.bam.bai"
    threads: 8
    shell:
        """
        mkdir -p results/{wildcards.sm}/align && \
        RG=$(samtools view -H {input.data} | grep "^@RG" | sed 's/\t/\\\\t/g') && \
        minimap2 -t {threads} --MD --secondary=no -R "$RG" -Y -y -a -x map-pb  {input.fa} <(samtools fastq -T "*" {input.data}) |\
        samtools sort -u > {output.aligned_bam} && \
        samtools index {output.aligned_bam}
        """

rule sequence_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.bam"
    params:
        regions= get_input_regs
    output:
        read_metrics="results/{sm}/qc/{sm}.detailed_seq_metrics.tbl.gz",
        summary_metrics="results/{sm}/qc/{sm}.summary_seq_metrics.tbl.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sequence_metrics.py"

rule plot_seq_qc:
    input:
        summary_metrics = "results/{sm}/qc/{sm}.summary_seq_metrics.tbl.gz"
    params:
        region = "{region}"
    output:
        deam_rate = "results/{sm}/qc/{sm}.{region}.deam_rate.pdf",
        mut_rate = "results/{sm}/qc/{sm}.{region}.mut_rate.pdf",
        strandtype = "results/{sm}/qc/{sm}.{region}.strandtype.pdf",
        bias = "results/{sm}/qc/{sm}.{region}.bias.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_sequence_metrics.py"


rule targeting_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.bam"
    params:
        regions= get_input_regs
    output:
        detailed="results/{sm}/qc/{sm}.detailed_targeting_metrics.tbl.gz",
        summary="results/{sm}/qc/{sm}.summary_targeting_metrics.tbl"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/target_metrics.py"

rule plot_targeting_qc:
    input:
        targeting_metrics= "results/{sm}/qc/{sm}.summary_targeting_metrics.tbl"
    params:
        regions= get_input_regs
    output:
        plot= "results/{sm}/qc/{sm}.targeting_plot.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_targeting_metrics.py"



# Later, consensus generation. For now, this should not be included by default
# Consensus sequence bam generation- filter for just full-length and C>T/G>A designated strands. Use MSA to generate a consensus. Map consensus.
# Targeting metrics for consensus files
# Mutation rate metrics for consensus files
# add rule to save input parameters to smk folder