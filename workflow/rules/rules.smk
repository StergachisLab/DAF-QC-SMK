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

# Targeting metrics- per target % of reads, all targets % of reads, metrics for full length reads
# Mutation rate metrics, output C,G, or other strand designation
# Deduplication
# Targeting metrics for deduplicated files
# Mutation rate metrics for deduplicated files

# Later, consensus generation. For now, this should not be included by default
# Consensus sequence bam generation- filter for just full-length and C>T/G>A designated strands. Use MSA to generate a consensus. Map consensus. 