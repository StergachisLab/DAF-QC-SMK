# Basic rules
# TODO add option for alignment and check if bam is aligned
rule align:
    input:
        fa= REF,
        data= get_input_bam
    conda:
        "../envs/cmd.yaml"
    output:
        aligned_bam="results/{sm}/align/{sm}.mapped.bam",
        index="results/{sm}/align/{sm}.mapped.bam.bai"
    threads: 8
    shell:
        """
        mkdir -p results/{wildcards.sm}/align && \
        minimap2 -t {threads} --MD -N 0 -Y -y -a -x map-pb {input.fa} <(samtools fastq -T "*" {input.data}) | samtools sort -u > {output.aligned_bam} && \
        samtools index {output.aligned_bam}
        """

rule deduplicate
    input:
        "results/{sm}/align/{sm}.mapped.bam"
    conda:
        "../envs/cmd.yaml"
    output:
        dedup_bam="results/{sm}/align/{sm}.dupmark.bam"
        dup_index="results/{sm}/align/{sm}.dupmark.bam"
    threads: 8
    shell:
    """
    mkdir -p results/{wildcards.sm}/align && \
    pbmarkdup -j {threads} 
    


# Targeting metrics- per target % of reads, all targets % of reads, metrics for full length reads
# Mutation rate metrics, output C,G, or other strand designation
# Deduplication
# Targeting metrics for deduplicated files
# Mutation rate metrics for deduplicated files

# Later, consensus generation. For now, this should not be included by default
# Consensus sequence bam generation- filter for just full-length and C>T/G>A designated strands. Use MSA to generate a consensus. Map consensus. 


# first step: get align option working. Get it to take file from table.