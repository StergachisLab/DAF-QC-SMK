rule filter_primary:
    input:
        "results/{sm}/align/{sm}.mapped.bam"
    conda:
        "../envs/cmd.yaml"
    output:
        prim_bam=temp("temp/{sm}/dedup/{sm}.primary.bam"),
        prim_index=temp("temp/{sm}/dedup/{sm}.primary.bam.bai")
    threads: 8
    shell:
        """
        mkdir -p temp/{wildcards.sm}/dedup && \
        samtools view -b -F 2306 {input} > {output.prim_bam} && \
        samtools index {output.prim_bam}
        """
