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
    

if PLATFORM == "pacbio":
# TODO add option for alignment and check if bam is aligned
# Output this as a CRAM file instead
    rule align:
        input:
            fa= REF,
            data= "temp/{sm}/align/{sm}.dupmark.bam"
        conda:
            "../envs/cmd.yaml"
        output:
            aligned_bam="results/{sm}/align/{sm}.mapped.reads.bam",
            index="results/{sm}/align/{sm}.mapped.reads.bam.bai"
        threads: 8
        shell:
            """
            mkdir -p results/{wildcards.sm}/align && \
            minimap2 -t {threads} --MD --secondary=no -Y -y -a -x map-pb  {input.fa} <(samtools fastq -T "*" {input.data}) |\
            samtools sort > {output.aligned_bam} && \
            samtools index {output.aligned_bam}
            """

elif PLATFORM == "ont":

    rule align_ont:
        input:
            fa=REF,
            data=get_input_bam
        conda:
            "../envs/cmd.yaml"
        output:
            aligned_bam="results/{sm}/align/{sm}.mapped.reads.bam",
            index="results/{sm}/align/{sm}.mapped.reads.bam.bai"
        threads: 8
        shell:
            """
            mkdir -p results/{wildcards.sm}/align && \
            if [[ "{IS_FASTQ}" == "True" ]]; then
                minimap2 -t {threads} --MD --secondary=no -Y -y -a -x map-ont  {input.fa} {input.data} |\
                samtools sort > {output.aligned_bam}
            else
                minimap2 -t {threads} --MD --secondary=no -Y -y -a -x map-ont  {input.fa} <(samtools fastq -T "*" {input.data}) |\
                samtools sort > {output.aligned_bam} 
            fi && \
            samtools index {output.aligned_bam}
            """

rule targeting_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.reads.bam"
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

rule sequence_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.{type}.bam",
        targeting_data="results/{sm}/qc/{sm}.detailed_targeting_metrics.tbl.gz"
    params:
        regions= get_input_regs,
        chimera_cutoff = CHIMERA_CUTOFF
    output:
        read_metrics="results/{sm}/qc/{sm}.detailed_seq_metrics.{type}.tbl.gz",
        summary_metrics="results/{sm}/qc/{sm}.summary_seq_metrics.{type}.tbl.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sequence_metrics.py"

rule plot_seq_qc:
    input:
        summary_metrics = "results/{sm}/qc/{sm}.summary_seq_metrics.{type}.tbl.gz"
    params:
        region = "{region}"
    output:
        deam_rate = "results/{sm}/qc/{sm}.{region}.deam_rate.{type}.pdf",
        mut_rate = "results/{sm}/qc/{sm}.{region}.mut_rate.{type}.pdf",
        strandtype = "results/{sm}/qc/{sm}.{region}.strandtype.{type}.pdf",
        bias = "results/{sm}/qc/{sm}.{region}.bias.{type}.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_sequence_metrics.py"

rule filter_bam:
    input:
        bam="results/{sm}/align/{sm}.mapped.reads.bam",
        seq_metrics="results/{sm}/qc/{sm}.detailed_seq_metrics.reads.tbl.gz"
    output:
        filtered_bam="temp/{sm}/align/{sm}.filtered.bam",
        index="temp/{sm}/align/{sm}.filtered.bam.bai"
    conda:
        "../envs/cmd.yaml"
    shell:
        """
        samtools view -F 2306 -b -N <(zcat {input.seq_metrics} | awk -F'\t' '$6=="CT" || $6=="GA" {{print $1}}') {input.bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """


rule decorate_strands:
    input:
        bam="temp/{sm}/align/{sm}.filtered.bam",
        seq_metrics="results/{sm}/qc/{sm}.detailed_seq_metrics.reads.tbl.gz"
    output:
        decorated_bam="results/{sm}/align/{sm}.decorated.bam"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/decorate_strands.py"

rule index_decorated:
    input:
        "results/{sm}/align/{sm}.decorated.bam"
    output:
        "results/{sm}/align/{sm}.decorated.bam.bai"
    conda:
        "../envs/cmd.yaml"
    shell:
        """
        samtools index {input}
        """



rule deduplication_metrics:
    input:
       bam="temp/{sm}/align/{sm}.filtered.bam"
    params:
        regions= get_input_regs
    output:
        deduplication_metrics="results/{sm}/qc/{sm}.deduplication_metrics.tbl.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/deduplication_metrics.py"



rule plot_deduplication_metrics:
    input:
        summary_metrics="results/{sm}/qc/{sm}.deduplication_metrics.tbl.gz"
    params:
        region = "{region}"
    output:
        duplication_reads = "results/{sm}/qc/{sm}.{region}.duplication_reads.pdf",
        duplication_groups = "results/{sm}/qc/{sm}.{region}.duplication_groups.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_deduplication_metrics.py"



rule build_consensus:
    input:
        bam="temp/{sm}/align/{sm}.filtered.bam"
    params:
        consensus_min_reads=CONSENSUS_MIN_READS
    output:
        bam="temp/{sm}/align/{sm}.consensus.bam"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_consensus.py"

# output as temp file, needs to be aligned

# rule align consensus
# output is aligned consensus.

# Once consensus functionality is added, add to all rule

# As an optional note, can add script on github to decorate the consensus and merge it with the filtered file.

# Later, consensus generation. For now, this should not be included by default
# Consensus sequence bam generation- filter for just full-length and C>T/G>A designated strands. Use MSA to generate a consensus. Map consensus.
# Targeting metrics for consensus files
# Mutation rate metrics for consensus files
# add rule to save input parameters to smk folder
# add logs where needed
# remove read group line from align rule, make align universal for reads and consensus
# make targeting plot % on targets more clear, fix percent vs fraction issue.