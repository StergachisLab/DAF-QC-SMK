# Basic rules
rule deduplicate:
    input:
        data= get_input_file
    params:
        dup_end_len=DUP_END_LENGTH,
        min_id=MIN_ID_PERC
    conda:
        "../envs/cmd.yaml"
    output:
        dedup_bam=temp("temp/{sm}/align/{sm}.reads.bam"),
        dup_index=temp("temp/{sm}/align/{sm}.reads.bam.bai")
    log:
        "logs/{sm}/dedup/{sm}.dedup.log"
    threads: 8
    benchmark:
        "benchmark/{sm}.benchmark.txt"
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
            data= "temp/{sm}/align/{sm}.{type}.bam"
        conda:
            "../envs/cmd.yaml"
        output:
            aligned_bam="results/{sm}/align/{sm}.mapped.{type}.bam",
            index="results/{sm}/align/{sm}.mapped.{type}.bam.bai"
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
            data=get_input_file
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
        regions= get_input_regs,
        end_tolerance=END_TOLERANCE
    output:
        detailed="results/{sm}/qc/reads/{sm}.detailed_targeting_metrics.tbl.gz",
        summary="results/{sm}/qc/reads/{sm}.summary_targeting_metrics.tbl"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/target_metrics.py"



rule plot_targeting_qc:
    input:
        targeting_metrics= "results/{sm}/qc/reads/{sm}.summary_targeting_metrics.tbl"
    params:
        regions= get_input_regs
    output:
        plot= "results/{sm}/qc/reads/plots/{sm}.targeting_plot.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_targeting_metrics.py"

rule sequence_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.{type}.bam",
        targeting_data=get_targeting_data
    params:
        regions= get_input_regs,
        chimera_cutoff = CHIMERA_CUTOFF,
        min_deamination_count = MIN_DEAMINATION_COUNT
    output:
        read_metrics="results/{sm}/qc/{type}/{sm}.detailed_seq_metrics.{type}.tbl.gz",
        summary_metrics="results/{sm}/qc/{type}/{sm}.summary_seq_metrics.{type}.tbl.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sequence_metrics.py"

rule plot_seq_qc:
    input:
        summary_metrics = "results/{sm}/qc/{type}/{sm}.summary_seq_metrics.{type}.tbl.gz"
    params:
        region = "{region}"
    output:
        deam_rate = "results/{sm}/qc/{type}/plots/{sm}.{region}.deam_rate.{type}.pdf",
        mut_rate = "results/{sm}/qc/{type}/plots/{sm}.{region}.mut_rate.{type}.pdf",
        strandtype = "results/{sm}/qc/{type}/plots/{sm}.{region}.strandtype.{type}.pdf",
        bias = "results/{sm}/qc/{type}/plots/{sm}.{region}.bias.{type}.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_sequence_metrics.py"



rule filter_bam:
    input:
        bam="results/{sm}/align/{sm}.mapped.reads.bam",
        seq_metrics="results/{sm}/qc/reads/{sm}.detailed_seq_metrics.reads.tbl.gz"
    params:
        sample_size=DECORATED_SAMPLESIZE
    output:
        filtered_bam="temp/{sm}/align/{sm}.filtered.bam",
        index="temp/{sm}/align/{sm}.filtered.bam.bai",
        sample_bam="temp/{sm}/align/{sm}.sample.bam",
        sample_index="temp/{sm}/align/{sm}.sample.bam.bai"
    conda:
        "../envs/cmd.yaml"
    shell:
        """
        samtools view -F 2306 -b -N <(zcat {input.seq_metrics} | awk -F'\t' '$8=="CT" || $8=="GA" {{print $1}}') {input.bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        read_count=$(samtools view -c {output.filtered_bam})
        read_frac=$(awk -v s={params.sample_size} -v b=$read_count 'BEGIN {{if (b>s) {{printf "%.6f\\n", s/b}} else {{printf "%.6f\\n", 1.0}}}}')
        if (( $(echo "$read_frac != 1.000000" | bc -l) )); then
            samtools view -b -s $read_frac {output.filtered_bam} > {output.sample_bam}
            samtools index {output.sample_bam}
        else
            ln -s $(realpath {output.filtered_bam}) {output.sample_bam}
            ln -s $(realpath {output.index}) {output.sample_index}
        fi
        """


rule decorate_strands:
    input:
        bam="temp/{sm}/align/{sm}.sample.bam",
        seq_metrics="results/{sm}/qc/reads/{sm}.detailed_seq_metrics.reads.tbl.gz"
    output:
        decorated_bam="results/{sm}/align/{sm}.decorated.reads.bam"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/decorate_strands.py"


rule index_decorated:
    input:
        "results/{sm}/align/{sm}.decorated.reads.bam"
    output:
        "results/{sm}/align/{sm}.decorated.reads.bam.bai"
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
        deduplication_metrics="results/{sm}/qc/reads/{sm}.deduplication_metrics.tbl.gz"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/deduplication_metrics.py"



rule plot_deduplication_metrics:
    input:
        summary_metrics="results/{sm}/qc/reads/{sm}.deduplication_metrics.tbl.gz"
    params:
        region = "{region}",
        consensus_min_reads = CONSENSUS_MIN_READS
    output:
        duplication_groups = "results/{sm}/qc/reads/plots/{sm}.{region}.duplication_groups.pdf"
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
        bam=temp("temp/{sm}/align/{sm}.consensus.bam")
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_consensus.py"


rule make_dashboard:
    input:
        pdfs=get_qc_plot_names
    params:
        sample_name="{sm}",
        regions=get_input_regs
    output:
        dashboard="results/{sm}/qc/{sm}.dashboard.html"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_dashboard.py"




# Targeting metrics for consensus files?
# add rule to save input parameters to smk folder
# add logs where needed
# fix start end coordinates in strand metrics to clarify fiber coordinates vs target coordinates
# As an optional note, can add script on github to decorate the consensus and merge it with the filtered file.
# add total coverage of targets
# add option to require a certain amount of CT or GA to call a strand
# add proportion or count of reads that are above cutoff on deduplication plot
# subsample for decorated and sequencing metrics? 