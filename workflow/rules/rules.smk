# Basic rules
rule deduplicate:
    input:
        data=get_input_file,
    params:
        dup_end_len=DUP_END_LENGTH,
        min_id=MIN_ID_PERC,
    conda:
        "../envs/cmd.yaml"
    output:
        dedup_bam="temp/{sm}/align/{sm}.reads.bam",
        dup_index="temp/{sm}/align/{sm}.reads.bam.bai",
        pb_index="temp/{sm}/align/{sm}.reads.bam.pbi",
    resources:
        mem_mb=20 * 1024,
        runtime=1440,
    log:
        "logs/{sm}/dedup/{sm}.dedup.log",
    threads: 40
    benchmark:
        repeat("benchmark/{sm}.dedup_benchmark.txt", 1 if BENCHMARK else 0)
    shell:
        """
        mkdir -p temp/{wildcards.sm}/align logs/{wildcards.sm}/dedup && \
        pbmarkdup -j {threads} --end-length {params.dup_end_len}  --min-id-perc {params.min_id} --log-file {log} {input.data} {output.dedup_bam} && \
        samtools index {output.dedup_bam}
     """


if PLATFORM == "pacbio":

    # TODO add option for alignment and check if bam is aligned
    # Output this as a CRAM file instead
    rule align:
        input:
            fa=REF,
            data="temp/{sm}/align/{sm}.{type}.bam",
        conda:
            "../envs/cmd.yaml"
        output:
            aligned_bam="results/{sm}/align/{sm}.mapped.{type}.bam",
            index="results/{sm}/align/{sm}.mapped.{type}.bam.bai",
        threads: 16
        resources:
            mem_mb=20 * 1024,
            runtime=60 * 4,
        benchmark:
            repeat("benchmark/{sm}.align_benchmark.{type}.txt", 1 if BENCHMARK else 0)
        shell:
            """
            mkdir -p results/{wildcards.sm}/align && \
            minimap2 -t {threads} --MD --secondary=no -Y -y -a -x map-pb  {input.fa} <(samtools fastq -T "*" {input.data}) |\
            rb add-rg -t 8 -u {input.data} |\
            samtools sort > {output.aligned_bam} && \
            samtools index {output.aligned_bam}
            """

elif PLATFORM == "ont":

    rule align_ont:
        input:
            fa=REF,
            data=get_input_file,
        conda:
            "../envs/cmd.yaml"
        output:
            aligned_bam="results/{sm}/align/{sm}.mapped.reads.bam",
            index="results/{sm}/align/{sm}.mapped.reads.bam.bai",
        threads: 16
        resources:
            mem_mb=20 * 1024,
            runtime=60 * 4,
        benchmark:
            repeat("benchmark/{sm}.align_benchmark_ont.txt", 1 if BENCHMARK else 0)
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
        data="results/{sm}/align/{sm}.mapped.reads.bam",
    params:
        regions=get_input_regs,
        end_tolerance=END_TOLERANCE,
    output:
        detailed="results/{sm}/qc/reads/{sm}.detailed_targeting_metrics.tbl.gz",
        summary="results/{sm}/qc/reads/{sm}.summary_targeting_metrics.tbl",
    threads: 8
    benchmark:
        repeat("benchmark/{sm}.targetqc_benchmark.txt", 1 if BENCHMARK else 0)
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/target_metrics.py"


rule plot_targeting_qc:
    input:
        targeting_metrics="results/{sm}/qc/reads/{sm}.summary_targeting_metrics.tbl",
    params:
        regions=get_input_regs,
    output:
        plot="results/{sm}/qc/reads/plots/{sm}.targeting_plot.pdf",
        metrics_txt="results/{sm}/qc/reads/plots/{sm}.targeting_plot.txt",
    benchmark:
        repeat("benchmark/{sm}.targetplot_benchmark.txt", 1 if BENCHMARK else 0)
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_targeting_metrics.py"


rule sequence_qc:
    input:
        data="results/{sm}/align/{sm}.mapped.{type}.bam",
        targeting_data=get_targeting_data,
    params:
        regions=get_input_regs,
        chimera_cutoff=CHIMERA_CUTOFF,
        min_deamination_count=MIN_DEAMINATION_COUNT,
    output:
        read_metrics="results/{sm}/qc/{type}/{sm}.detailed_seq_metrics.{type}.tbl.gz",
        summary_metrics="results/{sm}/qc/{type}/{sm}.summary_seq_metrics.{type}.tbl.gz",
    threads: 8
    resources:
        mem_mb=20 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.seqqc_benchmark.{type}.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/sequence_metrics.py"


rule plot_seq_qc:
    input:
        summary_metrics="results/{sm}/qc/{type}/{sm}.summary_seq_metrics.{type}.tbl.gz",
    params:
        region="{region}",
    output:
        deam_rate="results/{sm}/qc/{type}/plots/{sm}.{region}.deam_rate.{type}.pdf",
        mut_rate="results/{sm}/qc/{type}/plots/{sm}.{region}.mut_rate.{type}.pdf",
        strandtype="results/{sm}/qc/{type}/plots/{sm}.{region}.strandtype.{type}.pdf",
        bias="results/{sm}/qc/{type}/plots/{sm}.{region}.bias.{type}.pdf",
        deam_rate_table="results/{sm}/qc/{type}/plots/{sm}.{region}.deam_rate.{type}.txt",  # for testing
        mut_rate_table="results/{sm}/qc/{type}/plots/{sm}.{region}.mut_rate.{type}.txt",
        bias_table="results/{sm}/qc/{type}/plots/{sm}.{region}.bias.{type}.txt",
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.seqplot_benchmark.{type}.{region}.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/plot_sequence_metrics.py"


rule filter_bam:
    input:
        bam="results/{sm}/align/{sm}.mapped.reads.bam",
        seq_metrics="results/{sm}/qc/reads/{sm}.detailed_seq_metrics.reads.tbl.gz",
    params:
        sample_size=DECORATED_SAMPLESIZE,
    output:
        filtered_bam=temp("temp/{sm}/align/{sm}.filtered.bam"),
        index=temp("temp/{sm}/align/{sm}.filtered.bam.bai"),
        sample_bam=temp("temp/{sm}/align/{sm}.sample.bam"),
        sample_index=temp("temp/{sm}/align/{sm}.sample.bam.bai"),
    threads: 8
    conda:
        "../envs/cmd.yaml"
    benchmark:
        repeat("benchmark/{sm}.filterbam_benchmark.txt", 1 if BENCHMARK else 0)
    shell:
        """
        samtools view -@ {threads} -F 2306 -b -N <(zcat {input.seq_metrics} | awk -F'\t' '$8=="CT" || $8=="GA" {{print $1}}') {input.bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        read_count=$(samtools view -@ {threads} -c {output.filtered_bam})
        read_frac=$(awk -v s={params.sample_size} -v b=$read_count 'BEGIN {{if (b>s) {{printf "%.6f\\n", s/b}} else {{printf "%.6f\\n", 1.0}}}}')
        if (( $(echo "$read_frac != 1.000000" | bc -l) )); then
            samtools view -@ {threads} -b -s $read_frac {output.filtered_bam} > {output.sample_bam}
            samtools index {output.sample_bam}
        else
            ln -s $(realpath {output.filtered_bam}) {output.sample_bam}
            ln -s $(realpath {output.index}) {output.sample_index}
        fi
        """


rule decorate_strands:
    input:
        bam="temp/{sm}/align/{sm}.sample.bam",
        bai="temp/{sm}/align/{sm}.sample.bam.bai",
        seq_metrics="results/{sm}/qc/reads/{sm}.detailed_seq_metrics.reads.tbl.gz",
        filtered_bam="temp/{sm}/align/{sm}.filtered.bam",
        filtered_bai="temp/{sm}/align/{sm}.filtered.bam.bai",
    output:
        decorated_bam="results/{sm}/align/{sm}.decorated.reads.bam",
    threads: 8
    resources:
        mem_mb=20 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.decorate_benchmark.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/decorate_strands.py"


rule index_decorated:
    input:
        "results/{sm}/align/{sm}.decorated.reads.bam",
    output:
        "results/{sm}/align/{sm}.decorated.reads.bam.bai",
    conda:
        "../envs/cmd.yaml"
    shell:
        """
        samtools index {input}
        """


rule deduplication_metrics:
    input:
        bam="temp/{sm}/align/{sm}.filtered.bam",
        bai="temp/{sm}/align/{sm}.filtered.bam.bai",
    params:
        regions=get_input_regs,
    output:
        deduplication_metrics="results/{sm}/qc/reads/{sm}.deduplication_metrics.tbl.gz",
    threads: 8
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.dedupqc_benchmark.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/deduplication_metrics.py"


rule plot_deduplication_metrics:
    input:
        summary_metrics="results/{sm}/qc/reads/{sm}.deduplication_metrics.tbl.gz",
    params:
        region="{region}",
        consensus_min_reads=CONSENSUS_MIN_READS,
    output:
        duplication_groups="results/{sm}/qc/reads/plots/{sm}.{region}.duplication_groups.pdf",
        duplication_table="results/{sm}/qc/reads/plots/{sm}.{region}.duplication_groups.txt",
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.dedupplot_benchmark.{region}.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/plot_deduplication_metrics.py"


rule build_consensus:
    input:
        bam="temp/{sm}/align/{sm}.filtered.bam",
        bai="temp/{sm}/align/{sm}.filtered.bam.bai",
    params:
        consensus_min_reads=CONSENSUS_MIN_READS,
    output:
        bam=temp("temp/{sm}/align/{sm}.consensus.bam"),
    threads: 8
    resources:
        mem_mb=20 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.consensus_benchmark.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/build_consensus.py"


rule make_dashboard:
    input:
        pdfs=get_qc_plot_names,
    params:
        sample_name="{sm}",
        regions=get_input_regs,
    output:
        dashboard="results/{sm}/qc/{sm}.dashboard.html",
    conda:
        "../envs/python.yaml"
    benchmark:
        repeat("benchmark/{sm}.dashboard_benchmark.txt", 1 if BENCHMARK else 0)
    script:
        "../scripts/create_dashboard.py"
