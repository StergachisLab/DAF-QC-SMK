def get_ref():
    if "ref" not in config:
        raise ValueError("FIRE: ref parameter is missing in config.yaml")
    ref = config["ref"]
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)


def get_manifest():
    manifest = config.get("manifest")
    if manifest is None:
        raise ValueError("manifest parameter is missing in config.yaml")
    if not os.path.exists(manifest):
        raise ValueError(f"Manifest file {manifest} does not exist")
    manifest = pd.read_csv(config["manifest"], sep=r"\s+", comment="#").set_index(
        "sample"
    )
    return manifest


def get_dup_end_length():
    if "dup_end_length" not in config:
        dup_end_length = 0
    else:
        dup_end_length = config["dup_end_length"]
    return dup_end_length


def get_min_id_perc():
    if "min_id_perc" not in config:
        min_id_perc = 99.2
    else:
        min_id_perc = config["min_id_perc"]
    return min_id_perc


def get_chimera_cutoff():
    if "chimera_cutoff" not in config:
        chimera_cutoff = 0.9
    else:
        chimera_cutoff = config["chimera_cutoff"]
    return chimera_cutoff


def get_min_deamination_count():
    if "min_deamination_count" not in config:
        min_deamination_count = 50
    else:
        min_deamination_count = config["min_deamination_count"]
    return min_deamination_count


def get_end_tolerance():
    if "end_tolerance" not in config:
        end_tolerance = 30
    else:
        end_tolerance = config["end_tolerance"]
    return end_tolerance


def get_consensus_min_reads():
    if "consensus_min_reads" not in config:
        consensus_min_reads = 3
    else:
        consensus_min_reads = config["consensus_min_reads"]
    return consensus_min_reads


def get_platform():
    if "platform" not in config:
        platform = "pacbio"
    else:
        platform = config["platform"]
    return platform


def get_consensus():
    if "consensus" not in config:
        consensus = True
    else:
        consensus = config["consensus"]
    return consensus


def get_is_fastq():
    if "is_fastq" not in config:
        is_fastq = False
    else:
        is_fastq = config["is_fastq"]
    return is_fastq


def get_decorated_samplesize():
    if "decorated_samplesize" not in config:
        sample_size = 5000
    else:
        sample_size = config["decorated_samplesize"]
    return sample_size


def get_input_file(wc):
    return MANIFEST.loc[wc.sm, "file"]


def get_input_regs(wc):
    regions = MANIFEST.loc[wc.sm, "regs"].strip().split(",")
    return regions


def get_strand_qc_plot_paths():
    outputs = []
    types = ["reads"] if CONSENSUS is False else ["reads", "consensus"]
    for sample in MANIFEST.index:
        regions = MANIFEST.loc[sample, "regs"].strip().split(",")
        for region in regions:
            reg_label = region.replace(":", "_").replace("-", "_")
            prefix = f"results/{sample}/qc"
            for plot in ["deam_rate", "bias", "mut_rate", "strandtype"]:
                for item in types:
                    outputs.append(
                        f"{prefix}/{item}/plots/{sample}.{reg_label}.{plot}.{item}.pdf"
                    )
    return outputs


def get_deduplication_qc_plot_paths():
    outputs = []
    for sample in MANIFEST.index:
        regions = MANIFEST.loc[sample, "regs"].strip().split(",")
        for region in regions:
            reg_label = region.replace(":", "_").replace("-", "_")
            prefix = f"results/{sample}/qc/reads/plots/{sample}.{reg_label}"
            outputs.append(f"{prefix}.duplication_groups.pdf")
    return outputs


def get_targeting_data(wc):
    if wc.type == "reads":
        return f"results/{wc.sm}/qc/reads/{wc.sm}.detailed_targeting_metrics.tbl.gz"
    else:
        return []


def get_qc_plot_names(wc):
    types = ["reads"] if CONSENSUS is False else ["reads", "consensus"]
    plots = []

    sample = wc.sm
    regions = MANIFEST.loc[sample, "regs"].strip().split(",")

    plots.append(f"results/{sample}/qc/reads/plots/{sample}.targeting_plot.pdf")

    for region in regions:
        reg_label = region.replace(":", "_").replace("-", "_")
        prefix = f"results/{sample}/qc"

        for plot in ["deam_rate", "bias", "mut_rate", "strandtype"]:
            for item in types:
                plots.append(
                    f"{prefix}/{item}/plots/{sample}.{reg_label}.{plot}.{item}.pdf"
                )

        if PLATFORM != "ont":
            plots.append(
                f"{prefix}/reads/plots/{sample}.{reg_label}.duplication_groups.pdf"
            )

    return plots
