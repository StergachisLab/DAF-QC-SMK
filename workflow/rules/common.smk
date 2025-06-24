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
        dup_end_length = 8000
    else:
        dup_end_length=config["dup_end_length"]
    return dup_end_length

def get_min_id_perc():
    if "min_id_perc" not in config:
        min_id_perc = 99.2
    else:
        min_id_perc=config["min_id_perc"]
    return min_id_perc

def get_chimera_cutoff():
    if "chimera_cutoff" not in config:
        chimera_cutoff = 0.9
    else:
        chimera_cutoff=config["chimera_cutoff"]
    return chimera_cutoff

def get_input_bam(wc):
    return MANIFEST.loc[wc.sm, "bam"]

def get_input_regs(wc):
    regions = MANIFEST.loc[wc.sm, "regs"].strip().split(",")
    return regions

def get_strand_qc_plot_paths():
    outputs = []
    for sample in MANIFEST.index:
        regions = MANIFEST.loc[sample, "regs"].strip().split(",")
        for region in regions:
            reg_label = region.replace(":","_").replace("-", "_")
            prefix = f"results/{sample}/qc/{sample}.{reg_label}"
            for plot in ["deam_rate", "bias", "mut_rate", "strandtype"]:
                for type in ["reads"]:
                    outputs.append(f"{prefix}.{plot}.{type}.pdf")
    return outputs

def get_deduplication_qc_plot_paths():
    outputs = []
    for sample in MANIFEST.index:
        regions = MANIFEST.loc[sample, "regs"].strip().split(",")
        for region in regions:
            reg_label = region.replace(":","_").replace("-", "_")
            prefix = f"results/{sample}/qc/{sample}.{reg_label}"
            for level in ["reads", "groups"]:
                outputs.append(f"{prefix}.duplication_{level}.pdf")
    return outputs