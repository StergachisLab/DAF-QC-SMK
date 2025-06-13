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

def get_input_bam(wc):
    return MANIFEST.loc[wc.sm, "bam"]
