import pysam
import pandas as pd
import os


def parse_region(region):
    """
    Parses a genomic region in the format chr:start-end.
    Returns a tuple of chromosome, start, and end positions.
    """
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    return chrom, start, end


def strand_metrics(read, cutoff=0.9, min_deamination_count=50):
    """
    Determines the strand acted upon by DddA based on the proportion of C->T & G->A,
    and counts deamination for each two base pair context (AC, CC, GC, TC).
    For G->A deamination strands, the complement is recorded.

    Args:
        read (pysam.AlignedSegment): A read from a BAM file.
        cutoff (float): Proportion threshold to determine strand type.

    Returns a tuple containing:
        - strand: "CT", "GA", "chimera", "undetermined", or "none"
        - doublets: Dictionary with counts of deamination doublets
        - mutation_count: Number of mutations in the read
        - deamination_pos: List of positions in the read where deamination occurred
    """
    # Determine strand
    ct = 0
    ga = 0
    non_ref = 0

    seq = read.query_sequence
    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)

    for pos in pair:
        if pos[0] is None or pos[1] is None:  # indel, ignore
            continue
        else:
            strand_base = seq[pos[0]]
            ref_base = pos[2].upper()
            if strand_base != ref_base:
                non_ref += 1
                if ref_base == "C" and strand_base == "T":
                    ct += 1
                elif ref_base == "G" and strand_base == "A":
                    ga += 1

    if ct + ga == 0:
        strand = "none"
    elif ct + ga / non_ref <= cutoff or ct + ga < min_deamination_count:
        strand = "undetermined"
    elif ct / (ct + ga) >= cutoff:
        strand = "CT"
    elif ga / (ct + ga) >= cutoff:
        strand = "GA"
    else:
        strand = "chimera"

    # Count and record position of deaminations for CT and GA strands
    if strand in ["CT", "GA"]:
        mutation_count = non_ref - ct if strand == "CT" else non_ref - ga

        deamination_pos = []  # deaminated positions in read coordinates
        doublets = {
            "AC": [0, 0],
            "CC": [0, 0],
            "GC": [0, 0],
            "TC": [0, 0],
            "OC": [0, 0], # OC represents no base/indel before C
        }

        doublet_dict = {"GA": "TC", "GC": "GC", "GG": "CC", "GT": "AC"}

        for i, pos in enumerate(pair):
            if pos[0] is None or pos[1] is None:  # indel, ignore
                continue

            ref_base = pos[2].upper()

            if ref_base != strand[0]:  # ignore non C/G bases
                continue

            strand_base = seq[pos[0]]

            if strand == "CT":
                doublet = (
                    pair[i - 1][2].upper() + ref_base
                    if pair[i - 1][2] is not None and i > 0
                    else "OC"
                )
            else:
                doublet = (
                    doublet_dict[ref_base + pair[i + 1][2].upper()]
                    if i + 1 < len(pair) and pair[i + 1][2] is not None
                    else "OC"
                )

            # Record total counts for each pair type
            doublets[doublet][0] += 1

            # Record deamination counts and positions
            if strand_base != ref_base and strand_base == strand[1]:
                deamination_pos.append(pos[0])
                doublets[doublet][1] += 1

    else:
        doublets = {
            "AC": [None, None],
            "CC": [None, None],
            "GC": [None, None],
            "TC": [None, None],
            "OC": [None, None],
        }  # OC represents no base/indel before C
        mutation_count = None
        deamination_pos = None

    return strand, doublets, mutation_count, deamination_pos


def strand_metrics_table(
    psfile, chrom, start, end, cutoff=0.9, min_deamination_count=50, include_readnames=None
):
    """
    Iterates through reads in a specified region of a BAM file,
    calculates strand metrics, and returns a DataFrame with the results.

    Args:
        psfile (pysam.AlignmentFile): The BAM file to read.
        chrom (str): Chromosome name.
        start (int): Start position of the region.
        end (int): End position of the region.
        chimera_cutoff (float): Proportion threshold to determine strand type.
        min_deamination_count (int): Minimum number of deaminations to designate a strand.
        include_readnames (list, optional): List of read names to include in the analysis.
            If None, all reads in the region are included.
    """

    read_collector = []

    # check that this region has reads
    if psfile.count(chrom, start, end) == 0:
        return pd.DataFrame()

    for read in psfile.fetch(chrom, start, end):
        if include_readnames is not None and read.query_name not in include_readnames:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        strand, doublets, mutation_count, deam_pos = strand_metrics(
            read, cutoff=cutoff, min_deamination_count=min_deamination_count
        )

        duplicate = read.get_tag("du") if read.has_tag("du") else "None"

        read_data = {
            "read_name": read.query_name,
            "chrom": chrom,
            "reg_start": start,
            "reg_end": end,
            "strand_start": read.reference_start,
            "strand_end": read.reference_end,
            "length": len(read.query_sequence),
            "strand": strand,
            "duplicate": duplicate,
            "mutation_count": mutation_count,
            "deamination_positions": ",".join(map(str, deam_pos)) if deam_pos is not None else ""
        }

        keys = ["AC", "CC", "GC", "TC", "OC"]

        if strand in ["CT", "GA"]:
            for key in keys:
                read_data[f"{key}_count"] = doublets[key][0]
                read_data[f"{key}_deam"] = doublets[key][1]
                read_data[f"{key}_deam_rate"] = (
                    doublets[key][1] / doublets[key][0]
                    if doublets[key][0] > 0
                    else None
                )

        read_collector.append(read_data)

    reads_table = pd.DataFrame(read_collector)
    
    reads_table["total_count"] = reads_table[[f"{key}_count" for key in keys]].sum(axis=1)
    reads_table["total_deam"] = reads_table[[f"{key}_deam" for key in keys]].sum(axis=1)

    reads_table["all_deam_rate"] = (
        reads_table["total_deam"] / reads_table["total_count"]
    ).where(reads_table["total_count"] > 0)
    reads_table["mutation_rate"] = (
        reads_table["mutation_count"] / reads_table["length"]
    ).where(reads_table["length"] > 0)

    return reads_table


def aggregate_strand_metrics(table):
    """
    Aggregates single fiber strand metrics from a DataFrame by chromosome, start, end, and strand designation.
    """

    aggregate_collector = []

    for group in table.groupby(["chrom", "reg_start", "reg_end", "strand"]):
        mutation_rate = group[1]["mutation_rate"].dropna().tolist()
        all_deam_rate = group[1]["all_deam_rate"].dropna().tolist()
        AC_deam_rates = group[1]["AC_deam_rate"].dropna().tolist()
        CC_deam_rates = group[1]["CC_deam_rate"].dropna().tolist()
        GC_deam_rates = group[1]["GC_deam_rate"].dropna().tolist()
        TC_deam_rates = group[1]["TC_deam_rate"].dropna().tolist()
        OC_deam_rates = group[1]["OC_deam_rate"].dropna().tolist()
        count = len(group[1])

        aggregate_template = {
            "chrom": group[0][0],
            "reg_start": group[0][1],
            "reg_end": group[0][2],
            "strand": group[0][3],
            "count": count,
            "mutation_rate": ",".join(map(str, mutation_rate)) if mutation_rate else "",
            "all_deam_rate": ",".join(map(str, all_deam_rate)) if all_deam_rate else "",
            "AC_deam_rate": ",".join(map(str, AC_deam_rates)) if AC_deam_rates else "",
            "CC_deam_rate": ",".join(map(str, CC_deam_rates)) if CC_deam_rates else "",
            "GC_deam_rate": ",".join(map(str, GC_deam_rates)) if GC_deam_rates else "",
            "TC_deam_rate": ",".join(map(str, TC_deam_rates)) if TC_deam_rates else "",
            "OC_deam_rate": ",".join(map(str, OC_deam_rates)) if OC_deam_rates else "",
        }

        aggregate_collector.append(aggregate_template)

        aggregate_table = pd.DataFrame(aggregate_collector)

    return aggregate_table


## Main script execution
bam_path = snakemake.input.data
regions = snakemake.params.regions
targeting_metrics = snakemake.input.targeting_data
chimera_cutoff = snakemake.params.chimera_cutoff
min_deamination_count = snakemake.params.min_deamination_count
read_metrics = snakemake.output.read_metrics
summary_metrics = snakemake.output.summary_metrics
threads = snakemake.threads

# for testing
# bam_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/align/htt_test.mapped.reads.bam"
# regions = ["chr4:3073138-3075853", "chr3:179228176-179236561"]
# targeting_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test.detailed_targeting_metrics.tbl.gz"
# chimera_cutoff = 0.9
# min_deamination_count = 50 # minimum number of deaminations to designate a strand (otherwise undetermined)
# read_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.read_metrics_optimized.tsv.gz"
# summary_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.summary_metrics_optimized.tsv.gz"

os.makedirs(os.path.dirname(read_metrics), exist_ok=True)

pybam = pysam.AlignmentFile(bam_path, "rb", threads=threads)
tables = []

if len(targeting_metrics) > 0:
    targeting_df = pd.read_csv(targeting_metrics, sep="\t", compression="gzip")

    if targeting_df["full_length_reads"].isnull().all():
        raise ValueError(
            "No full length reads found in targeting metrics file. \n"
            "Check that the correct targeted regions are specified in your manifest table.\n"
            "Set end_tolerance to a larger value in your config file if the regions are not exact matches to the BAM file."
        )

    for region in regions:
        chrom, start, end = parse_region(region)

        full_length_reads = targeting_df[
            (targeting_df["chrom"] == chrom)
            & (targeting_df["start"] == start)
            & (targeting_df["end"] == end)
        ]["full_length_reads"].values

        if (
            full_length_reads is None
            or len(full_length_reads) == 0
            or pd.isna(full_length_reads[0])
        ):
            continue

        full_length_reads = full_length_reads[0].split(",")

        reg_table = strand_metrics_table(
            pybam,
            chrom,
            start,
            end,
            cutoff=chimera_cutoff,
            min_deamination_count=min_deamination_count,
            include_readnames=full_length_reads,
        )

        tables.append(reg_table)

else:
    for region in regions:
        chrom, start, end = parse_region(region)
        reg_table = strand_metrics_table(
            pybam, chrom, start, end, cutoff=chimera_cutoff, min_deamination_count=min_deamination_count
        )
        tables.append(reg_table)

table = pd.concat(tables, ignore_index=True)

aggregate_table = aggregate_strand_metrics(table)

table.to_csv(read_metrics, sep="\t", index=False, compression="gzip")
aggregate_table.to_csv(summary_metrics, sep="\t", index=False, compression="gzip")