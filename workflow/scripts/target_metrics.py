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


def region_metrics_table(bam, chrom, start, end, tolerance=30):
    """
    Analyze reads for full-length alignments

    Filters for reads that align at the expected ends (+/- tolerace)
    and requires that reads do not have more than 30 bp of soft clipping.
    Filters out non-primary alignments.

    Args:
        bam (pysam.AlignmentFile): BAM file to analyze.
        chrom (str): Chromosome name.
        start (int): Start position of the region.
        end (int): End position of the region.
        tolerance (int): Tolerance for start and end positions (bp), default is 30.

    Returns:
        Dataframe containing the following columns:
            - chrom: Chromosome name.
            - start: Start position of the region.
            - end: End position of the region.
            - full_length_reads: List of read names that are full-length.
            - non_full_length_reads: List of read names that overlap the region but are not full-length (either due to mapping or soft clipping).
            - non_primary_reads: List of read names that overlap the region but are not primary alignments.
    """

    full_length_reads = []
    non_full_length_reads_map = []
    non_full_length_reads_clip = []
    non_primary_reads = []

    for read in bam.fetch(chrom, start, end):
        # skip reads that are unmapped, secondary, or supplementary
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            non_primary_reads.append(read.query_name)
            continue

        map_start, map_end = read.reference_start, read.reference_end

        # check for full length alignments
        map_start_diff = map_start - start
        map_end_diff = map_end - end

        if abs(map_start_diff) > tolerance or abs(map_end_diff) > tolerance:
            non_full_length_reads_map.append(read.query_name)
            continue

        # check for soft clipping
        cigar = read.cigartuples
        if not cigar:
            soft_clip_5, soft_clip_3 = 0, 0
        else:
            soft_clip_5 = (
                cigar[0][1] if cigar[0][0] == 4 else 0
            )  # 4 indicates soft clipping
            soft_clip_3 = cigar[-1][1] if cigar[-1][0] == 4 else 0

        if soft_clip_5 <= tolerance and soft_clip_3 <= tolerance:
            full_length_reads.append(read.query_name)
        else:
            non_full_length_reads_clip.append(read.query_name)

    # Output results as dataframe
    read_data = pd.DataFrame(
        {
            "chrom": chrom,
            "start": start,
            "end": end,
            "full_length_reads": [full_length_reads],
            "non_full_length_reads": [
                non_full_length_reads_map + non_full_length_reads_clip
            ],
            "non_primary_reads": [non_primary_reads],
        }
    )

    return read_data


def count_total_fibers(bam):
    """
    Count the total number of fibers in a BAM file.
    This counts primary reads and unmapped reads, representing each fiber once.
    """
    bam.reset()

    def is_primary(read):
        return not (read.is_secondary or read.is_supplementary)

    fiber_count = bam.count(
        read_callback=is_primary, until_eof=True
    )  # Counts primary and unmapped reads

    return fiber_count


def aggregate_strand_metrics(table, total_fibers):
    """
    Counts the number of full-length, non-full-length, and non-primary reads and adds them to a summary table.
    """

    aggregate_table = pd.DataFrame(
        {
            "chrom": table["chrom"],
            "start": table["start"],
            "end": table["end"],
            "#_full_length_reads": table["full_length_reads"].apply(len),
            "#_non_full_length_reads": table["non_full_length_reads"].apply(len),
            "#_non_primary_reads": table["non_primary_reads"].apply(len),
            "total_fibers in bam(primary+unmapped)": total_fibers,
        }
    )

    return aggregate_table


# Script excecution

bam_path = snakemake.input.data
regions = snakemake.params.regions
end_tolerance = snakemake.params.end_tolerance
output_detailed = snakemake.output.detailed
output_summary = snakemake.output.summary
threads = snakemake.threads

# for testing and troubleshooting
# bam_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test_2/align/htt_test_2.mapped.reads.bam"
# regions=["chr4:3073138-3075853", "chr3:179228176-179236561"]
# end_tolerance=30
# output_detailed="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/detailed-tar.tsv.gz"
# output_summary="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/summary-tar.tsv.gz"

os.makedirs(os.path.dirname(output_detailed), exist_ok=True)

pybam = pysam.AlignmentFile(bam_path, "rb", threads=threads)

tables = []

for region in regions:
    chrom, start, end = parse_region(region)
    reg_table = region_metrics_table(pybam, chrom, start, end, tolerance=end_tolerance)

    tables.append(reg_table)
table = pd.concat(tables, ignore_index=True)

total_fibers = count_total_fibers(pybam)

aggregate_table = aggregate_strand_metrics(table, total_fibers)

# save tables
for column in ["full_length_reads", "non_full_length_reads", "non_primary_reads"]:
    table[column] = table[column].apply(
        lambda x: ",".join(x) if isinstance(x, list) else ""
    )

table.to_csv(output_detailed, sep="\t", index=False, compression="gzip")
aggregate_table.to_csv(output_summary, sep="\t", index=False)