import pysam
import pandas as pd
import os


bam_path = snakemake.input.data
regions = snakemake.params.regions
end_tolerance = snakemake.params.end_tolerance
output_detailed = snakemake.output.detailed
output_summary = snakemake.output.summary

# for testing and troubleshooting
#bam_path="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test_2/align/htt_test_2.mapped.reads.bam"
#regions=["chr4:3073138-3075853", "chr3:179228176-179236561"]
#output_detailed="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/detailed-tar.tsv.gz"
#output_summary="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/summary-tar.tsv.gz"

def parse_region(region):
    # region should be in chr:start-end format
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    return chrom, start, end


def region_metrics_table(bam, chrom, start, end, tolerance=30):
    # Filters for reads that align at the expected ends (+/- tolerace)
    # and requires that reads do not have more than 30 bp of soft clipping.
    # Filters out non-primary alignments

    full_length_reads = []
    non_full_length_reads_map = []
    non_full_length_reads_clip = []
    non_primary_reads = []
    read_count = 0
#    map_start, map_end = None, None

    for read in bam.fetch(chrom, start, end):
        read_count += 1
        
        map_start, map_end = read.reference_start, read.reference_end
        # skipping reads that are unmapped, secondary, or supplementary
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            non_primary_reads.append(read.query_name)
            continue

        # check for full length alignments



        map_start_diff = map_start - start
        map_end_diff = map_end - end

        if abs(map_start_diff) > tolerance or abs(map_end_diff) > tolerance:
            non_full_length_reads_map.append(read.query_name)
            continue

        # check for soft clipping
        cigar = read.cigartuples
        if cigar[0][0] == 4:  # 4 indicates soft clipping
            soft_clip_5 = cigar[0][1]
        else:
            soft_clip_5 = 0
        if cigar[-1][0] == 4:  # 4 indicates soft clipping
            soft_clip_3 = cigar[-1][1]
        else:
            soft_clip_3 = 0

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
    # This returns the total count of fibers in the BAM file by counting primary reads and unmapped reads, i.e. representing each fiber once
    def is_primary(read):
        return not (read.is_secondary or read.is_supplementary)

    primary_count = bam.count(read_callback=is_primary)
    unmapped_count = bam.count(until_eof=True)
    total_count = primary_count + unmapped_count
    return total_count


def aggregate_strand_metrics(table, total_fibers):
    # Aggregate metrics across all regions
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


os.makedirs(os.path.dirname(output_detailed), exist_ok=True)

# bam = pyft.Fiberbam(bam_path)
pybam = pysam.AlignmentFile(bam_path, "rb")


for region in regions:
    chrom, start, end = parse_region(region)


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
