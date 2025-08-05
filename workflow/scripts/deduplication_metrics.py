import pysam
import pandas as pd
from collections import Counter


def parse_region(region):
    # region should be in chr:start-end format
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    return chrom, start, end


def count_du_tags(pybam, chrom, start, end):
    du_counts = Counter()
    for read in pybam.fetch(chrom, start, end):
        du_tag = read.get_tag("du") if read.has_tag("du") else read.query_name
        du_counts[du_tag] += 1

    sorted_du_counts = du_counts.most_common()
    if not sorted_du_counts:
        return None, None

    du_tags = ",".join([tag for tag, count in sorted_du_counts])
    counts = ",".join([str(count) for tag, count in sorted_du_counts])

    return du_tags, counts


# Main script execution
regions = snakemake.params.regions
bamfile = snakemake.input.bam
output = snakemake.output.deduplication_metrics
threads = snakemake.threads

# testing
# regions="chr4:3073138-3075853"
# bamfile="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/temp/htt_test/align/htt_test.filtered.bam"
# output="/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/htt_test_deduplication.tsv.gz"

bam = pysam.AlignmentFile(bamfile, "rb", threads=threads)

region_collector = []

for region in regions:
    chrom, start, end = parse_region(region)

    du_tags, counts = count_du_tags(bam, chrom, start, end)

    region_collector.append(
        {
            "chrom": chrom,
            "start": start,
            "end": end,
            "du_tags": du_tags,
            "counts": counts,
        }
    )

# create a dataframe from the region_collector
du_dataframe = pd.DataFrame(region_collector)

# save the dataframe to a csv file
du_dataframe.to_csv(output, sep="\t", index=False, compression="gzip")
# close the bam file
bam.close()
