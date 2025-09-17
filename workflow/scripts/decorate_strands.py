import pysam
import pandas as pd

bam_name = snakemake.input.bam
seq_metrics = snakemake.input.seq_metrics
output_bam = snakemake.output.decorated_bam
threads = snakemake.threads

# bam_name = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/temp/htt_test/align/htt_test.filtered.bam"
# seq_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/htt_test.detailed_seq_metrics.reads.tbl.gz"
# output_bam = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/align/htt_test.decorated.bam"

# For each read, check in table. Get strand, replace with ambiguity code at designated location

# write corrected reads to new BAM
bam = pysam.AlignmentFile(bam_name, "rb", threads=threads)
table = pd.read_csv(seq_metrics, sep="\t")
corrected_bam = pysam.AlignmentFile(output_bam, "wb", template=bam, threads=threads)
for read in bam.fetch():
    if read.is_secondary is True or read.is_supplementary is True:
        continue

    read_id = read.query_name
    # Determine strand and deamination positions
    deam_dictionary = {"CT": "Y", "GA": "R"}
    deam_pos = (
        table.loc[table["read_name"] == read_id, "deamination_positions"]
        .values[0]
        .strip("[]")
        .split(",")
    )
    deam_pos = [int(x) for x in deam_pos if x != ""]
    strand = table.loc[table["read_name"] == read_id, "strand"].values[0]
    read_seq = list(read.query_sequence)
    for pos in deam_pos:
        assert read_seq[pos] == strand[1], (
            f"Position {pos} in read {read_id} does not match expected base {strand[0]}, found {read_seq[pos]}"
        )
        read_seq[pos] = deam_dictionary[strand]
    new_seq = "".join(read_seq)
#    deam_pos = [pos + 1 for pos in deam_pos]  # convert to 1-based positions
    MD = read.get_tag("MD")

    quals = read.query_qualities
    read.query_sequence = new_seq
    read.query_qualities = quals

#    read.set_tags(
#        [
##            ("da", deam_pos),
##            ("fd", deam_pos[0] if deam_pos else 0, "i"),
##            ("ld", deam_pos[-1] if deam_pos else 0, "i"),
#            ("st", strand),
#            ("MD", MD),
#        ]
#    )
    read.set_tag("st", strand)
    read.set_tag("MD", MD)
    corrected_bam.write(read)
corrected_bam.close()
bam.close()
