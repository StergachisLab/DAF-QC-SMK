import pysam
import pyabpoa as pa
from collections import defaultdict


def extract_du_seqs(pysam_file):
    reads = defaultdict(list)
    readgroups = defaultdict(str)
    for read in pysam_file.fetch():
        seq = read.query_sequence
        rg = read.get_tag("RG") if read.has_tag("RG") else ""
        if read.has_tag("du"):
            du = read.get_tag("du")
            reads[du].append(seq)
            if rg and rg not in readgroups[du]:
                readgroups[du] = rg
        else:
            # Default reads, those that give the du tag name, are not marked with the du tag, so this is necessary
            du = read.query_name
            reads[du].insert(0, seq)
            if rg and rg not in readgroups[du]:
                readgroups[du] = rg

    return reads, readgroups


def abpoa_MSA(seq_list):
    # abPOA implementation.
    """Takes list of sequences as input and returns abPOA consensus and MSA.
    Uses abPOA built-in consensus generator"""
    aligner = pa.msa_aligner()
    alignment = aligner.msa(seq_list, out_cons=True, out_msa=True)

    consensus = alignment.cons_seq[0]
    msa_alignment = alignment.msa_seq

    return consensus, msa_alignment


def collect_consensus(dup_dict, min_read_count=3):

    du_names = dup_dict.keys()
    dups_consensus = []

    for du in du_names:
        dups_read_count = len(dup_dict[du])
        
        if dups_read_count < min_read_count:
            continue

        if dups_read_count == 1:
            # If only one read, use it as consensus
            consensus_seq = dup_dict[du][0]
            dups_consensus.append([du, consensus_seq, dups_read_count])
            continue

        elif dups_read_count == 2:
            # Use the pbmarkdup default read, which is the first in the list
            consensus_seq = dup_dict[du][0]
            dups_consensus.append([du, consensus_seq, dups_read_count])
            continue

        sequences = dup_dict[du]

        consensus_seq, rep_MSA = abpoa_MSA(sequences)
        dups_consensus.append([du, consensus_seq, dups_read_count])

    return dups_consensus


def consensus_dfm_to_bam(dups_consensus, output_bam, rg_dict):
    for du, consensus_seq, dups_read_count in dups_consensus:
        read_name = f"{du}_consensus"
        zm_tag = int(du.split("/")[1])
        rg = rg_dict[du] if du in rg_dict and len(rg_dict[du]) > 0 else None

        a = pysam.AlignedSegment()
        a.query_name = read_name
        a.query_sequence = consensus_seq
        a.flag = 4  # unmapped flag
        a.is_unmapped = True
        a.set_tag("du", du)
        a.set_tag("zm", zm_tag)
        a.set_tag("YC", "240,187,201") # Colors reads, useful for visualization
        a.set_tag("dc", dups_read_count) # Duplication count
        if rg:
            a.set_tag("RG", rg)

        output_bam.write(a)


def dedup_bam(pysam_file, output_bam, min_read_count=3):
    # First create a dfm from the input table to group reads by pbmarkdup du group tags
    dup_dict, rg_dict = extract_du_seqs(pysam_file)
    # Use MSA to get consensus sequences (using default read numbers here)
    dups_consensus = collect_consensus(dup_dict, min_read_count=min_read_count)
    # Create an unaligned bam file from consensus sequences
    consensus_dfm_to_bam(dups_consensus, output_bam, rg_dict)


# Main script execution
bam_file = snakemake.input.bam
output_bam = snakemake.output.bam
min_read_count = snakemake.params.consensus_min_reads
threads = snakemake.threads

# For testing
# bam_file = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/temp/htt_test/align/htt_test.filtered.bam"
# output_bam = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/htt_test.consensus.bam"
# min_read_count = 3
# threads = 8

input_bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)
output_bam = pysam.AlignmentFile(output_bam, "wb", template=input_bam, threads=threads)

dedup_bam(input_bam, output_bam, min_read_count=min_read_count)

input_bam.close()
output_bam.close()