def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


# Basic rules needed
# Align bam, if align is true
# Targeting metrics- per target % of reads, all targets % of reads, metrics for full length reads
# Mutation rate metrics, output C,G, or other strand designation
# Deduplication
# Targeting metrics for deduplicated files
# Mutation rate metrics for deduplicated files

# Later, consensus generation. For now, this should not be included by default
# Consensus sequence bam generation- filter for just full-length and C>T/G>A designated strands. Use MSA to generate a consensus. Map consensus. 
