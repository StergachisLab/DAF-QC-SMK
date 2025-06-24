#import pyft
import pysam
import numpy as np
import pandas as pd

bam_path = snakemake.input.data
regions = snakemake.params.regions

#label = "1M"

#chrom = 'chr3'
#bam_path=f"/gscratch/stergachislab/bohaczuk/data/DAF_processing/chris_PIK3CA_plasmidsaurus/GM12878_{label}_DddA_with_MD.bam"
#regions=["chr4:3073071-3076052"]
#regions=["chr3:179229957-179234791", "chr3:179229957-179234791"]
#chrom='chr3'
#start= 179229957
#end= 179234791

def parse_region(region):
	# region should be in chr:start-end format
	chrom, positions = region.split(':')
	start, end = map(int, positions.split('-'))
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

    for read in bam.fetch(chrom, start, end):
        read_count +=1
# skipping reads that are unmapped, secondary, or supplementary        
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            non_primary_reads.append(read.query_name)
            continue

# check for full length alignments

        map_start, map_end = read.reference_start, read.reference_end

        map_start_diff =  map_start - start
        map_end_diff = map_end - end
#        print("map start", map_start, "map end", map_end, "start", start, "end", end, "diffs", map_start_diff, map_end_diff)

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

#        print("clipping", soft_clip_5, soft_clip_3)

        if soft_clip_5 <= tolerance and soft_clip_3 <= tolerance:
            full_length_reads.append(read.query_name)
#            print('full')
        else:
            non_full_length_reads_clip.append(read.query_name)

# Output results as dataframe
    read_data = pd.DataFrame({
        'chrom': chrom,
        'start': start,
        'end': end,
        'full_length_reads': [full_length_reads],
        'non_full_length_reads_map': [non_full_length_reads_map],
        'non_full_length_reads_clip': [non_full_length_reads_clip],
        'non_primary_reads': [non_primary_reads],
        'read_count': read_count})
    
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
    aggregate_table = table.copy()
    # Get proportions of reads in each category
    chroms = table['chrom']
    starts = table['start']
    ends = table['end']
    full_length_reads = table['full_length_reads'].apply(len)
    partial_reads = table['non_full_length_reads_map'].apply(len) + table['non_full_length_reads_clip'].apply(len)
    non_primary_reads = table['non_primary_reads'].apply(len)


    aggregate_table = pd.DataFrame({
        'chrom': chroms,
        'start': starts,
        'end': ends,
        '#_full_length_reads': full_length_reads,
        '#_partial_reads': partial_reads,
        '#_non_primary_reads': non_primary_reads,
        'total_fibers in bam(primary+unmapped)': total_fibers
    })

    return aggregate_table


#bam = pyft.Fiberbam(bam_path)
pybam = pysam.AlignmentFile(bam_path, 'rb')


for region in regions:
    chrom, start, end = parse_region(region)
     

tables=[]

for region in regions:
    chrom, start, end = parse_region(region)
    reg_table=region_metrics_table(pybam, chrom, start, end)

    tables.append(reg_table)
table=pd.concat(tables, ignore_index=True)

total_fibers = count_total_fibers(pybam)

aggregate_table = aggregate_strand_metrics(table, total_fibers)

# save tables
for column in ['full_length_reads', 'non_full_length_reads_map', 'non_full_length_reads_clip', 'non_primary_reads']:
    table[column] = table[column].apply(lambda x: ','.join(x) if isinstance(x, list) else '')


table.to_csv(snakemake.output.detailed, sep='\t', index=False, compression='gzip')
aggregate_table.to_csv(snakemake.output.summary, sep='\t', index=False)


'''

full_length_reads, non_full_length_reads_map, non_full_length_reads_clip, non_primary_reads, read_count = get_full_length_reads(bam, chrom, start, end)

print(f"Total reads processed in region: {read_count}")
print(f"Full length reads: {len(full_length_reads)}")
print(f"Non-full length reads (map): {len(non_full_length_reads_map)}")
print(f"Non-full length reads (clip): {len(non_full_length_reads_clip)}")
print(f"Non-primary reads: {len(non_primary_reads)}")
print(f"Total fibers in region: {len(full_length_reads) + len(non_full_length_reads_map) + len(non_full_length_reads_clip)}")
print(f"Total fibers in bam (inc. unmapped): {count_total_fibers(bam)}")
print(f"Targeting efficiency (full length): {len(full_length_reads) / count_total_fibers(bam) * 100:.2f}%")
print(f"Targeting efficiency (all) : {(len(full_length_reads) + len(non_full_length_reads_map) + len(non_full_length_reads_clip)) / count_total_fibers(bam) * 100:.2f}%")
print(f"Percent of target mapped reads that are full-length: {len(full_length_reads) / (len(full_length_reads) + len(non_full_length_reads_map)+ len(non_full_length_reads_clip)) * 100:.2f}%")







# save tables
for column in ['mutation_rate', 'all_deam_rate', 'AC_deam_rate', 'CC_deam_rate', 'GC_deam_rate', 'TC_deam_rate', 'OC_deam_rate']:
    table[column] = table[column].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else '')
    aggregate_table[column] = aggregate_table[column].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else '')

table.to_csv(snakemake.output.read_metrics, sep='\t', index=False, compression='gzip')
aggregate_table.to_csv(snakemake.output.summary_metrics, sep='\t', index=False, compression='gzip')


'''