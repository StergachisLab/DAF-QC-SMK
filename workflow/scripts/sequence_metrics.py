import pysam
import numpy as np
import pandas as pd
import os
#import time

# Inputs: aligned bam file with MD tags, regions (string), output file location
# outputs: csv or text files of data for each region, list of CT and GA reads, possibly: bam file with ambiguity codes

bam_path = snakemake.input.data
regions = snakemake.params.regions
targeting_metrics = snakemake.input.targeting_data
chimera_cutoff = snakemake.params.chimera_cutoff
min_deamination_count = snakemake.params.min_deamination_count
read_metrics = snakemake.output.read_metrics
summary_metrics = snakemake.output.summary_metrics

# for testing
#bam_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/align/htt_test.mapped.reads.bam"
#regions = ["chr4:3073138-3075853", "chr3:179228176-179236561"]
#targeting_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test.detailed_targeting_metrics.tbl.gz"
#chimera_cutoff = 0.9
#min_deamination_count = 50 # minimum number of deaminations to designate a strand (otherwise undetermined)
#read_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.read_metrics_optimized.tsv.gz"
#summary_metrics = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.summary_metrics_optimized.tsv.gz"



def parse_region(region):
    # region should be in chr:start-end format
    chrom, positions = region.split(':')
    start, end = map(int, positions.split('-'))
    return chrom, start, end



def strand_metrics(read, cutoff=0.9):
        # based on the proportion of C->T & G->A determine the strand acted upon by DddA
    # only counting single base substitutions

    c = 0
    g = 0
    non_ref = 0

#    start_strand = time.time()

    seq = read.query_sequence
    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)


    for pos in pair:
        if pos[0] is None or pos[1] is None: # indel, ignore
            continue
        else:
            strand_base= seq[pos[0]]
#            qi = pos[0]
            ref_base = pos[2].upper()
            if strand_base != ref_base:
                non_ref += 1
                if ref_base == 'C' and strand_base == 'T':
                    c += 1
                elif ref_base == 'G' and strand_base == 'A':
                    g += 1

    
    if c+g == 0:
        strand='none'
    elif c+g/non_ref <= cutoff or c+g < min_deamination_count:
        strand = 'undetermined'
    elif c/(c+g) >= cutoff:
        strand = 'CT'
    elif g/(c+g) >= cutoff:
        strand = 'GA'
    else:
        strand = 'chimera'
    
#    print("designate_strand", time.time()-start_strand)

    if strand in ["CT", "GA"]:
        deamination_pos = [] # deaminated positions in read coordinates
        doublets = {'AC':[0,0], 'CC': [0,0], 'GC':[0,0], 'TC':[0,0], 'OC':[0,0]} # OC represents no base/indel before C
#        doublet_dict = {'GA': 'TC', 'GC': 'GC', 'GG': 'CC', 'GT': 'AC', 'GO': 'OC',
#                    'AC': 'AC', 'CC': 'CC', 'GC':'GC', 'TC':'TC', 'OC':'OC'}
        doublet_dict = {'GA': 'TC', 'GC': 'GC', 'GG': 'CC', 'GT': 'AC'}
        
        if strand == "CT":
            mutation_count = non_ref - c # 
        else:
            mutation_count = non_ref - g

        for i, pos in enumerate(pair):
            if pos[0] == None or pos[1] == None: # indel, ignore
                continue

            ref_base = pos[2].upper()

            if ref_base != strand[0]: # ignore non C/G bases
                continue

            strand_base = seq[pos[0]]

            if strand == 'CT':
                doublet = pair[i-1][2].upper() + ref_base if pair[i-1][2] is not None and i>0 else 'OC'
            elif strand == 'GA':
#                doublet = ref_base + pair[i+1][2].upper() if i + 1 < len(pair) and pair[i+1][2] is not None else 'GO'
                doublet = doublet_dict[ref_base + pair[i+1][2].upper()] if i + 1 < len(pair) and pair[i+1][2] is not None else 'OC'


    # Record total counts for each pair type
#           doublets[doublet_dict[doublet]][0] += 1
            doublets[doublet][0] += 1



            if strand_base != ref_base and strand_base == strand[1]:
                deamination_pos.append(pos[0])
                doublets[doublet][1] += 1

    else:
        doublets = {'AC':[None,None], 'CC': [None,None], 'GC':[None,None], 'TC':[None,None], 'OC':[None,None]} # OC represents no base/indel before C
        mutation_count = None
        deamination_pos = None

    return strand, doublets, mutation_count, deamination_pos



def strand_metrics_table (psfile, chrom, start, end, chimera_cutoff=0.9, include_readnames=None):
    # psfile should be pysam alignment file
    # include_readnames is a list of read names to include, if None all reads are included
    
    read_collector=[]
    
	# If a read filter is applied, 

    # check that this region has reads
    if psfile.count(chrom, start, end) == 0:
        return pd.DataFrame()

    for read in psfile.fetch(chrom, start, end):
        if include_readnames is not None and read.query_name not in include_readnames:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        


        strand, doublets, mutation_count, deam_pos = strand_metrics(read, chimera_cutoff)

        duplicate = read.get_tag('du') if read.has_tag('du') else "None"

        read_data = {           
            'read_name': read.query_name, 
            'chr': chrom, 
            'reg_st': start, 
            'reg_end': end,
            'strand_st': read.reference_start,
            'strand_end': read.reference_end,
            'length': len(read.query_sequence),
            'strand': strand,
            'duplicate' : duplicate,
            'mutation_count' : mutation_count,
            'deamination_positions' : ','.join(map(str, deam_pos)) if deam_pos is not None else ''
        }
        
        keys = ['AC', 'CC', 'GC', 'TC', 'OC']

        if strand in ['CT', 'GA']:
            for key in keys:
                read_data[f'{key}_count'] = doublets[key][0]
                read_data[f'{key}_deam'] = doublets[key][1]
                read_data[f'{key}_deam_rate'] = doublets[key][1]/doublets[key][0] if doublets[key][0] > 0 else None


        read_collector.append(read_data)

    reads_table = pd.DataFrame(read_collector)

    reads_table['total_deam'] = reads_table[[f'{key}_deam' for key in keys]].sum(axis=1)
    reads_table['total_count'] = reads_table[[f'{key}_count' for key in keys]].sum(axis=1)
    reads_table['all_deam_rate'] = (
        reads_table['total_deam'] / reads_table['total_count']
    ).where(reads_table['total_count'] > 0)
    reads_table['mutation_rate'] = (
        reads_table['mutation_count'] / reads_table['length']
    ).where(reads_table['length'] > 0)


    return reads_table


def aggregate_strand_metrics(table):
#    aggregate_template = {'chrom': chrom, 'start': start, 'end': end, 'strand':None, 'count':None, 'mutation_rate':[], 'all_deam_rate':[],
#                          'AC_deam_rate':[], 'CC_deam_rate':[], 'GC_deam_rate':[], 'TC_deam_rate':[], 'OC_deam_rate':[]}
    aggregate_collector = []

    for group in table.groupby(['chr', 'reg_st', 'reg_end', 'strand']):
        mutation_rate = group[1]['mutation_rate'].dropna().tolist()
        all_deam_rate = group[1]['all_deam_rate'].dropna().tolist()
        AC_deam_rates = group[1]['AC_deam_rate'].dropna().tolist()
        CC_deam_rates = group[1]['CC_deam_rate'].dropna().tolist()
        GC_deam_rates = group[1]['GC_deam_rate'].dropna().tolist()
        TC_deam_rates = group[1]['TC_deam_rate'].dropna().tolist()
        OC_deam_rates = group[1]['OC_deam_rate'].dropna().tolist()
        count = len(group[1])

        aggregate_template = {
            'chrom': group[0][0], 
            'reg_start': group[0][1], 
            'reg_end': group[0][2], 
            'strand': group[0][3], 
            'count': count, 
            'mutation_rate': ','.join(map(str, mutation_rate)) if mutation_rate else '', 
            'all_deam_rate': ','.join(map(str, all_deam_rate)) if all_deam_rate else '',
            'AC_deam_rate': ','.join(map(str, AC_deam_rates)) if AC_deam_rates else '', 
            'CC_deam_rate': ','.join(map(str, CC_deam_rates)) if CC_deam_rates else '',
            'GC_deam_rate': ','.join(map(str, GC_deam_rates)) if GC_deam_rates else '', 
            'TC_deam_rate': ','.join(map(str, TC_deam_rates)) if TC_deam_rates else '', 
            'OC_deam_rate': ','.join(map(str, OC_deam_rates)) if OC_deam_rates else ''
        }

        aggregate_collector.append(aggregate_template)

        aggregate_table = pd.DataFrame(aggregate_collector)




    return aggregate_table


#start_time = time.time()

os.makedirs(os.path.dirname(read_metrics), exist_ok=True)


pybam=pysam.AlignmentFile(bam_path, 'rb')





tables=[]

if len(targeting_metrics) > 0:

    targeting_df = pd.read_csv(targeting_metrics, sep='\t', compression='gzip')
    for region in regions:
        
        chrom, start, end = parse_region(region)

        full_length_reads = targeting_df[(targeting_df['chrom'] == chrom) & (targeting_df['start'] == start) & (targeting_df['end'] == end)]['full_length_reads'].values

        if full_length_reads is None or full_length_reads[0] is np.nan or len(full_length_reads) == 0:
            continue

        full_length_reads = full_length_reads[0].split(',')
    #    print (f"Region {region} has {len(full_length_reads)} full length reads")
        reg_table=strand_metrics_table(pybam, chrom, start, end, chimera_cutoff=chimera_cutoff, include_readnames=full_length_reads,)

        tables.append(reg_table)



else:    
    for region in regions:
        chrom, start, end = parse_region(region)
        reg_table=strand_metrics_table(pybam, chrom, start, end, chimera_cutoff=chimera_cutoff)
        tables.append(reg_table)
table=pd.concat(tables, ignore_index=True)

#table_time = time.time()
#print("tables", table_time-start_time)


aggregate_table=aggregate_strand_metrics(table)



table.to_csv(read_metrics, sep='\t', index=False, compression='gzip')
aggregate_table.to_csv(summary_metrics, sep='\t', index=False, compression='gzip')

#end_time = time.time()
#elapsed_time = end_time - start_time
#print(f"Elapsed time: {elapsed_time:.2f} seconds")