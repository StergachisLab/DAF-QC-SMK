import pysam
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Inputs: aligned bam file with MD tags, regions (string), output file location
# outputs: csv or text files of data for each region, list of CT and GA reads, possibly: bam file with ambiguity codes

bam_path = snakemake.input.data
regions = snakemake.params.regions





#label = "100k"

#bam_path= "/mmfs1/gscratch/stergachislab/bohaczuk/analysis/chdi-hd/primer-design/25.3.7_vega/HTTLNA_hg38_with_md.bam"
#bam_path=f"/gscratch/stergachislab/bohaczuk/data/DAF_processing/chris_PIK3CA_plasmidsaurus/GM12878_{label}_DddA_with_MD.bam"
#region="chr4:3073071-3076052"
#chrom='chr3'
#start= 179228176
#end= 179236561
#save = False
#ref = "/mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa"
#chrom = 'chr4'
#start = 3073071
#end = 3076052

def parse_region(region):
	# region should be in chr:start-end format
	chrom, positions = region.split(':')
	start, end = map(int, positions.split('-'))
	return chrom, start, end


def determine_da_strand(read, cutoff=0.9): # Modified from Elliott's script
    # based on the proportion of C->T & G->A determine the strand acted upon by DddA
    # only counting single base substitutions

    c = 0
    g = 0
    total = 0

    seq = read.query_sequence
    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)

    for pos in pair:
        if pos[0] == None or pos[1] == None: # indel, ignore
            pass
        else:
            strand_base= seq[pos[0]]
#            qi = pos[0]
            ref_base = pos[2].upper()
            if strand_base != ref_base:
                total += 1
                change = ref_base + strand_base
                if change == "CT":
                    c += 1
                elif change == 'GA':
                    g += 1

    if c+g == 0:
        return('none')
    elif c+g/total <= cutoff:
        return('undetermined')
    elif c/(c+g) >= cutoff:
        return('CT')
    elif g/(c+g) >= cutoff:
        return('GA')
    else:
        return('chimera')
    

def strand_metrics(read, strand):
    # Requires CT or GA strand designation to calculate metrics
    read_name = read.query_name
    deamination_pos = [] # deaminated positions in read coordinates
    mutation_count = 0 # Non C>T (top) or G>A (bottom) mutation count
#    ref_S_count = 0 # total C (top) or G (bottom) reference positions
    # TODO, delete ref_S_count if doublet approach works
    doublets = {'AC':[0,0], 'CC': [0,0], 'GC':[0,0], 'TC':[0,0], 'OC':[0,0]} # OC represents no base/indel before C
#    GA_to_CT_dict = {'GA': 'TC', 'GC': 'GC', 'GG': 'CC', 'GT': 'AC', 'GO': 'OC'}
    doublet_dict = {'GA': 'TC', 'GC': 'GC', 'GG': 'CC', 'GT': 'AC', 'GO': 'OC',
                    'AC': 'AC', 'CC': 'CC', 'GC':'GC', 'TC':'TC', 'OC':'OC'}

    seq = read.query_sequence
    seq_len = len(seq)
    pair=read.get_aligned_pairs(matches_only=False, with_seq=True)
        

    for i, pos in enumerate(pair):



        if pos[0] == None or pos[1] == None: # indel, ignore
            continue
        
        strand_base = seq[pos[0]]
        ref_base = pos[2].upper()

        if ref_base != strand[0]: # ignore non C/G bases
            if strand_base == ref_base:
                continue
            else:
# Count mutations at non C (for CT) or G (for GA) positions
                mutation_count += 1
                continue



        if strand == 'CT':
#            doublet = seq[pos[0]-1] + seq[pos[0]] if pos[0] > 0 else 'OC'
            doublet = pair[i-1][2].upper() + ref_base if pair[i-1][2] is not None and i>0 else 'OC'
        elif strand == 'GA':
#            doublet = seq[pos[0]] + seq[pos[0]+1] if pos[0] < seq_len - 1 else 'GO'
            doublet = ref_base + pair[i+1][2].upper() if i + 1 < len(pair) and pair[i+1][2] is not None else 'GO'

# Record total counts for each pair type
        doublets[doublet_dict[doublet]][0] += 1        

        if strand_base != ref_base:
            change = ref_base + strand_base
            if change == strand:
# Record deamination event at each pair type
                deamination_pos.append(pos[0])
                doublets[doublet_dict[doublet]][1] += 1
            else:
# Count mutations at C (for CT) or G (for GA) positions that are not deaminations
                mutation_count += 1


#        ref_S_count += 1
                    
#    deamination_rate = sum([doublets[key][1] for key in doublets])/sum([doublets[key][0] for key in doublets])
#    print('new',deamination_rate)


#    deamination_rate_old = len(deamination_pos)/ref_S_count
#    print('old',deamination_rate_old)
#    mutation_rate = mutation_count/seq_len

#    return deamination_rate, mutation_rate, seq_len, deamination_pos, doublets


    return doublets, mutation_count
# output read name, strand designation, AC_count, AC_deam, CC_count, CC_deam, GC_count, GC_deam, TC_count, TC_deam, OC_count, OC_deam, mutation_count, read_len


def strand_metrics_table (psfile, chrom, start, end):
    # psfile should be pysam alignment file
    
    read_collector=[]

    # check that this region has reads
    if psfile.count(chrom, start, end) == 0:
        return pd.DataFrame()

    for read in psfile.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary:
            continue

        strand = determine_da_strand(read)

        read_data = {           
            'read_name': read.query_name, 
            'chr': chrom, 
            'st': start, 
            'end': end, 
            'length': len(read.query_sequence),
            'strand': strand,
        }
        


        if strand in ['CT', 'GA']:
            doublets, mutation_count = strand_metrics(read, strand)
            read_data['mutation_count'] = mutation_count
            for key in ['AC', 'CC', 'GC', 'TC', 'OC']:
                read_data[f'{key}_count'] = doublets[key][0]
                read_data[f'{key}_deam'] = doublets[key][1]
                read_data[f'{key}_deam_rate'] = doublets[key][1]/doublets[key][0] if doublets[key][0] > 0 else None
        else:
            read_data['mutation_count'] = None
            for key in ['AC', 'CC', 'GC', 'TC', 'OC']:
                read_data[f'{key}_count'] = None
                read_data[f'{key}_deam'] = None
                read_data[f'{key}_deam_rate'] = None


        read_collector.append(read_data)

    reads_table = pd.DataFrame(read_collector)

    keys = ['AC', 'CC', 'GC', 'TC', 'OC']
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

    for group in table.groupby(['chr', 'st', 'end', 'strand']):
        mutation_rate = group[1]['mutation_rate'].dropna().tolist()
        all_deam_rate = group[1]['all_deam_rate'].dropna().tolist()
        AC_deam_rates = group[1]['AC_deam_rate'].dropna().tolist()
        CC_deam_rates = group[1]['CC_deam_rate'].dropna().tolist()
        GC_deam_rates = group[1]['GC_deam_rate'].dropna().tolist()
        TC_deam_rates = group[1]['TC_deam_rate'].dropna().tolist()
        OC_deam_rates = group[1]['OC_deam_rate'].dropna().tolist()
        count = len(group[1])

        
        aggregate_template = {'chrom': group[0][0], 'start': group[0][1], 'end': group[0][2], 
                            'strand': group[0][3], 'count':count, 'mutation_rate': mutation_rate, 'all_deam_rate': all_deam_rate,
                            'AC_deam_rate': AC_deam_rates, 'CC_deam_rate': CC_deam_rates,
                            'GC_deam_rate': GC_deam_rates, 'TC_deam_rate': TC_deam_rates, 'OC_deam_rate': OC_deam_rates}
        aggregate_collector.append(aggregate_template)

        aggregate_table = pd.DataFrame(aggregate_collector)




    return aggregate_table




pybam=pysam.AlignmentFile(bam_path, 'rb')

tables=[]

for region in regions:
    chrom, start, end = parse_region(region)
    reg_table=strand_metrics_table(pybam, chrom, start, end)

    tables.append(reg_table)
table=pd.concat(tables, ignore_index=True)


aggregate_table=aggregate_strand_metrics(table)


# save tables
for column in ['mutation_rate', 'all_deam_rate', 'AC_deam_rate', 'CC_deam_rate', 'GC_deam_rate', 'TC_deam_rate', 'OC_deam_rate']:
    table[column] = table[column].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else '')
    aggregate_table[column] = aggregate_table[column].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else '')

table.to_csv(snakemake.output.read_metrics, sep='\t', index=False, compression='gzip')
aggregate_table.to_csv(snakemake.output.summary_metrics, sep='\t', index=False, compression='gzip')

