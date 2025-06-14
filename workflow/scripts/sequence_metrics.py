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

for region in regions.split(','):
    chrom, start, end = parse_region(region)
    reg_table=strand_metrics_table(pybam, chrom, start, end)

    tables.append(reg_table)
table=pd.concat(tables, ignore_index=True)


aggregate_table=aggregate_strand_metrics(table)


# save tables

table.to_csv(snakemake.output.read_metrics, sep='\t', index=False, compression='gzip')
aggregate_table.to_csv(snakemake.output.summary_metrics, sep='\t', index=False, compression='gzip')


'''

pybam=pysam.AlignmentFile(bam_path, 'rb')

regions.split(',')
for region in regions.split(','):
    chrom, start, end = parse_region(region)


#ref_fasta=Fasta(ref)

read_counter, aggregate_metrics= aggregate_strand_metrics(pybam, chrom, start, end)


# plot histogram of deamination rates
deamination_rates = [x[0] for x in aggregate_metrics]
median_deamination = np.median(deamination_rates)
deamination_10 = np.percentile(deamination_rates, 10)
deamination_90 = np.percentile(deamination_rates, 90)


fig= plt.figure(figsize=(10, 6))
weights = [1/len(deamination_rates)] * len(deamination_rates)  # Normalize histogra
plt.hist(deamination_rates, bins=50, color='blue', alpha=0.7, weights=weights)
plt.xlim(0, 1)  # Set x-axis limit to 0-1

# add median and percentiles to the plot
plt.axvline(median_deamination, color='black', linestyle='dashed', linewidth=1, label=f'Median: {median_deamination:.2f}')
plt.axvline(deamination_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {deamination_10:.2f}')
plt.axvline(deamination_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {deamination_90:.2f}')

# add text label next to median and quartile lines
# Get y-axis limits for positioning text
y_min, y_max = plt.ylim()

x_min, x_max = plt.xlim()

# Add text labels
plt.text(median_deamination , y_max + y_max/80, '50%', rotation=90, va='bottom')
plt.text(deamination_10  , y_max + y_max/80, '10%', rotation=90, va='bottom')
plt.text(deamination_90 , y_max + y_max/80, '90%', rotation=90, va='bottom')

plt.text(x_max - x_max/4, y_max * 0.9, f'Median: {median_deamination:.2f}', color='black', fontsize=10)
plt.text(x_max - x_max/4, y_max * 0.8, f'10th Percentile: {deamination_10:.2f}', color='black', fontsize=10)
plt.text(x_max - x_max/4, y_max * 0.7, f'90th Percentile: {deamination_90:.2f}', color='black', fontsize=10)


plt.title('Deamination Rate ' + label)
plt.xlabel('Deamination Rate')
plt.ylabel('Frequency')
if save == True:
    plt.savefig(f'testing/deamination_rate_histogram_{label}.pdf', format='pdf')
plt.show()

# add median of deamination rates to the plot
# 10th and 90th percentile 
# x axis consistent (0-1)

# plot histogram of mutation rates
mutation_rates = [x[1] for x in aggregate_metrics]
median_mutation = np.median(mutation_rates)
mutation_10 = np.percentile(mutation_rates, 10)
mutation_90 = np.percentile(mutation_rates, 90)

upper_limit = 0.02
mutation_rates = [x if x <= upper_limit else upper_limit for x in mutation_rates]
  # Collapse higher rates into last bin
fig= plt.figure(figsize=(10, 6))
weights = [1/len(mutation_rates)] * len(mutation_rates)  # Normalize histogram
plt.hist(mutation_rates, bins=50, color='red', alpha=0.7, weights=weights)
plt.xlim(0, upper_limit)  # Set x-axis limit to 0-0.02

plt.axvline(median_mutation, color='black', linestyle='dashed', linewidth=1, label=f'Median: {median_mutation:.6f}')
plt.axvline(mutation_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {mutation_10:.6f}')
plt.axvline(mutation_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {mutation_90:.6f}')

# add text label next to median and quartile lines
# Get y-axis limits for positioning text
y_min, y_max = plt.ylim()
x_min, x_max = plt.xlim()

# Add text labels
plt.text(median_mutation, y_max + y_max/80, '50%', rotation=90, va='bottom')
plt.text(mutation_10, y_max + y_max/80, '10%', rotation=90, va='bottom')
plt.text(mutation_90, y_max + y_max/80, '90%', rotation=90, va='bottom')

plt.text(x_max - x_max/4, y_max * 0.9, f'Median: {median_mutation:.6f}', color='black', fontsize=10)
plt.text(x_max - x_max/4, y_max * 0.8, f'10th Percentile: {mutation_10:.6f}', color='black', fontsize=10)
plt.text(x_max - x_max/4, y_max * 0.7, f'90th Percentile: {mutation_90:.6f}', color='black', fontsize=10)



plt.title('Non-reference Variant Rate ' + label)
plt.xlabel('Mutation Rate')
plt.ylabel('Frequency')
if save == True:
    plt.savefig(f'testing/mutation_rate_histogram_{label}.pdf', format='pdf')
plt.show()

#mutation rate (0-0.02), set as default. Collapse higher into last bin
# nonref variant rate


# Calculate proportion of CT, GA, chimeric, undetermined, and none reads
total_reads = sum(read_counter.values())
proportions = {key: value / total_reads for key, value in read_counter.items()}
print(f"Read Proportions {label}:")
for key, value in proportions.items():
    print(f"{key}: {value:.4%}, {read_counter[key]} reads")
if save == True:
    with open(f'testing/read_proportions_{label}.txt', 'w') as f:
        f.write(f"Read Proportions {label}:\n")
        for key, value in proportions.items():
            f.write(f"{key}: {value:.4%}, {read_counter[key]} reads\n")



# plot read proportions as a stacked plot
# Define the order with CT and GA at bottom
ordered_keys = ['CT', 'GA', 'chimera', 'none', 'undetermined']
values = [proportions[key] for key in ordered_keys]

# Create figure
fig, ax = plt.subplots(figsize=(6, 4))


# Define colors for each category
colors = {'CT': 'red',
          'GA': 'green', 
          'chimera': '#ff7f0e',  # orange
          'undetermined': '#d62728',  # red
          'none': '#9467bd'  # purple
         }


# Create stacked bar
bottom = 0
for i, (cat, value) in enumerate(zip(ordered_keys, values)):
    ax.bar(0, value, bottom=bottom, label=f'{cat} ({value:.4%})', width=0.75, color=colors[cat])
       
    bottom += value

ax.set_xlim(-1, 1)
ax.set_ylim(0, 1)
ax.set_ylabel('Proportion of Reads')
ax.set_title(f'Read Classification Proportions {label}')
ax.set_xticks([])
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Get handles and labels and reverse them
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
if save == True:
    plt.savefig(f'testing/read_proportions_{label}.pdf', format='pdf')
plt.show()




# plot deamination of doublets
doublet_dict= [x[4] for x in aggregate_metrics]

AC_values = [read['AC'][1]/read['AC'][0] for read in doublet_dict if read['AC'][0] > 0]
CC_values = [read['CC'][1]/read['CC'][0] for read in doublet_dict if read['CC'][0] > 0]
GC_values = [read['GC'][1]/read['GC'][0] for read in doublet_dict if read['GC'][0] > 0]
TC_values = [read['TC'][1]/read['TC'][0] for read in doublet_dict if read['TC'][0] > 0]

#for read in doublet_dict:
#    for key in read:
#        read[key] = read[key][1]/ read[key][0] if read[key][0] > 0 else None

#AC_values = [read['AC'] for read in doublet_dict if read['AC'] is not None]
#CC_values = [read['CC'] for read in doublet_dict if read['CC'] is not None]
#GC_values = [read['GC'] for read in doublet_dict if read['GC'] is not None]
#TC_values = [read['TC'] for read in doublet_dict if read['TC'] is not None]

AC_stats = np.median(AC_values), np.percentile(AC_values, 10), np.percentile(AC_values, 90)
CC_stats = np.median(CC_values), np.percentile(CC_values, 10), np.percentile(CC_values, 90)
GC_stats = np.median(GC_values), np.percentile(GC_values, 10), np.percentile(GC_values, 90)
TC_stats = np.median(TC_values), np.percentile(TC_values, 10), np.percentile(TC_values, 90)
#OC_values = [read['OC'] for read in doublet_dict]

non_TC_values = AC_values + CC_values + GC_values

non_TC_weights = [1/len(non_TC_values)] * len(non_TC_values)
TC_weights = [1/len(TC_values)] * len(TC_values)
# plot all values on the same histogram
plt.figure(figsize=(10, 6))
plt.hist(non_TC_values, bins=50, alpha=0.5, label='Non-TC', color='blue', weights=non_TC_weights)
#plt.hist(AC_values, bins=50, alpha=0.5, label='AC', color='blue', weights=weights)
#plt.hist(CC_values, bins=50, alpha=0.5, label='CC', color='orange', weights=weights)
#plt.hist(GC_values, bins=50, alpha=0.5, label='GC', color='green', weights=weights)
plt.hist(TC_values, bins=50, alpha=0.5, label='TC', color='red', weights=TC_weights)
#plt.hist(OC_values, bins=50, alpha=0.5, label='OC', color='purple', weights=weights)

plt.xlim(0, 1)  # Set x-axis limit to 0-1
plt.xlabel('Deamination Rate')
plt.ylabel('Frequency')
plt.title('Deamination Rate by Doublet Type ' + label)
plt.legend()

y_min, y_max = plt.ylim()
x_min, x_max = plt.xlim()

plt.text(.6*x_max, 0.9 * y_max, f'AC: {AC_stats[0]:.2f} ({AC_stats[1]:.2f}, {AC_stats[2]:.2f})\n')
plt.text(.6*x_max, 0.8 * y_max, f'CC: {CC_stats[0]:.2f} ({CC_stats[1]:.2f}, {CC_stats[2]:.2f})\n')
plt.text(.6*x_max, 0.7 * y_max, f'GC: {GC_stats[0]:.2f} ({GC_stats[1]:.2f}, {GC_stats[2]:.2f})\n')
plt.text(.6*x_max, 0.6 * y_max, f'TC: {TC_stats[0]:.2f} ({TC_stats[1]:.2f}, {TC_stats[2]:.2f})\n')

#if save == True:
plt.savefig(f'testing/deamination_rate_by_doublet_{label}.pdf', format='pdf')
plt.tight_layout()



plt.show()


'''
'''
# Collect reads (testing only)
read_names = {'CT': [], 'GA': [], 'chimera': [], 'undetermined': [], 'none': []}
for read in pybam.fetch(chrom, start, end):
    if not read.is_secondary and not read.is_supplementary:
        strand = determine_da_strand(read)
        read_names[strand].append(read.query_name)


with open(f'testing/CT_names{label}.txt', 'w') as ctfile:
    for name in read_names['CT']:
        ctfile.write(name + '\n')
with open(f'testing/GA_names{label}.txt', 'w') as gafile:
    for name in read_names['GA']:
        gafile.write(name + '\n')
with open(f'testing/chimera_names{label}.txt', 'w') as chimfile:
    for name in read_names['chimera']:
        chimfile.write(name + '\n')
with open(f'testing/undetermined_names{label}.txt', 'w') as undfile:
    for name in read_names['undetermined']:
        undfile.write(name + '\n')
with open(f'testing/none_names{label}.txt', 'w') as nonefile:
    for name in read_names['none']:
        nonefile.write(name + '\n')


'''