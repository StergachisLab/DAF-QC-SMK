import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


summary_path = snakemake.input.summary_metrics
region = snakemake.params.region
deam_rate = snakemake.output.deam_rate
mut_rate = snakemake.output.mut_rate
strandtype = snakemake.output.strandtype
bias = snakemake.output.bias




# for testing
#summary_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test.summary_seq_metrics.reads.tbl.gz"
#summary_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.summary_metrics.tsv.gz"
#region = "chr4_3073138_3075853"
#deam_rate = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.deamination_rate_histogram.pdf"
#mut_rate = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.mutation_rate_histogram.pdf"
#strandtype = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.read_classification_proportions.pdf"
#bias = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test_manual.deamination_rate_by_doublet_type.pdf"



summary_metrics = pd.read_csv(summary_path, sep='\t', header=0)


for column in ['mutation_rate', 'all_deam_rate', 'AC_deam_rate', 'CC_deam_rate', 'GC_deam_rate', 'TC_deam_rate', 'OC_deam_rate']:
    summary_metrics[column] = summary_metrics[column].apply(lambda x: [float(i) for i in x.split(',')] if isinstance(x, str) else [])

summary_metrics['count'] = summary_metrics['count'].astype(int)


chrom, start, end = region.split("_")
start = int(start)
end = int(end)


label = f"{chrom}:{start}-{end}"
# Filter the summary metrics for the specified region
region_df = summary_metrics[(summary_metrics['chrom'] == chrom) & (summary_metrics['reg_start'] == start) & (summary_metrics['reg_end'] == end)]


if region_df.empty:
    print(f"No data found for region {chrom}:{start}-{end}. Skipping plots.")
    dummy_plot = plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, 'No data available for this region', fontsize=12, ha='center', va='center')
    plt.axis('off')
    plt.savefig(deam_rate, format='pdf')
    plt.savefig(mut_rate, format='pdf')
    plt.savefig(strandtype, format='pdf')
    plt.savefig(bias, format='pdf')
    
else:
    # plot histogram of deamination rates
    
    # Set axis limits
    y_limit = 0.1
    x_limit = 1

    CT_deam_rates = region_df[region_df['strand'] == 'CT']['all_deam_rate'].tolist()[0]
    CT_deam_rates = [float(x) for x in CT_deam_rates if not pd.isna(x)]
    GA_deam_rates = region_df[region_df['strand'] == 'GA']['all_deam_rate'].tolist()[0]
    GA_deam_rates = [float(x) for x in GA_deam_rates if not pd.isna(x)]
    all_deam_rates = CT_deam_rates + GA_deam_rates


    median_deamination = np.median(all_deam_rates)
    deamination_10 = np.percentile(all_deam_rates, 10)
    deamination_90 = np.percentile(all_deam_rates, 90)
    median_CT = np.median(CT_deam_rates)
    median_GA = np.median(GA_deam_rates)


    fig= plt.figure(figsize=(10, 6))
    weights = [1/len(all_deam_rates)] * len(all_deam_rates)  # Normalize histogram
    bin_specs = [x/100 for x in range(0, 101)] # create consistent bin widths
    counts, bins, patches=plt.hist(all_deam_rates, bins=bin_specs, color='blue', alpha=0.7, weights=weights)
    plt.xlim(0, x_limit)
    plt.ylim(0, y_limit)

    
    # Draw red line above histogram if bin exceeds y axis limit
    for i, count in enumerate(counts):
        if count > y_limit:
            bin_left = bins[i]
            bin_right = bins[i+1]
            plt.plot([bin_left, bin_right], [y_limit-y_limit*0.005, y_limit-y_limit*0.005], 
                    color='red', linewidth=2)

    # add median and percentiles to the plot
    plt.axvline(median_deamination, color='black', linestyle='dashed', linewidth=1, label=f'Median: {median_deamination:.2f}')
    plt.axvline(deamination_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {deamination_10:.2f}')
    plt.axvline(deamination_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {deamination_90:.2f}')

    # add text label next to median and quartile lines
    # Add text labels
    plt.text(median_deamination + 0.005*x_limit, y_limit -0.1*y_limit, '50%', rotation=90, va='bottom')
    plt.text(deamination_10 + 0.005*x_limit , y_limit -0.1*y_limit, '10%', rotation=90, va='bottom')
    plt.text(deamination_90 + 0.005*x_limit, y_limit -0.1*y_limit, '90%', rotation=90, va='bottom')

    plt.text(x_limit - x_limit/4, y_limit * 0.9, f'Median: {median_deamination:.2f}', color='black', fontsize=10)
    plt.text(x_limit - x_limit/4, y_limit * 0.8, f'10th Percentile: {deamination_10:.2f}', color='black', fontsize=10)
    plt.text(x_limit - x_limit/4, y_limit * 0.7, f'90th Percentile: {deamination_90:.2f}', color='black', fontsize=10)
    plt.text(x_limit -x_limit/4, y_limit * 0.6, f'Median CT: {median_CT:.2f}', color='red', fontsize=10)
    plt.text(x_limit - x_limit/4, y_limit * 0.5, f'Median GA: {median_GA:.2f}', color='green', fontsize=10)


    plt.title(f'Deamination Rate {label}')
    plt.xlabel('Deamination Rate')
    plt.ylabel('Frequency')

    plt.savefig(deam_rate, format='pdf')
#    plt.show()


    # plot histogram of mutation rates

    # Set axis limits
    y_limit = 0.5
    x_limit = 0.02

    mutation_rates = region_df[region_df['strand'].isin(['CT', 'GA'])]['mutation_rate'].tolist()[0]
    mutation_rates = [float(x) for x in mutation_rates if not pd.isna(x)]
    median_mutation = np.median(mutation_rates)
    mutation_10 = np.percentile(mutation_rates, 10)
    mutation_90 = np.percentile(mutation_rates, 90)


    mutation_rates = [x if x <= x_limit else x_limit for x in mutation_rates]
    # Collapse higher rates into last bin
    fig= plt.figure(figsize=(10, 6))
    weights = [1/len(mutation_rates)] * len(mutation_rates)  # Normalize histogram
    bin_specs = [x/2000 for x in range(0, 41)] # create consistent bin widths
    counts, bins, patches=plt.hist(mutation_rates, bins=bin_specs, color='blue', alpha=0.7, weights=weights)
    plt.xlim(0, x_limit)  # Set x-axis limit to 0-0.02
    plt.ylim(0, y_limit)  # Set y-axis limit to 0-0.2

    # Draw red line above histogram if bin exceeds y axis limit
    for i, count in enumerate(counts):
        if count > y_limit:
            bin_left = bins[i]
            bin_right = bins[i+1]
            plt.plot([bin_left, bin_right], [y_limit-y_limit*0.005, y_limit-y_limit*0.005], 
                    color='red', linewidth=2)

    plt.axvline(median_mutation, color='black', linestyle='dashed', linewidth=1, label=f'Median: {median_mutation:.6f}')
    plt.axvline(mutation_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {mutation_10:.6f}')
    plt.axvline(mutation_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {mutation_90:.6f}')

    # add text label next to median and quartile lines
    plt.text(median_mutation, y_limit + y_limit/80, '50%', rotation=90, va='bottom')
    plt.text(mutation_10, y_limit + y_limit/80, '10%', rotation=90, va='bottom')
    plt.text(mutation_90, y_limit + y_limit/80, '90%', rotation=90, va='bottom')

    plt.text(x_limit - x_limit/4, y_limit * 0.9, f'Median: {median_mutation:.6f}', color='black', fontsize=10)
    plt.text(x_limit - x_limit/4, y_limit * 0.8, f'10th Percentile: {mutation_10:.6f}', color='black', fontsize=10)
    plt.text(x_limit - x_limit/4, y_limit * 0.7, f'90th Percentile: {mutation_90:.6f}', color='black', fontsize=10)



    plt.title('Non-reference Variant Rate ' + label)
    plt.xlabel('Mutation Rate')
    plt.ylabel('Frequency')
    plt.savefig(mut_rate, format='pdf')
#    plt.show()


    # Calculate proportion of CT, GA, chimeric, undetermined, and none reads
    # Set axis limits
    y_limit = 1.0

    total_reads = region_df['count'].sum()
    proportions = {label: region_df[region_df['strand'] == label]['count'].to_list()[0] / total_reads for label in region_df['strand'].unique()}

    ordered_keys = ['CT', 'GA', 'chimera', 'undetermined', 'none']
    ordered_keys = [key for key in ordered_keys if key in proportions]
    values = [proportions[key] for key in ordered_keys if key in proportions]
   
    # Create figure
    fig, ax = plt.subplots(figsize=(10,6))


    # Define colors for each category
    colors = {'CT': 'red',
            'GA': 'green', 
            'chimera': '#ff7f0e',  # orange
            'undetermined': 'blue',  # blue
            'none': '#9467bd'  # purple
            }


    # Create stacked bar
    bottom = 0
    for i, (cat, value) in enumerate(zip(ordered_keys, values)):
        print(f'{cat}: {value:.4%}')  # Debugging output to check proportions
        
        ax.bar(0, value, bottom=bottom, label=f'{cat} ({value:.4%})', width=0.75, color=colors[cat])
        
        bottom += value

    ax.set_ylim(0, 1)
    ax.set_ylabel('Proportion of Reads')
    ax.set_title(f'Read Classification Proportions {label}')
    ax.set_xticks([])
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Get handles and labels and reverse them
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(strandtype, format='pdf')
#    plt.show()



    # plot deamination of doublets
    # Set axis limits
    y_limit = 0.1
    x_limit = 1


    AC_values = region_df[region_df['strand'].isin(['CT','GA'])]['AC_deam_rate'].tolist()[0]
    AC_values = [float(x) for x in AC_values if not pd.isna(x)]
    CC_values = region_df[region_df['strand'].isin(['CT','GA'])]['CC_deam_rate'].tolist()[0]
    CC_values = [float(x) for x in CC_values if not pd.isna(x)]
    GC_values = region_df[region_df['strand'].isin(['CT','GA'])]['GC_deam_rate'].tolist()[0]
    GC_values = [float(x) for x in GC_values if not pd.isna(x)]
    TC_values = region_df[region_df['strand'].isin(['CT','GA'])]['TC_deam_rate'].tolist()[0]
    TC_values = [float(x) for x in TC_values if not pd.isna(x)]
    non_TC_values = AC_values + CC_values + GC_values



    AC_stats = np.median(AC_values), np.percentile(AC_values, 10), np.percentile(AC_values, 90)
    CC_stats = np.median(CC_values), np.percentile(CC_values, 10), np.percentile(CC_values, 90)
    GC_stats = np.median(GC_values), np.percentile(GC_values, 10), np.percentile(GC_values, 90)
    TC_stats = np.median(TC_values), np.percentile(TC_values, 10), np.percentile(TC_values, 90)



    non_TC_weights = [1/len(non_TC_values)] * len(non_TC_values)
    TC_weights = [1/len(TC_values)] * len(TC_values)

    # plot all values on the same histogram
    plt.figure(figsize=(10, 6))
    bin_specs = [x/100 for x in range(0, 101)] # create consistent bin widths
    counts, bins, patches=plt.hist(non_TC_values, bins=bin_specs, alpha=0.5, label='Non-TC', color='blue', weights=non_TC_weights)
    # Draw red line above histogram if bin exceeds y axis limit
    for i, count in enumerate(counts):
        if count > y_limit:
            bin_left = bins[i]
            bin_right = bins[i+1]
            plt.plot([bin_left, bin_right], [y_limit-y_limit*0.005, y_limit-y_limit*0.005], 
                    color='red', linewidth=2)
    #plt.hist(AC_values, bins=50, alpha=0.5, label='AC', color='blue', weights=weights)
    #plt.hist(CC_values, bins=50, alpha=0.5, label='CC', color='orange', weights=weights)
    #plt.hist(GC_values, bins=50, alpha=0.5, label='GC', color='green', weights=weights)
    counts, bins, patches=plt.hist(TC_values, bins=50, alpha=0.5, label='TC', color='red', weights=TC_weights)
    # Draw red line above histogram if bin exceeds y axis limit
    for i, count in enumerate(counts):
        if count > y_limit:
            bin_left = bins[i]
            bin_right = bins[i+1]
            plt.plot([bin_left, bin_right], [y_limit-y_limit*0.005, y_limit-y_limit*0.005], 
                    color='red', linewidth=2)
    #plt.hist(OC_values, bins=50, alpha=0.5, label='OC', color='purple', weights=weights)

    plt.xlim(0, x_limit)  # Set x-axis limit to 0-1
    plt.ylim(0, y_limit)
    plt.xlabel('Deamination Rate')
    plt.ylabel('Frequency')
    plt.title('Deamination Rate by Doublet Type ' + label)
    plt.legend()


    plt.text(.6*x_limit, 0.9 * y_limit, f'AC: {AC_stats[0]:.2f} ({AC_stats[1]:.2f}, {AC_stats[2]:.2f})\n')
    plt.text(.6*x_limit, 0.8 * y_limit, f'CC: {CC_stats[0]:.2f} ({CC_stats[1]:.2f}, {CC_stats[2]:.2f})\n')
    plt.text(.6*x_limit, 0.7 * y_limit, f'GC: {GC_stats[0]:.2f} ({GC_stats[1]:.2f}, {GC_stats[2]:.2f})\n')
    plt.text(.6*x_limit, 0.6 * y_limit, f'TC: {TC_stats[0]:.2f} ({TC_stats[1]:.2f}, {TC_stats[2]:.2f})\n')


    plt.tight_layout()
    plt.savefig(bias, format='pdf')

#    plt.show()
