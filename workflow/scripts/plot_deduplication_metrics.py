import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

summary_path = snakemake.input.summary_metrics
region = snakemake.params.region
consensus_min_reads = snakemake.params.consensus_min_reads
#output_reads = snakemake.output.duplication_reads
output_groups = snakemake.output.duplication_groups



#summary_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/PS01031/qc/reads/PS01031.deduplication_metrics.tbl.gz"
#region = "chr5_34760400_34763141"
#consensus_min_reads = 3

summary_metrics = pd.read_csv(summary_path, sep='\t', header=0)

summary_metrics['du_tags'] = summary_metrics['du_tags'].apply(lambda x: x.split(',') if isinstance(x, str) else None)
summary_metrics['values'] = summary_metrics['values'].apply(lambda x: [float(i) for i in x.split(',')] if isinstance(x, str) else None)



chrom, start, end = region.split("_")
start = int(start)
end = int(end)

#print(chrom, start, end)


label = f"{chrom}:{start}-{end}"
# Filter the summary metrics for the specified region
region_df = summary_metrics[(summary_metrics['chrom'] == chrom) & (summary_metrics['start'] == start) & (summary_metrics['end'] == end)]


if region_df['values'].to_list()[0] is None:
    print(f"No data found for region {chrom}:{start}-{end}. Skipping plots.")
    dummy_plot = plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, 'No data available for this region', fontsize=12, ha='center', va='center')
    plt.axis('off')
#    plt.savefig(output_reads, format='pdf')
    plt.savefig(output_groups, format='pdf')
#    plt.show()

    
else:

    # plot histogram of duplicate counts by group and by read

    du_values = region_df['values'].explode().dropna().astype(int).tolist()

    # get the sum of du_values
    total_duplicates = sum(du_values)
    dups_threshold = sum(x for x in du_values if x >= consensus_min_reads)
    perc_dups_threshold = 100 * dups_threshold / total_duplicates if total_duplicates > 0 else 0

    # Plot duplicates by group, i.e. only duplicates (group size >1) are considered, and each group has equal weight
    dup_only = [x for x in du_values if x>1]

    group_median_duplicates = np.median(dup_only) if dup_only else 0
    group_duplicates_10 = np.percentile(dup_only, 10) if dup_only else 0
    group_duplicates_90 = np.percentile(dup_only, 90) if dup_only else 0

    upper_limit = 150
    dup_only = [x if x <= upper_limit else upper_limit for x in dup_only]

    fig = plt.figure(figsize=(10, 6))
    weights = [1/len(dup_only)] * len(dup_only)
    plt.hist(dup_only, bins=35, color='blue', alpha=0.7, weights=weights)
    plt.xlim(0, upper_limit)
    plt.ylim(0, 1.0)

    plt.axvline(group_median_duplicates, color='black', linestyle='dashed', linewidth=1, label=f'Median: {group_median_duplicates:.2f}')
    plt.axvline(group_duplicates_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {group_duplicates_10:.2f}')
    plt.axvline(group_duplicates_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {group_duplicates_90:.2f}')

    # add text label next to median and quartile lines
    # Get y-axis limits for positioning text
    y_min, y_max = plt.ylim()
    x_min, x_max = plt.xlim()
    # Add text labels
    plt.text(group_median_duplicates, y_max + y_max/80, '50%', rotation=90, va='bottom')
    plt.text(group_duplicates_10, y_max + y_max/80, '10%', rotation=90, va='bottom')
    plt.text(group_duplicates_90, y_max + y_max/80, '90%', rotation=90, va='bottom')
    plt.text(x_max - x_max/4, y_max * 0.9, f'Median: {group_median_duplicates:.2f}', color='black', fontsize=10)
    plt.text(x_max - x_max/4, y_max * 0.8, f'10th Percentile: {group_duplicates_10:.2f}', color='black', fontsize=10)
    plt.text(x_max - x_max/4, y_max * 0.7, f'90th Percentile: {group_duplicates_90:.2f}', color='black', fontsize=10)
    plt.text(x_max - x_max/4, y_max * 0.6, f'{total_duplicates} reads ({perc_dups_threshold:.0f}%) are\nin groups of {consensus_min_reads}+ duplicates', color='black', fontsize=10)

    plt.title(f'Duplicates by Group {label}')
    plt.xlabel('Duplicate Counts')
    plt.ylabel('Proportion of groups')
    plt.savefig(output_groups, format='pdf')
#    plt.show()



'''

    median_duplicates = np.median(dups_only) if dups_only else 0
    duplicates_10 = np.percentile(dups_only, 10) if dups_only else 0
    duplicates_90 = np.percentile(dups_only, 90) if dups_only else 0


    upper_limit = 150
    dups_only  = [x if x <= upper_limit else upper_limit for x in dups_only]



    fig = plt.figure(figsize=(10, 6))
    weights = dups_only  # Normalize histogram
    plt.hist(dups_only, bins=50, color='blue', alpha=0.7, weights=weights)
    plt.xlim(0, upper_limit)

    plt.axvline(median_duplicates, color='black', linestyle='dashed', linewidth=1, label=f'Median: {median_duplicates:.2f}')
    plt.axvline(duplicates_10, color='black', linestyle='dashed', linewidth=1, label=f'10th Percentile: {duplicates_10:.2f}')
    plt.axvline(duplicates_90, color='black', linestyle='dashed', linewidth=1, label=f'90th Percentile: {duplicates_90:.2f}')

    # add text label next to median and quartile lines
    # Get y-axis limits for positioning text
    y_min, y_max = plt.ylim()

    x_min, x_max = plt.xlim()

    # Add text labels
    plt.text(median_duplicates, y_max + y_max/80, '50%', rotation=90, va='bottom')
    plt.text(duplicates_10, y_max + y_max/80, '10%', rotation=90, va='bottom')
    plt.text(duplicates_90, y_max + y_max/80, '90%', rotation=90, va='bottom')

    plt.text(x_max - x_max/4, y_max * 0.9, f'Median: {median_duplicates:.2f}', color='black', fontsize=10)
    plt.text(x_max - x_max/4, y_max * 0.8, f'10th Percentile: {duplicates_10:.2f}', color='black', fontsize=10)
    plt.text(x_max - x_max/4, y_max * 0.7, f'90th Percentile: {duplicates_90:.2f}', color='black', fontsize=10)

    plt.title(f'Duplicate Counts {label}')
    plt.xlabel('Duplicate Counts')
    plt.ylabel('Frequency')

#    plt.savefig(output_path, format='pdf')
    plt.show()

'''