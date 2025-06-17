import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

metrics_path = snakemake.input.targeting_metrics
regions = snakemake.params.regions
#metrics_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/htt_test.summary_targeting_metrics.tbl"
#regions = ["chr4:3073138-3075853","chr3:179228176-179236561"]
regions.append("total")



metrics = pd.read_csv(metrics_path, sep='\t')

metrics['proportion_full_length'] = metrics['#_full_length_reads'] / metrics['total_fibers in bam(primary+unmapped)']
metrics['proportion_partial'] = metrics['#_partial_reads'] / metrics['total_fibers in bam(primary+unmapped)']
metrics['regions'] = metrics['chrom'] + ':' + metrics['start'].astype(str) + '-' + metrics['end'].astype(str)

# Add a total row for the whole file
total_full_length = metrics['proportion_full_length'].sum()
total_partial = metrics['proportion_partial'].sum()
metrics.loc[len(metrics), :] = {'regions': 'total',
    'proportion_full_length': total_full_length,
    'proportion_partial': total_partial
}



plt.figure(figsize=(10, 6))
plt.bar(metrics['regions'], metrics['proportion_full_length'], label='Full Length Reads', alpha=0.7)
plt.bar(metrics['regions'], metrics['proportion_partial'], label='Partial Reads', alpha=0.7, bottom=metrics['proportion_full_length'])
plt.ylabel('Proportion of Reads')
plt.title('Targeting Metrics by Region')
plt.legend()
# Remove outside frame
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# add text to right side of plot
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
ystart = 0.9 * ymax
plt.text(0.95*xmax, ystart, "Proportions of total fibers", fontweight='bold')
for i, row in metrics.iterrows():
    ystart= ystart - 0.2
    target_perc_full_length = row['proportion_full_length']/ (row['proportion_full_length'] + row['proportion_partial']) if (row['proportion_full_length'] + row['proportion_partial']) > 0 else 0
    plt.text(0.95*xmax, ystart, (f"Region: {row['regions']}\n"
                                 f"Full Length: {row['proportion_full_length']:.2f}\n"
                                 f"Partial: {row['proportion_partial']:.2f}\n"
                                 f"{target_perc_full_length:.2f} % full length at region"))
    #plt.text(i, row['proportion_partial'] + 0.01, f"{row['proportion_partial']:.2f}", ha='center', va='bottom')

plt.savefig(snakemake.output.plot, format='pdf', bbox_inches='tight')

#plt.show()

# need to add bar for total in whole file