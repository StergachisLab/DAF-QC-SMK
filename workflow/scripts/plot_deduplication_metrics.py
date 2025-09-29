import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


matplotlib.rcParams['pdf.fonttype'] = 42


def mark_cutoff(count, bins, y_limit):
    """
    Draw a red line above the histogram if any bin exceeds the y-axis limit.
    """
    for i, count in enumerate(count):
        if count > y_limit:
            bin_left = bins[i]
            bin_right = bins[i + 1]
            plt.plot(
                [bin_left, bin_right],
                [y_limit - y_limit * 0.005, y_limit - y_limit * 0.005],
                color="red",
                linewidth=2,
            )


# Main script execution
summary_path = snakemake.input.summary_metrics
region = snakemake.params.region
consensus_min_reads = snakemake.params.consensus_min_reads
output_groups = snakemake.output.duplication_groups
output_table = snakemake.output.duplication_table

# For testing
# summary_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/PS01031/qc/reads/PS01031.deduplication_metrics.tbl.gz"
# region = "chr5_34760400_34763141"
# consensus_min_reads = 3

summary_metrics = pd.read_csv(summary_path, sep="\t", header=0)

summary_metrics["du_tags"] = summary_metrics["du_tags"].apply(
    lambda x: x.split(",") if isinstance(x, str) else None
)
summary_metrics["counts"] = summary_metrics["counts"].apply(
    lambda x: [int(i) for i in x.split(",")] if isinstance(x, str) else None
)


chrom, start, end = region.split("_")
start = int(start)
end = int(end)

label = f"{chrom}:{start}-{end}"

# Filter the summary metrics for the specified region
region_df = summary_metrics[
    (summary_metrics["chrom"] == chrom)
    & (summary_metrics["start"] == start)
    & (summary_metrics["end"] == end)
]

if region_df["counts"].to_list()[0] is None:
    print(f"No data found for region {chrom}:{start}-{end}. Skipping plots.")
    plt.figure(figsize=(10, 6))
    plt.text(
        0.5,
        0.5,
        "No data available for this region",
        fontsize=12,
        ha="center",
        va="center",
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_groups, format="pdf")
#    plt.show()
    with open(output_table, "w") as f:
        f.write("")

else:
    # plot histogram of duplicate counts by group and by read
    du_values = region_df["counts"].explode().dropna().tolist()

    # get the sum of du_values
    total_reads = sum(du_values)
    pass_threshold = sum(x for x in du_values if x >= consensus_min_reads)
    perc_dups_threshold = 100 * pass_threshold / total_reads if total_reads > 0 else 0

    # Plot duplicates by group, i.e. only duplicates group size > min cutoff are considered, and each group has equal weight
    dup_only = [x for x in du_values if x >= consensus_min_reads]

    if len(dup_only) == 0:
        print(f"No duplicate groups found with size >= {consensus_min_reads}. Skipping plots.")
        plt.figure(figsize=(10, 6))
        plt.text(
            0.5,
            0.5,
            f"No duplicate groups with size >= {consensus_min_reads}",
            fontsize=12,
            ha="center",
            va="center",
        )
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(output_groups, format="pdf")
        with open(output_table, "w") as f:
            f.write("")

    else:

        group_median_duplicates = np.median(dup_only) if dup_only else 0
        group_duplicates_10 = np.percentile(dup_only, 10) if dup_only else 0
        group_duplicates_90 = np.percentile(dup_only, 90) if dup_only else 0

        x_limit = 500
        y_limit = 0.5
        dup_only = [x if x <= x_limit else x_limit for x in dup_only]

        fig = plt.figure(figsize=(10, 6))
        weights = [1 / len(dup_only)] * len(dup_only)
        bin_specs = range(0, x_limit + 1, 4)  # create consistent bin widths
        counts, bins, patches = plt.hist(
            dup_only, bins=bin_specs, color="blue", alpha=0.7, weights=weights
        )
        plt.xlim(consensus_min_reads, x_limit)
        plt.ylim(0, y_limit)

        # Draw red line above histogram if bin exceeds y axis limit
        mark_cutoff(counts, bins, y_limit)

        plt.axvline(
            group_median_duplicates,
            color="black",
            linestyle="dashed",
            linewidth=1,
            label=f"Median: {group_median_duplicates:.2f}",
        )
        plt.axvline(
            group_duplicates_10,
            color="black",
            linestyle="dashed",
            linewidth=1,
            label=f"10th Percentile: {group_duplicates_10:.2f}",
        )
        plt.axvline(
            group_duplicates_90,
            color="black",
            linestyle="dashed",
            linewidth=1,
            label=f"90th Percentile: {group_duplicates_90:.2f}",
        )

        # Add text labels
        plt.text(
            group_median_duplicates, y_limit + y_limit / 80, "50%", rotation=90, va="bottom"
        )
        plt.text(
            group_duplicates_10, y_limit + y_limit / 80, "10%", rotation=90, va="bottom"
        )
        plt.text(
            group_duplicates_90, y_limit + y_limit / 80, "90%", rotation=90, va="bottom"
        )
        
        plt.title(f"Duplicate groups ≥ {consensus_min_reads}, {label}")
        plt.xlabel("Duplicate Counts")
        plt.ylabel("Proportion of groups")
        plt.savefig(output_groups, format="pdf")
    #    plt.show()

        with open(output_table, "w") as f:
            f.write(
                f"Median duplicates per group: {group_median_duplicates:.2f}\n"
                f"10th Percentile duplicates per group: {group_duplicates_10:.2f}\n"
                f"90th Percentile duplicates per group: {group_duplicates_90:.2f}\n"
                f"Reads in groups of {consensus_min_reads}+ duplicates: {pass_threshold} ({perc_dups_threshold:.2f}%)\n"
                f"Reads in groups of < {consensus_min_reads} duplicates: {total_reads - pass_threshold} ({100 - perc_dups_threshold:.2f}%)\n"
            )


