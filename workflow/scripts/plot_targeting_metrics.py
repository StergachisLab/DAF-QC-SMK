import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


matplotlib.rcParams['pdf.fonttype'] = 42


def calculate_metrics(metrics_path):
    """
    Load the metrics file and calculate the proportion of full-length and non-full-length reads
    relative to the total fibers in the BAM file and fiber counts within each region.
    """

    metrics = pd.read_csv(metrics_path, sep="\t")

    metrics["proportion_full_length"] = (
        metrics["#_full_length_reads"]
        / metrics["total_fibers in bam(primary+unmapped)"]
    )

    metrics["proportion_non_full_length"] = (
        metrics["#_non_full_length_reads"]
        / metrics["total_fibers in bam(primary+unmapped)"]
    )

    metrics["proportion_total"] = (
        metrics["proportion_full_length"] + metrics["proportion_non_full_length"]
    )

    metrics["percent_full_length"] = metrics.apply(
        lambda x: x["proportion_full_length"] / x["proportion_total"] * 100
        if x["proportion_total"] > 0
        else 0,
        axis=1,
    )

    metrics["regions"] = (
        metrics["chrom"]
        + ":"
        + metrics["start"].astype(str)
        + "-"
        + metrics["end"].astype(str)
    )

    return metrics


def calculate_summary(metrics):
    """
    Aggregate metrics across all regions.
    """

    total_metrics = {
        "on_target": metrics["proportion_total"].sum(),
        "full_length": metrics["proportion_full_length"].sum(),
        "non_full_length": metrics["proportion_non_full_length"].sum(),
        "full_count": metrics["#_full_length_reads"].sum(),
        "off_target": 1
        - metrics["proportion_full_length"].sum()
        - metrics["proportion_non_full_length"].sum(),
    }

    return total_metrics


# Plot the metrics
metrics_path = snakemake.input.targeting_metrics
regions = snakemake.params.regions
output_plot = snakemake.output.plot
output_metrics = snakemake.output.metrics_txt

# for testing and troubleshooting
# metrics_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/results/htt_test/qc/reads/htt_test.summary_targeting_metrics.tbl"
# regions = ["chr4:3073138-3075853","chr3:179228176-179236561"]

COLORS = {
    "full_length": "blue",
    "non_full_length": "orange",
    "off_target": "gray",
}

metrics = calculate_metrics(metrics_path)

total_metrics = calculate_summary(metrics)

# Create the bar plot
plt.figure(figsize=(10, 6))
plt.bar(
    metrics["regions"],
    metrics["proportion_full_length"],
    label="Full Length Reads",
    alpha=0.7,
    color=COLORS["full_length"],
)

plt.bar(
    metrics["regions"],
    metrics["proportion_non_full_length"],
    label="Non Full-Length Reads",
    alpha=0.7,
    bottom=metrics["proportion_full_length"],
    color=COLORS["non_full_length"],
)

plt.bar(
    "All Targets", total_metrics["full_length"], alpha=0.7, color=COLORS["full_length"]
)

plt.bar(
    "All Targets",
    total_metrics["non_full_length"],
    alpha=0.7,
    color=COLORS["non_full_length"],
    bottom=total_metrics["full_length"],
)

plt.bar(
    "Off-Target",
    total_metrics["off_target"],
    label="Off-target Reads",
    alpha=0.7,
    color=COLORS["off_target"],
)

# Format the plot and save
plt.ylabel("Proportion of Reads")
plt.title("Targeting Metrics by Region")
plt.legend(bbox_to_anchor=(1, 0), loc="lower left")
plt.xticks(rotation=45, ha="right")
plt.gca().spines["top"].set_visible(False)  # removes top frame
plt.gca().spines["right"].set_visible(False)

plt.savefig(output_plot, format="pdf", bbox_inches="tight")

# Save the metrics to a txt file
with open(output_metrics, "w") as f:
    f.write("Proportions of total fibers\n")
    for i, row in metrics.iterrows():
        f.write(
            f"All: {row['proportion_total']:.2f}\n"
            f"Full Length: {row['proportion_full_length']:.2f}\n"
            f"Non Full-Length: {row['proportion_non_full_length']:.2f}\n"
            f"{row['percent_full_length']:.0f} % full length at region\n"
            f"{row['#_full_length_reads']:.0f}X full-length coverage\n\n"
        )
    f.write(
        f"All Targets:\n"
        f"All: {total_metrics['on_target']:.2f}\n"
        f"Full Length: {total_metrics['full_length']:.2f}\n"
        f"Non Full-Length: {total_metrics['non_full_length']:.2f}\n"
        f"{total_metrics['full_count']} full length fibers\n\n"
    )
    f.write(f"Off-target: {total_metrics['off_target']:.2f}")