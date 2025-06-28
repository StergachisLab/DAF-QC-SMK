import pandas as pd
import matplotlib.pyplot as plt

metrics_path = snakemake.input.targeting_metrics
regions = snakemake.params.regions
output_plot = snakemake.output.plot

# for testing and troubleshooting
# metrics_path = "/mmfs1/gscratch/stergachislab/bohaczuk/scripts/DAF-QC-SMK/test/summary-tar.tsv.gz"
# regions = ["chr4:3073138-3075853","chr3:179228176-179236561"]

metrics = pd.read_csv(metrics_path, sep="\t")

metrics["proportion_full_length"] = (
    metrics["#_full_length_reads"] / metrics["total_fibers in bam(primary+unmapped)"]
)
metrics["proportion_non_full_length"] = (
    metrics["#_non_full_length_reads"]
    / metrics["total_fibers in bam(primary+unmapped)"]
)
metrics["regions"] = (
    metrics["chrom"]
    + ":"
    + metrics["start"].astype(str)
    + "-"
    + metrics["end"].astype(str)
)

# Add a total row for the whole file
total_full_length = metrics["proportion_full_length"].sum()
total_non_full_length = metrics["proportion_non_full_length"].sum()
metrics.loc[len(metrics), :] = {
    "regions": "All Targets",
    "proportion_full_length": total_full_length,
    "proportion_non_full_length": total_non_full_length,
    "proportion_off_target": 0,  # placeholder for plotting
}

# Calculate proportion Off-target
off_target = 1 - total_full_length - total_non_full_length


plt.figure(figsize=(10, 6))
plt.bar(
    metrics["regions"],
    metrics["proportion_full_length"],
    label="Full Length Reads",
    alpha=0.7,
)
plt.bar(
    metrics["regions"],
    metrics["proportion_non_full_length"],
    label="Non Full-Length Reads",
    alpha=0.7,
    bottom=metrics["proportion_full_length"],
)
plt.bar("Off-Target", off_target, label="Off-target Reads", alpha=0.7, color="gray")
plt.ylabel("Proportion of Reads")
plt.title("Targeting Metrics by Region")
plt.legend(bbox_to_anchor=(1, 0), loc="lower left")
# Remove outside frame
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)

# add text to right side of plot
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim(0, 1)
ystart = ymax
plt.text(0.95 * xmax, ystart, "Proportions of total fibers", fontweight="bold")
for i, row in metrics.iterrows():
    ystart = ystart - 0.05
    target_total_fibers = (
        row["proportion_full_length"] + row["proportion_non_full_length"]
    )
    target_perc_full_length = (
        100
        * row["proportion_full_length"]
        / (row["proportion_full_length"] + row["proportion_non_full_length"])
        if (row["proportion_full_length"] + row["proportion_non_full_length"]) > 0
        else 0
    )
    plt.text(0.95 * xmax, ystart, f"{row['regions']}", fontweight="bold")
    ystart = ystart - 0.15
    plt.text(
        0.95 * xmax,
        ystart,
        (
            f"All: {target_total_fibers:.2f}\n"
            f"Full Length: {row['proportion_full_length']:.2f}\n"
            f"Non Full-Length: {row['proportion_non_full_length']:.2f}\n"
            f"{target_perc_full_length:.0f} % full length at region"
        ),
    )
plt.text(0.95 * xmax, ystart - 0.05, f"Off-target: {off_target:.2f}", fontweight="bold")

plt.savefig(output_plot, format="pdf", bbox_inches="tight")
