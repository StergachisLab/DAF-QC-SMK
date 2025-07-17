# DAF-seq processing and QC pipeline

[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)

This is a Snakemake project template. The `Snakefile` is under `workflow`.


## Install

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.

You can then install the `pixi` environment by cloning this repository and running:

```bash
pixi install
```

## Usage

`pixi` handles the execution of the Snakemake workflows. You can run it with either of the following:


```
pixi shell
snakemake --configfile config/config.yaml
```


```bash
pixi run snakemake ...
```

And if you want to run this Snakemake from another directory you can do so with:

```bash
pixi run --manifest-path /path/to/snakemake/pixi.toml snakemake ...
```

where you update `/path/to/snakemake/pixi.toml` to the path of the `pixi.toml` you cloned.

And in place of `...` use all the normal Snakemake arguments for your workflow.


## Inputs
Table that includes sample name, bam/FASTA path, and targeted regions. For compatibility with pbmarkdup, PacBio bam file inputs should be either unaligned or filtered for primary reads only if a consensus will be generated.
See `config/config.tbl` for a template.
```
sample	file	regs
htt_test	/mmfs1/gscratch/stergachislab/bohaczuk/analysis/chdi-hd/primer-design/25.3.7_vega/HTTLNA_reset.bam	chr4:3073138-3075853,chr3:179228176-179236561
``` 



config.yaml # points to sample table, collects reference genome, and allows for custom parameters.

## Outputs
BAM files
results/{sample_name}/align/{sample_name}.mapped.reads.bam" # bam containing all primary, supplementary, and unaligned reads with PCR duplicates marked (du and ds tags)
results/{sample_name}/align/{sample_name}.decorated.bam" # bam containing full-length reads designated top/bottom with C->T as Y(top strand) and G-A as R (bottom strand)
results/{sample_name}/align/{sample_name}.consensus.bam # bam containing MSA consensus of full length, top/bottom reads with a minimum number of reads specified by dups_required (default:3)

Data files
results/{sample_name}/qc/{sample_name}.detailed_targeting_metrics.tbl.gz # for each region, contains read names designated as full-length, non-full length, unaligned, and supplementary/secondary alignment.
results/{sm}/qc/{sm}.summary_targeting_metrics.tbl # for each region, contains the number of reads in each category


results/{sample_name}/qc/{sample_name}.detailed_seq_metrics.{type}.tbl.gz # contains strand by strand deamination and mutation metrics for full length reads or consensus sequences
results/{sample_name}/qc/{sample_name}.summary_seq_metrics.{type}.tbl.gz # contains proportions of deaminations by region and strand type for full length reads or consensus sequences



Plots




## Acknowledgements

Thank you to Mitchell Vollger for providing the template for this Snakemake workflow and for common functions borrowed from here workflow/rules/common.smk


## TODO
Convert to CRAM files where possible
Note that input of a region with no reads will not fail but will omit that region from plots. Add warning that prints to file if this happens.
Eventually, call nucs, FIRE/accessibility?
test
