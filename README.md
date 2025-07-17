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
test_sample	/my/bam/path/file.bam	chr4:3073138-3075853,chr3:179228176-179236561
``` 

Config file that includes paths to sample table and reference genome and specifies custom parameters.
See `config.config.yaml` for a template.
```
ref: /my/ref/path/genome.fa # path to a reference fasta
manifest: config/config.tbl # table with samples to process
platform: pacbio # sequencing platform, either 'pacbio' or 'ont'. Default is pacbio

# Optional (both platforms)
chimera_cutoff: 0.9 # minimum fraction of (C->T|G->A)/(C->T+G->A) required for a top or bottom strand designation (i.e. non-chimeric)
min_deamination_count: 50 # minimum number of deaminations for top/bottom strand designation
end_tolerance: 30 # +/- end length tolerance (in bp) for classifying full-length reads
decorated_samplesize: 5000 # Approximate number of reads to output as decorated reads bam for visualization

# PacBio-specific options. Input should always be bam files.
dup_end_length: 0 # length to consider for deduplication. 0 will consider the whole read.
dup_min_id_perc: 99.2 # minimum identity percentage required to mark duplicates
consensus: True # whether to generate a consensus sequence for each read. Default is true
consensus_min_reads: 3 # minimum number of reads required to generate a consensus sequence

# ONT-specific options
is_fastq: False # for ONT files only, specifies whether the input files are fastq, otherwise they are assumed to be bam files
```

## Outputs
### BAM files

Aligned bam containing all primary, supplementary, and unaligned reads with PCR duplicates marked (du and ds tags)
```
results/{sample_name}/align/{sample_name}.mapped.reads.bam"
```

bam containing full-length reads designated top/bottom with C->T as Y(top strand) and G->A as R (bottom strand). Strand designation is stored in the ST tag, and deaminated positions are stored in the `DA` tags and `FD` and `LD` tags for first and last deaminated position, respectively. If `decorated_samplesize` is specified in the config file, this will contain a randomly selected sample of full length, top/bottom strand reads.
```
results/{sample_name}/align/{sample_name}.decorated.bam"
```

bam containing MSA consensus of full length, top/bottom reads with a minimum number of reads specified by dups_required (default:3). Consensus read names are constructed from a representative read name (i.e. `pbmarkdup` du tag) with "_consensus" appended. The <pending name> tag indicates the number of reads that were used to construct the consensus. 
```
results/{sample_name}/align/{sample_name}.consensus.bam 
```


### Data files

For each region, contains read names designated as full-length, non-full length, unaligned, and supplementary/secondary alignment.
results/{sample_name}/qc/{sample_name}.detailed_targeting_metrics.tbl.gz 
results/{sm}/qc/{sm}.summary_targeting_metrics.tbl # for each region, contains the number of reads in each category


results/{sample_name}/qc/{sample_name}.detailed_seq_metrics.{type}.tbl.gz # contains strand by strand deamination and mutation metrics for full length reads or consensus sequences
results/{sample_name}/qc/{sample_name}.summary_seq_metrics.{type}.tbl.gz # contains proportions of deaminations by region and strand type for full length reads or consensus sequences


### Plots




## Acknowledgements

Thank you to Mitchell Vollger for providing the template for this Snakemake workflow and for common functions borrowed from here workflow/rules/common.smk


## TODO
Convert to CRAM files where possible
Note that input of a region with no reads will not fail but will omit that region from plots. Add warning that prints to file if this happens.
Eventually, call nucs, FIRE/accessibility?
test
