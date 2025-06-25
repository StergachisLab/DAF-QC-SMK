# Configuration options

See config.yaml for an example


# Required input options
Reference fasta file
```
ref: /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.fa
```

Manifest of input sample(s), must have three white-space separated columns: sample name (sample), input bam/fastq file path (bam), and regions. See config.tbl for an example. The `bam` or `fastq` files should be unaligned. `fastq` is only accepted for ONT sequencing.
```
manifest: config/config.tbl
```

Sequencing platform. Options are 'pacbio' or 'ont'. Default is pacbio. Note that read deduplication and consensus sequence generation are currently only available for PacBio.
```
platform: pacbio
```

# Optional (both platforms)
For each single strand, minimum fraction of C->T + G->A mutations that must be either C->T or G->A for a read to be designated as non-chimeric. 
```
chimera_cutoff: 0.9 # minimum fraction of (C->T|G->A)/(C->T+G->A) required for a top or bottom strand designation (i.e. non-chimeric)
```

# PacBio-specific options
End length to consider for deduplication, sets pbmarkdup `--end-length` flag. Default is 8000
```
dup_end_length: 8000
```

Minimum ID percent to require for deduplication, sets pbmarkdup `--min-id-perc` flag. Default is 99.2
```
dup_min_id_perc: 99.2
```

Generate and output consensus sequence for each read. Default is true.
```
consensus: True
```

The minimum number of reads required for a duplicate group to generate a consensus sequence. Default is 3.
```
consensus_min_reads: 3
```

# ONT-specific options
Indicate whether the input file is `fastq`. Default is bam.
```
is_fastq: False
```
