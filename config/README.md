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
For each read, minimum fraction of C->T + G->A mutations that must be either C->T or G->A for a read to be designated as non-chimeric. Default is 0.9
```
chimera_cutoff: 0.9 # minimum fraction of (C->T|G->A)/(C->T+G->A) required for a top or bottom strand designation (i.e. non-chimeric)
```

For each read, minimum number of C->T + G->A mutations for top/bottom strand designation. Default is 50
```
min_deamination_count: 50 # minimum number of deaminations for top/bottom strand designation
```

+/- end length tolerance (in bp) for classifying full-length reads. Default is 30 bp
```
end_tolerance: 30 # +/- end length tolerance (in bp) for classifying full-length reads
```

Approximate number of filtered reads to output in decorated reads bam for visualization. If this number exceeds the read count post-filtering for full-length and top/bottom reads, all reads will be used. Default is 5000.
```
decorated_samplesize: 5000 # Approximate number of reads to output as decorated reads bam for visualization
```

# PacBio-specific options
End length to consider for deduplication, sets pbmarkdup `--end-length` flag. Default is 0, which uses the whole read
```
dup_end_length: 0
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
