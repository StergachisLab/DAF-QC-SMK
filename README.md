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

`pixi` handles the execution of the Snakemake workflows:

## I have been running it with the following. The options below should also work but have not yet been tested.
```
pixi shell
snakemake --configfile config/config.yaml
```


## Also probably works:
```bash
pixi run snakemake ...
```

And if you want to run this Snakemake from another directory you can do so with:

```bash
pixi run --manifest-path /path/to/snakemake/pixi.toml snakemake ...
```

where you update `/path/to/snakemake/pixi.toml` to the path of the `pixi.toml` you cloned.

And in place of `...` use all the normal Snakemake arguments for your workflow.


## Acknowledgements

Thank you to Mitchell Vollger for providing the template for this Snakemake workflow and for common functions borrowed from here workflow/rules/common.smk


## TODO
Convert to CRAM files where possible
Note that input of a region with no reads will not fail but will omit that region from plots. Add warning that prints to file if this happens.
