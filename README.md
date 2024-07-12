# Structural Phylome

From a given taxon sampling and a given seed species compute a normal sequence based phylome and a novel structure based phylome.

## Usage

### Data preparation

The first thing must be installing snakemake and you can easily do this with this conda command

```
conda create -c conda-forge -c bioconda -n snakemake snakemake hdf5 snakefmt snakedeploy
```

First of all you need to download the data, to do this you need a config file with X columns. You also need gsutil (follow these instructions for [gsutil installation](https://cloud.google.com/storage/docs/gsutil_install)) as it is quite difficult to install it is the only dependency that you'll need to manage alone. Once you have these you can run:

`snakemake -s workflow/download_data.smk --configfile config/test.yaml -p -j2 --sdm conda`

This first pipeline is necessary to get all input files. From the input table (with uniprot codes and taxid mostly) we can download all the pdbs from google and then consider only those entries with mean average quality > `params["low_confidence"]` value. These will be moved into the `high_cif` folder for each proteome.

You can change the directory where all these data are stored with `params["data_dir"]` parameter but I would just use the default one.

### Phylogeny pipeline

# Developments

## Ideas

- [ ] Use Qmaker with busco single copy or mcl single copy to run qmake and compare resulting model with GTR: NOT NECESSARY ANYMORE
- [ ] One day do the remote homology analysis

## TODOs

- [ ] Make some choices optional -> it would require considerable effort!
- [ ] Better Documentation

## Issues
