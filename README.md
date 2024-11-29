# Structural Phylome

From a given taxon sampling and a given seed species compute a normal sequence based phylome and a novel structure based phylome.

## Usage

### Installation

The first thing must be installing snakemake and you can easily do this with this conda command

```
conda create -c conda-forge -c bioconda -n snakemake snakemake hdf5 snakefmt snakedeploy
```

Then from now one I'd reccomend to manage further dependencies with conda using the `--sdm conda` flag in snakemake. Otherwise you can find in `workflow/envs/` different yaml files to recreate what is needed at each step.

The only dependency not automatically managed in the pipeline is gsutil: follow these instructions for [gsutil installation](https://cloud.google.com/storage/docs/gsutil_install).

### Data preparation

First of all you need to download the data, to do this you need an input file with 3 columns (Uniprot ID\tTaxid\tmnemonic), see `data/input/Test_set.txt`. You also need a config file where you can specify various parameters, see `config/test.yaml` for example. The analysis in the paper were done using `config/Hsapopi_test.yaml`.

```yaml
# Input file
dataset: 'example'
taxids: 'data/input/Test_set.txt'
species_tree: 'data/sptrees/Test.nwk'
seed: ['UP000000625']
root: 'ecoli'

# this is the number of seed genes to run the pipeline
test_seeds: 15
n_examples: 3
```

Importantly, global parameters that are likely to be used across different datasets are in `config/params.yaml`. Note that the values in the first yaml are prioritary to the ones in this last file!

Once you have these you can simply run:

```
snakemake -s workflow/download_data.smk --configfile config/test.yaml -p -j2 --sdm conda
```

This first pipeline is necessary to get all input files. From the input table we can download all the pdbs from google and then consider only those entries with mean average quality > `params["low_confidence"]` value. These will be moved into the `high_cif` folder for each proteome.

You can change the directory where all these data are stored with `params["data_dir"]` parameter but I would just use the default one.

### Phylogeny pipeline

To run the actual pipeline you just need as input the species tree and the previous uniprot table file. Once you have these you can decide the target sets and the methodological implementations that you'd like to explore:

```yaml
combinations: ["3Di_3Di", "aa_FM_", "aa_LG", "3Di_GTR", "3Di_FT", "3Di_LLM", "3Di_AF"]
#, "3Di_FTPY"]
# Subset of the first one that requires ML 
combinations_ML: ["3Di_3Di", "aa_LG", "3Di_GTR", "3Di_LLM", "3Di_AF"]
```
<!-- modes: ['blast', 'fs', 'common', 'union'] -->

Another important parameter is the `seed`. You can specify one or a list of many (although this feature is more or less untested) and the pipeline will create one output per seed species.

To run the pipeline you can simply run this command and monitor that everything is more or less running fine.

```
snakemake --configfile config/test.yaml -p -j2 -k --sdm conda
```

### Outputs

* `results/{dataset}/trees/{seed}_unrooted_trees.txt`: has 5 columns (gene ID, target set, alphabet, model and tree text). This file is easily parsable in R or Python to do further analyses. You can find some R scripts to do that in `workflow/scripts/` 
* `results/{dataset}/plots/{seed}*pdf`: here you will find various plot that should inform you about the homology search and tree reconstruction quality.
* `results/data`: here you will find very useful files for specific downstream analyses including the fastas and a cif file, a prediction file and the PAE file per protein. You will also find a `gffs` folder where uniprot gffs are stored and a `cath` folder where each proteome has a file that connects protein id to Pfam and Gene3D entries. In the `ids` folder you can also see various phylogenomic useful links in case you would like to benchmark with different homology sets as inferred by OMA, EggNOG or PhylomeDB (among others).
<!-- * `results/{dataset}/homology/{seed}_{method}_brh.tsv`: these files are not used in automatic downstream analyses but may be useful for some further exploration of the results. -->

# TODOs

* update README
* Why DL is different in union and foldseek

# DONE

* implement notebook instead of scripts to plot that can take both eggnoglike or phylomedb.
* Correlation between RF and average pLDDT of a tree
* Run OMA benchmark at eukaryotic level and compute TCS and our metrics: https://github.com/DessimozLab/fold_tree?tab=readme-ov-file#benchmarking-experiments YOU HAVE A PROPER CONDA ENVIRONMENT NAMED foldtree TO DO THIS! (Ran on my workstation as you need internet)
* Check that everything is alright in the two test seeds. 
* Do a proper eggnog and oma benchmark
* Agree on final set of models
* check overlap of blast and foldseek sets with various orthologous sets
* Better representation of CATH and Pfam analysis
* Switch from quicktree to FastME
* 3di + aa: finish rules for partition and add pipeline to the schema -> combination= {comb_part}
* use list of seeds as input.
* replace X with - when masking characters
* File structure is a bit complicated and there is a bit of redundancy, for example rerunning blast everytime or the metadata: done
