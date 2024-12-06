# Structural Phylome: A Tool for Structural Phylogenetic Analysis
[![Snakemake](https://img.shields.io/badge/snakemake-≥8-brightgreen.svg)](https://snakemake.github.io)

<!-- vscode-markdown-toc -->
* 1. [Usage](#Usage)
	* 1.1. [Installation](#Installation)
	* 1.2. [Data preparation](#Datapreparation)
	* 1.3. [Phylogeny pipeline](#Phylogenypipeline)
	* 1.4. [Outputs](#Outputs)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

This repo helps running different phylogenetic analyses, including workflows based on protein structures, given some seed sequences or predefined orthogroups. It may be useful to use and benchmark new structural phylogenetics method. It is designed to be easily expandable, so feel free to contribute with code or ideas for us to include!

The results from the first run of the pipeline are reported in this preprint: [Newly developed structure-based methods do not outperform standard sequence-based methods for large-scale phylogenomics](https://www.biorxiv.org/content/10.1101/2024.08.02.606352v1)

##  1. <a name='Usage'></a>Usage

###  1.1. <a name='Installation'></a>Installation

First you need to install snakemake and you can easily do this with this conda command:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake hdf5 snakefmt snakedeploy
```

Then from now one I'd reccomend to manage further dependencies with conda using the `--sdm conda` flag in snakemake. Otherwise you can find in `workflow/envs/` different yaml files to recreate what is needed at each step.

The only dependency not automatically managed in the pipeline is gsutil: follow these instructions for [gsutil installation](https://cloud.google.com/storage/docs/gsutil_install). Gsutil is used to download full UniProt proteomes from [AlphafoldDB](https://alphafold.ebi.ac.uk/).

###  1.2. <a name='Datapreparation'></a>Data preparation

To run the pipeline the user will need to prepare these files:

1. `metadata`: a file specifying the taxon sampling. All the protein structures and sequences from the species included in the file will be downloaded. **IMPORTANT:** the species must be present in UniProt, you can check [here](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README) if your species are present.
2. `seed_file`: A file with either one column with protein IDs of the seed species or two tab separated columns `orthogroup\tprotein_id`.
3. `species_tree`: The corresponding species tree in newick format
4. `configfile`: A `yaml` file with different parameters

The pipeline can be run in two distinct modes: *Phylome* and *OG*. For the first approach, the user only needs to input a list of protein IDs of the seed species (indicated by `seed: UniProt_id` in the `configfile`). Each protein will be aligned to the structures and sequences of the different taxas indicated in the `metadata` file. Alternatively, if the user already has defined orthogroups the homology search step is skipped and the different trees will be computed on these sets.

Importantly, global parameters that are likely to be used across different datasets are in `config/params.yaml`. Note that the values in the first custom yaml are prioritary to the ones in this `params.yaml`! However, it is mandatory that the `configfile` has these fields:


```yaml
# these will be the prefix of the output directory in results/homology
homology_dataset: 'hsap_euka'
# these will be the prefix of the output directory in results/phylogeny
phylo_dataset: 'hsap_1kseeds'
taxids: 'data/input/Hsapopi_set.txt'
species_tree: 'data/sptrees/homo_internal.spTree.nw'
# the seed uniprot id
seed: ['UP000005640']
root: 'Atha'

# this is the number of seed genes to run the pipeline
test_seeds: 'data/seeds/draft_seeds.txt'
```

Once the user has the 4 files, the data downloading can start:

```
snakemake -s workflow/download_data.smk --configfile config/test.yaml -p -j2 --sdm conda
```

This first pipeline is necessary to get all input files. From the input table we can download all the pdbs from google and then consider only those entries with mean average quality > `params["low_confidence"]` value. These proteins will be moved into the `high_cif` folder for each proteome.

You can change the directory where all these data are stored with `params["data_dir"]` parameter but I would just use the default one.

###  1.3. <a name='Phylogenypipeline'></a>Phylogeny pipeline

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

###  1.4. <a name='Outputs'></a>Outputs

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
