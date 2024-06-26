# Structural Phylome

From a given taxon sampling and a given seed species compute a normal sequence based phylome and a novel structure based phylome.

## Usage

### Data preparation

The first thing must be installing snakemake and you can easily do this with this conda command

```
conda create  -c conda-forge -c bioconda -n snakemake snakemake hdf5 snakefmt snakedeploy
```

First of all you need to download the data, to do this you need a config file with X columns. You also need gsutil (follow these instructions for [gsutil installation](https://cloud.google.com/storage/docs/gsutil_install)) as it is quite difficult to install it is the only dependency that you'll need to manage alone. Once you have these you can run:

`snakemake -s workflow/download_data.smk --configfile config/test.yaml -p -j2 --sdm conda`

This first pipeline is necessary to get all input files. From the input table (with uniprot codes and taxid mostly) we can download all the pdbs from google and then consider only those entries with mean average quality > `params["low_confidence"]` value. These will be moved into the `high_cif` folder for each proteome.

You can change the directory where all these data are stored with `params["data_dir"]` parameter but I would just use the default one.

### Phylogeny pipeline

## Conda envs

```
conda create -n sp_utils taxonkit seqkit csvtk
conda create -n sp_R r-tidyverse r-phytools r-adephylo r-ape r-phangorn r-patchwork r-wesanderson r-ggdist r-ggrepel r-ggrepel r-ggh4x bioconductor-ggtree bioconductor-biostrings bioconductor-ggmsa
conda create -n sp_homology foldseek blast
conda create -n sp_tree iqtree quicktree mafft trimal hmmer gxx
conda create -n sp_python python biopython pandas numpy fastme toytree tqdm scipy statsmodels cmake
```



# Issues

* RangerDTL 404 installation!

# TODOs

- [ ] Better Documentation
- [ ] Implement Conda

# Done

- [x] treestats file
- [x] Model TCS or DL score with alignment and tree statistics, this you could do it in a Rmd
- [x] SCOP/CATH analysis
- [x] Add union set
- [x] Add plot of distance to seed divided by targets, rooting?
- [x] trimming differences may be interesting!
- [x] measure if blast singletons have different distributions compared to common hits
- [x] update download data smk
- [x] implement ranger inside snakemake rule!
- [x] What to do with astral pro trees?
- [x] work on compare_trees.R and add treefile (model selection and ll information)
- [x] check runstats file
- [x] size astral pro plot by # quartets
- [x] add ALE step, chech if it works and add bl to species trees
- [x] Use Qmaker with busco single copy or mcl single copy to run qmake and compare resulting model with GTR: NOT NECESSARY ANYMORE
- [x] Entropy and how it correlates with lddt, also tool to mask
- [x] What to do if the alignment sucks? add intermediate step
- [x] add rscrpit analyze_runtimes
- [X] Some sequences fail if a sequence is all gap
- [x] add compare RF, TCS etc
- [x] analyze apro in the same script
- [x] Change 3 species
- [x] Substitution matrix
- [x] How to trim a structure? lddt, PAE, entropy
- [X] Add trimming of translated sequences and rule to see if structure ML alignment is shorter than X
- [X] get BRH and compare to sets of homologs
- [X] Query coverage and estimate False Positive rate, maybe with best reciprocal hit approach or similar, see powerpoint by edu
- [x] all v all -> results/dataset/homology/allvall
- [x] Get data frame of query and median lddt to then do the same analysis but filtering on qieries with good structures -> see rule get_lddt
- [x] When downloading, infer if euka or not, do not hard code in url
- [x] Datasets - Get 20 bac and 20 representative archaea from GTDB -> check notebook
- [x] check number of sequences and structures -> added rule get_nums
- [x] check how many repeated targets-hits there are -> fixed max hsps 1 and foldseek already does it i think
- [x] Add benchmark to evaluate ram and time in crucial commands, blast foldseek foldtree
- [x] Need to implement 6 alignment methods of phylome pipeline? if so do a script that does it and outputs the alg and metalig -> no
- [x] Choose thresholds
- [x] Foldseek all v all or compute the tree metrics of subset sequences -> allvalll
- [x] Datasets for tree inference have been selected
- [x] Blast command phylomeDB
