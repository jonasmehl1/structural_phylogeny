# Structural Phylome

From a given taxon sampling and a given seed species compute a normal sequence based phylome and a novel structure based phylome.

## Prokaryotic genomes

See notebook: **workflow/notebooks/prokaryote_sampling.html**

## Important softwares:

* [foldseek](https://github.com/steineggerlab/foldseek)
* [foldtree](https://github.com/DessimozLab/fold_tree)
* [3d-blast](http://3d-blast.life.nctu.edu.tw/dbsas.php)
* [quicktree](https://github.com/khowe/quicktree)
* [PDB_tool](https://github.com/realbigws/PDB_Tool): useful parsers and format converters
* [hhsuite](https://github.com/soedinglab/hh-suite/tree/master): some useful parsers
* [ProstT5](https://github.com/mheinzinger/ProstT5)

# Benchmarking

1. Compute Common ids trees: ml, foldtree, structml and compare RF and verticality
2. compute Foldseek ids: best method among foldtree or structml
3. compute Blast ids
4. Then compare trees from 2 and 3 and see if removing singletons makes the tree better

### Verticality

- Ranger-DTL: plot gene family binned by size (1-5,5-10...) against # of events
- TCS score as used in foldtree (Moi et al 2023)

nw_rename A1A5D9_blast_aa.ufboot new_map > test.nwk
then you can run ALE, however ALE cannot be used with foldtree

### Random Forest classifier of which method could be better

Train a RF model on sequence and structure features with verticality as outcome to predict if a gene family could be better analyzed with seq or struc

# Issues

* Rooting with MAD may be useful, check for what?
* no core cut in foldtree

# TODOs

- [ ] trimming differences may be interesting!
- [ ] treestats file
- [ ] measure if blast singletons have different distributions compared to common hits
- [ ] Model TCS or DL score with alignment and tree statistics, this you could do it in a Rmd

# Done

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
