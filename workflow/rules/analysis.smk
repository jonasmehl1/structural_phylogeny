outdir=config['outdir']+config['dataset']
# type of diff matrixed for foldtree
# mattypes = ['fident', 'alntmscore', 'lddt']
modes = ['blast', 'fs', 'common']
methods = ['blast', 'fs']
alphabets = ['aa', '3Di']
alphabets_fident = ['aa', '3Di', 'fident']

rule plot_evalues:
    input:
        table=config['input_table_file'],
        groups=config['taxons_file'],
        taxidmap=rules.make_blastdb.output.mapid,
        blast=rules.blast.output,
        blast_brh=rules.blast_brh.output,
        self_blast=outdir+"/homology/allvall/{seed}_{seed}_blast.tsv",
        fs=rules.foldseek.output,
        fs_brh=rules.foldseek_brh.output,
        self_fs=outdir+"/homology/allvall/{seed}_{seed}_fs.tsv"
    params: 
        eval_both=config['eval_both'],
        max_seqs=config['max_seqs']
    output: 
        eda=outdir+"/plots/{seed}_homology.png",
        saturation=outdir+"/plots/{seed}_saturation.png"
    script: "../scripts/compare_sampling.R"
#     shell:'''
# Rscript workflow/scripts/compare_sampling.R -i {input.table} -m {input.groups} -t {input.taxidmap} \
# -b {input.blast} -f {input.fs} --bb {input.blast_brh} --fb {input.fs_brh} \
# --sb {input.self_blast} --sf {input.self_fs} -o {output.eda} --o2 {output.saturation} -e {params.eval_both} -s {params.max_seqs}
# '''


### Comparisons

rule get_lineage:
    input: config['input_table_file']
    output: temp(outdir+"/reco/{seed}_lineage.tsv")
    shell:'''
cut -f2,9 {input} | taxonkit reformat -I 1 | sed 's/;/,/g' > {output}
'''

rule get_verticality:
    input: 
        trees=rules.check_rooted_trees.output.trees,
        lineage=rules.get_lineage.output,
        taxidmap=rules.make_blastdb.output.mapid,
        sptree=config['species_tree']
    output:
        ranger=temp(outdir+"/reco/{seed}_ranger.tsv"),
        table=outdir+"/reco/{seed}_scores.tsv"
    script:"../../software/foldtree/compute_scores.py"


rule plot_trees:
    input: 
        scores=rules.get_verticality.output.table,
        trees=outdir+"/trees/{seed}_unrooted_trees.txt"
    output:
        outdir+"/plots/{seed}_trees.pdf"
    script: "../scripts/compare_trees.R"


rule get_runstats:
    input: seeds_rooted_trees
    output:
        aln=outdir+"/stats/{seed}_aln.stats",
        time=outdir+"/stats/{seed}_runtime.stats"
    shell:'''
basedir=$(dirname {input} | rev | cut -f2- -d'/' | rev | sort -u)
echo $basedir
seqkit stats $basedir/*/*alg* -a -b -T | cut -f1,4,6,12 > {output.aln}
datasetdir=$(echo $basedir | rev | cut -f3- -d'/' | rev)
echo $datasetdir

for file in $(find $datasetdir/benchmarks/ -type f -name "*txt"); do
    bn=$(basename $file ".txt")
    step=$(basename $(dirname $file))
    seed=$(echo $bn | cut -f1 -d'_')
    gene=$(echo $bn | cut -f2 -d'_')
    method=$(echo $bn | cut -f3 -d'_')
    alphabet=$(echo $bn | cut -f4 -d'_')
    echo -e "$bn\t$step\t$seed\t$gene\t$method\t$alphabet\t$(tail -1 $file)"
done > {output.time}
'''

rule plot_runstats:
    input: 
        aln=rules.get_runstats.output.aln,
        time=rules.get_runstats.output.time
    output: outdir+"/plots/{seed}_runtime.pdf"
    script: "../scripts/analyze_runtimes.R"

rule prepare_astral_pro:
    input:
        table=config['input_table_file'],
        taxidmap=rules.make_blastdb.output.mapid
    output:
        outdir+"/reco/gene_species.map"
    shell:'''
csvtk join -H -t -f 2 -L {input.taxidmap} {input.table} | cut -f1,10 > {output}
'''

rule support_astral_pro:
    input:
        sptree=config['species_tree'],
        genemap=rules.prepare_astral_pro.output,
        trees=rules.check_unrooted_trees.output.trees,
    output: 
        gt=outdir+"/reco/{seed}_{mode}_{alphabet}_apro_input.nwk",
        st=outdir+"/reco/{seed}_{mode}_{alphabet}_apro_support.nwk"
    params: config['root']
    shell:'''
awk '$2=="{wildcards.mode}" && $3=="{wildcards.alphabet}"' {input.trees} | cut -f4 > {output.gt}
astral-pro -c {input.sptree} -a {input.genemap} -u 2 -i {output.gt} -o {output.st} -C --root {params}
'''

rule run_astral_pro:
    input:
        gt=rules.support_astral_pro.output.gt,
        genemap=rules.prepare_astral_pro.output,
        trees=rules.check_unrooted_trees.output.trees
    output: outdir+"/reco/{seed}_{mode}_{alphabet}_apro_sptree.nwk"
    params: config['root']
    log: outdir+"/log/apro/{seed}_{mode}_{alphabet}_apro.log"
    threads: 4
    shell:'''
astral-pro -t {threads} -a {input.genemap} -u 1 -i {input.gt} -o {output} --root {params} 2> {log}
'''

rule plot_astral_pro:
    input:
        sptree=config['species_tree_labels'],
        groups=config['taxons_file'],
        trees=expand(outdir+"/reco/{seed}_{mode}_{alphabet}_apro_support.nwk", seed=config['seed'], mode=modes, alphabet=alphabets_fident)
        # sptrees=expand(outdir+"/reco/{seed}_{mode}_{alphabet}_apro_sptree.nwk", seed=config['seed'], mode=modes, alphabet=alphabets_fident)
    output: outdir+"/plots/{seed}_astral_pro.pdf"
    script: "../scripts/analyze_apro.R"
# reco_dir=$(dirname {input.trees} | sort -u)

