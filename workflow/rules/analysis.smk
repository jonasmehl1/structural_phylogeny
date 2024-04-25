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
        eda=outdir+"/plots/{seed}_homology.pdf",
        saturation=outdir+"/plots/{seed}_saturation.pdf"
    script: "../scripts/compare_sampling.R"
#     shell:'''
# Rscript workflow/scripts/compare_sampling.R -i {input.table} -m {input.groups} -t {input.taxidmap} \
# -b {input.blast} -f {input.fs} --bb {input.blast_brh} --fb {input.fs_brh} \
# --sb {input.self_blast} --sf {input.self_fs} -o {output.eda} --o2 {output.saturation} -e {params.eval_both} -s {params.max_seqs}
# '''

rule plot_evalues_trees:
    input:
        table=config['input_table_file'],
        groups=config['taxons_file'],
        taxidmap=rules.make_blastdb.output.mapid,
        blast=rules.blast.output,
        fs=rules.foldseek.output
        # ids=outdir+"/ids/{seed}_common.ids"
    params: 
        eval_both=config['eval_both'],
        coverage=config['coverage'],
        max_seqs=config['max_seqs']
    output: outdir+"/plots/{seed}_singletons.pdf"
    script: "../scripts/compare_sampling_trees.R"


### Comparisons

rule get_lineage:
    input: config['input_table_file']
    output: temp(outdir+"/reco/{seed}_lineage.tsv")
    shell:'''
cut -f2,9 {input} | taxonkit reformat -I 1 | sed 's/;/,/g' > {output}
'''

rule get_verticality:
    input: 
        trees=rules.get_unrooted_trees.output.trees,
        lineage=rules.get_lineage.output,
        taxidmap=rules.make_blastdb.output.mapid
        # sptree=config['species_tree']
    output:
        # ranger=temp(outdir+"/reco/{seed}_ranger.tsv"),
        outdir+"/reco/{seed}_scores.tsv"
    script:"../../software/foldtree/compute_scores.py"


rule plot_trees:
    input: 
        scores=rules.get_verticality.output,
        reco=rules.merge_Ranger.output,
        trees=outdir+"/trees/{seed}_unrooted_trees.txt",
        mltrees=outdir+"/trees/{seed}_mltrees.txt"
    output:
        model=outdir+"/plots/{seed}_trees.pdf",
        reco=outdir+"/plots/{seed}_discordance.pdf"
    script: "../scripts/compare_trees.R"


rule get_runstats:
    input: seeds_unrooted_trees
    output:
        aln=outdir+"/stats/{seed}_aln.stats",
        time=outdir+"/stats/{seed}_runtime.stats"
    threads: 12
    shell:'''
basedir=$(dirname {input} | rev | cut -f2- -d'/' | rev | sort -u)
seqkit stats $basedir/*/*alg* -j {threads} -a -b -T | cut -f1,4,6,12 > {output.aln}

datasetdir=$(echo $basedir | rev | cut -f3- -d'/' | rev)
echo -e "dirname\\tbasename\\ts\\th:m:s\\tmax_rss\\tmax_vms\\tmax_uss\\tmax_pss\\tio_in\\tio_out\\tmean_load\\tcpu_time" > {output.time}
find $datasetdir/benchmarks/ -type f -name "*txt" -printf "%h\\t%f\\t" -exec tail -1 {{}} \; >> {output.time}
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
        trees=rules.get_unrooted_trees.output.trees,
    output: 
        gt=outdir+"/reco/{seed}_{mode}_{alphabet}_{model}_apro_input.nwk",
        st=outdir+"/reco/{seed}_{mode}_{alphabet}_{model}_apro_support.nwk"
    log: outdir+"/log/apro/{seed}_{mode}_{alphabet}_{model}_apro.log"
    params: config['root']
    shell:'''
awk '$2=="{wildcards.mode}" && $3=="{wildcards.alphabet}" && $4=="{wildcards.model}"' {input.trees} | \
cut -f5 > {output.gt}
astral-pro -c {input.sptree} -a {input.genemap} -u 2 -i {output.gt} -o {output.st} -C --root {params} 2> {log}
'''

# rule run_astral_pro:
#     input:
#         gt=rules.support_astral_pro.output.gt,
#         genemap=rules.prepare_astral_pro.output,
#         trees=rules.get_unrooted_trees.output.trees
#     output: outdir+"/reco/{seed}_{mode}_{alphabet}_{model}_apro_sptree.nwk"
#     params: config['root']
#     log: outdir+"/log/apro/{seed}_{mode}_{alphabet}_{model}_apro.log"
#     threads: 4
#     shell:'''
# astral-pro -t {threads} -a {input.genemap} -u 1 -i {input.gt} -o {output} --root {params} 2> {log}
# '''

rule plot_astral_pro:
    input:
        sptree=config['species_tree_labels'],
        groups=config['taxons_file'],
        trees=expand(outdir+"/reco/{seed}_{mode}_{comb}_apro_support.nwk", seed=config['seed'], mode=modes, comb=combinations)
        # sptrees=expand(outdir+"/reco/{seed}_{mode}_{alphabet}_apro_sptree.nwk", seed=config['seed'], mode=modes, alphabet=alphabets_fident)
    output: outdir+"/plots/{seed}_astral_pro.pdf"
    script: "../scripts/analyze_apro.R"
# reco_dir=$(dirname {input.trees} | sort -u)

