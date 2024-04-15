outdir=config['outdir']+config['dataset']

rule trim_aln:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    shell: '''
trimal -in {input} -out /dev/stdout -gappyout | seqtk seq -A | \
awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
'''
# trimal -in {input} -out {output} -cons {trimal_cons} -gt {trimal_gt}

rule iqtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.treefile",
        treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.treeline",
        ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.ufboot"
    params: 
        submat=config['subst_matrix_tree'],
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}.txt"
    threads: 4
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.treefile//')

if [ {wildcards.alphabet} = "3Di" ]; then
    model="--mset 3DI -mdef {params.submat}"
else
    model="-m LG+G4"
fi

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads} --boot-trees --quiet \
--mem 4G --cmin 4 --cmax 10 $model

best_model=$(grep "Model of substitution:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | sed 's/ //')
loglik=$(grep "Log-likelihood of the tree:" ${{tree_prefix}}.iqtree | cut -f2 -d':' | cut -f2 -d' ')

echo -e "{wildcards.i}\t$best_model\t$loglik\t$(cat {output.tree})" > {output.treeline}

mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.iqtree ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''


##### FOLDTREE #####

rule run_foldtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
    output:
        distmat=outdir+"/seeds/{seed}/{i}/{i}_{mode}_fident.txt",
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_fident.treefile"
    params: config['structure_dir']
    log: outdir+"/log/ft/{seed}_{i}_{mode}.log"
    benchmark: outdir+"/benchmarks/ft/{seed}_{i}_{mode}.txt"
    shell: '''
indir=$(dirname {input})
mkdir -p $indir/structs_{wildcards.mode}

for id in $(cat {input}); do
zcat {params}/*/high_cif/${{id}}-model_v4.cif.gz > $indir/structs_{wildcards.mode}/${{id}}.cif
done

python ./software/foldtree/foldtree.py -i $indir/structs_{wildcards.mode} -o $indir/{wildcards.i}_{wildcards.mode} \
-t $TMPDIR/{wildcards.i}_{wildcards.mode}_ft -c $indir/{wildcards.i}_{wildcards.mode}_core \
--corecut --correction --kernel fident > {log}

rm -r $indir/structs_{wildcards.mode}
rm -r $TMPDIR/{wildcards.i}_{wildcards.mode}_ft
rm -r $indir/{wildcards.i}_{wildcards.mode}_core
rm $indir/{wildcards.i}_{wildcards.mode}_fident.txt_fastme_stat.txt $indir/{wildcards.i}_{wildcards.mode}_fident.txt.tmp
rm $indir/{wildcards.i}_{wildcards.mode}_allvall.tsv
rm $indir/{wildcards.i}_{wildcards.mode}_core_allvall.tsv
'''


rule root_tree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}.treefile"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}.treefile.rooted"
    log: outdir+"/log/mad/{seed}_{i}_{mode}_{algorithm}.log"
    shell: '''
./software/mad {input} > {log}
sed -i \'2,$d\' {output}
'''