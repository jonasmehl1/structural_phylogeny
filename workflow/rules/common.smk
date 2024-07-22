# outdir=config['outdir']+config['dataset']

rule trim_aln:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    conda: "../envs/sp_tree.yaml"
    shell: '''
trimal -in {input} -out /dev/stdout -gappyout | seqtk seq -A | \
awk '!/^[X-]+$/' | seqtk seq -L 1 -l 60 > {output}
'''
# trimal -in {input} -out {output} -cons {trimal_cons} -gt {trimal_gt}

rule iqtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: 
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
        # treeline=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.treeline",
        # ufboot=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.ufboot",
        iqtree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.iqtree"
    wildcard_constraints:
        model="GTR|LG|3Di"
    params: 
        submat=config['subst_matrix_tree'],
        ufboot=config['UF_boot']
    log: outdir+"/log/iqtree/{seed}_{i}_{mode}_{alphabet}_{model}.log"
    benchmark: outdir+"/benchmarks/iqtree/{seed}_{i}_{mode}_{alphabet}_{model}.txt"
    threads: 4
    conda: "../envs/sp_tree.yaml"
    shell: '''
tree_prefix=$(echo {output.tree} | sed 's/.nwk//')
model={wildcards.model}

if [ "$model" == "GTR" ]; then
    model="GTR20"
elif [ "$model" == "3Di" ]; then
    model="3DI -mdef {params.submat}"
fi

iqtree2 -s {input} --prefix $tree_prefix -B {params.ufboot} -T {threads}  --quiet \
--mem 4G --cmin 4 --cmax 10 --mset $model

mv ${{tree_prefix}}.treefile {output.tree}
mv ${{tree_prefix}}.log {log}
rm -f ${{tree_prefix}}.model.gz ${{tree_prefix}}.splits.nex ${{tree_prefix}}.contree ${{tree_prefix}}.ckp.gz
'''
# --boot-trees

##### FOLDTREE #####

rule foldmason:
    input: 
        fa=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg",
        db=rules.make_foldseekdb.output
    output: 
        fa=temp(outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg.tmp"),
        html=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.html"
    conda: "../envs/sp_utils.yaml"
    shell: '''
cat <(seqkit grep -p "AF-{wildcards.i}-F1" {input.fa}) <(seqkit grep -v -p "AF-{wildcards.i}-F1" {input.fa}) |\
seqkit replace -p $ -r -model_v4.cif > {output.fa}
foldmason msa2lddtreport {input.db} {output.fa} {output.html}
'''

rule foldseek_allvall_tree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
    output: temp(outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt")
    params: config['data_dir']+"structures/"
    log: outdir+"/log/foldseek/{seed}_{i}_{mode}.log"
    benchmark: outdir+"/benchmarks/foldseek/{seed}_{i}_{mode}.txt"
    conda: "../envs/sp_homology.yaml"
    shell:'''
structdir=$(dirname {input})/structs_{wildcards.mode}
mkdir -p $structdir

for id in $(cat {input}); do
zcat {params}*/high_cif/${{id}}-model_v4.cif.gz > $structdir/${{id}}.cif
done

foldseek easy-search $structdir $structdir {output} $TMPDIR/{wildcards.i} \
--format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf > {log}

rm -r $structdir
'''

rule foldseek_distmat:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_allvall.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di_fident.txt"
    conda: "../envs/sp_python.yaml"
    script: "../scripts/foldtree/foldseekres2distmat_simple.py"

rule foldtree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_fident.txt"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FT.nwk"
    wildcard_constraints:
        alphabet="3Di"
    benchmark: outdir+"/benchmarks/foldtree/{seed}_{i}_{mode}_{alphabet}_FT.txt"
    conda: "../envs/sp_tree.yaml"
    shell:'''
quicktree -i m {input} | paste -s -d '' > {output}
'''

rule foldtree_py:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids"
    output:
        distmat=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FTPY.txt",
        tree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_FTPY.nwk"
    wildcard_constraints:
        alphabet="3Di"
    params: config['data_dir']+"structures/"
    log: outdir+"/log/ft/{seed}_{i}_{mode}_{alphabet}_FTPY.log"
    benchmark: outdir+"/benchmarks/ft/{seed}_{i}_{mode}_{alphabet}_FTPY.txt"
    conda: "../envs/sp_python.yaml"
    shell: '''
indir=$(dirname {input})
mkdir -p $indir/structs_{wildcards.mode}

for id in $(cat {input}); do
zcat {params}/*/high_cif/${{id}}-model_v4.cif.gz > $indir/structs_{wildcards.mode}/${{id}}.cif
done

python workflow/scripts/foldtree/foldtree.py -i $indir/structs_{wildcards.mode} -o $indir/{wildcards.i}_{wildcards.mode} \
-t $TMPDIR/{wildcards.i}_{wildcards.mode}_ft --outtree {output.tree} -c $indir/{wildcards.i}_{wildcards.mode}_core \
--corecut --correction --kernel fident > {log}

rm -r $indir/structs_{wildcards.mode}
rm -r $indir/{wildcards.i}_{wildcards.mode}_core
rm {output.distmat}_fastme_stat.txt {output.distmat}.tmp
rm $indir/{wildcards.i}_{wildcards.mode}_allvall.tsv
rm $indir/{wildcards.i}_{wildcards.mode}_core_allvall.tsv
'''
# rm -r $TMPDIR/{wildcards.i}_{wildcards.mode}_ft


rule quicktree:
    input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}.alg.clean"
    output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_QT.nwk"
    wildcard_constraints:
        alphabet="aa"
    params: config["distboot"]
    log: outdir+"/log/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.log"
    benchmark: outdir+"/benchmarks/quicktree/{seed}_{i}_{mode}_{alphabet}_QT.txt"
    conda: "../envs/sp_tree.yaml"
    shell:'''
esl-reformat stockholm {input} | quicktree -boot {params} -in a -out t /dev/stdin | paste -s -d '' > {output} 2> {log}
'''

rule RangerDTL:
    input:
        sptree=config['species_tree'],
        genetree=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}.nwk",
        taxidmap=rules.make_taxidmap_sp.output
    output: 
        full=outdir+"/seeds/{seed}/{i}/{i}_{mode}_{alphabet}_{model}_ranger.txt"
    conda: "../envs/sp_tree.yaml"
    shell:'''
echo -e "$(cat {input.sptree})\\n$(nw_rename {input.genetree} {input.taxidmap})" | \
Ranger-DTL.linux -i /dev/stdin -o {output.full} -T 2000 -q
'''

# rule root_tree:
#     input: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}_{model}.nwk"
#     output: outdir+"/seeds/{seed}/{i}/{i}_{mode}_{algorithm}_{model}.nwk.rooted"
#     log: outdir+"/log/mad/{seed}_{i}_{mode}_{algorithm}_{model}.log"
#     shell: '''
# mad {input} > {log}
# sed -i \'2,$d\' {output}
# '''

