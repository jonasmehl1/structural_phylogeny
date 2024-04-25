outdir=config['outdir']+config['dataset']

input_table = pd.read_csv(config['input_table_file'], header=None, sep='\t')
input_table.columns = ['uniprot', 'taxid', 'count1', 'count2', 'count3', 'genome', 'source', 'species', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])

rule make_foldseekdb:
    input:
        all_struct=expand(config['structure_dir']+'{code}/high_cif', code=codes),
    output:
        outdir+"/db/all_seqs_fsdb"
    shell:'''
structsdir=$(dirname {output})/structs
mkdir -p $structsdir

for id in {input}; do
bn=$(dirname $id | rev | cut -f1 -d'/' | rev)

if [ -L $structsdir/$bn ]; then
    unlink $structsdir/$bn
fi
ln -s -f -r $id $structsdir/$bn
done

foldseek createdb $structsdir {output}
'''

rule make_foldseekdb_seq:
    input: rules.make_foldseekdb.output
    output: outdir+"/db/all_seqs_fsdb_ss.fa"
    shell:'''
foldseek lndb {input}_h {input}_ss_h
foldseek convert2fasta {input}_ss {output}
sed -i 's/-model_v4.cif.gz//g' {output}
'''

rule make_foldseekdb_single:
    input: config['structure_dir']+'{code}/high_cif',
    output: outdir+"/db/single_dbs/{code}_fsdb"
    threads:
        8
    shell:'''
foldseek createdb --threads {threads} {input} {output}
'''

rule foldseek:
    input:
        q=config['structure_dir']+'{seed}/high_cif',
        db=rules.make_foldseekdb.output
    output: 
        outdir+"/homology/{seed}_fs.tsv"
    params: config['target_seqs']
    log: outdir+"/log/homology/{seed}_fs.log"
    benchmark: outdir+"/benchmarks/homology/{seed}_fs.txt"
    threads:
        24
    shell:'''
foldseek easy-search {input.q} {input.db} {output} $TMPDIR/{wildcards.seed} --threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif.gz//g' {output}
'''


rule foldseek_allvall:
    input:
        q=config['structure_dir']+'{seed}/high_cif',
        t=config['structure_dir']+'{code}/high_cif'
    output:
        q_t=outdir+"/homology/allvall/{seed}_{code}_fs.tsv",
        t_q=outdir+"/homology/allvall/{code}_{seed}_fs.tsv"
    params: config['max_seqs_brh']
    log: outdir+"/log/homology/{code}_{seed}_fs.log"
    benchmark: outdir+"/benchmarks/homology/{code}_{seed}_fs.txt"
    threads:
        4
    shell:'''
foldseek easy-search {input.q} {input.t} {output.q_t} $TMPDIR/{wildcards.seed}_{wildcards.code} \
--threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif.gz//g' {output.q_t}

foldseek easy-search {input.t} {input.q} {output.t_q} $TMPDIR/{wildcards.code}_{wildcards.seed} \
--threads {threads} --max-seqs {params} \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,rmsd,prob,qcov,tcov > {log}
sed -i 's/-model_v4.cif.gz//g' {output.t_q}
'''

rule foldseek_brh:
    input: 
        a=expand(outdir+"/homology/allvall/{seed}_{code}_fs.tsv", seed=config['seed'], code=codes),
        b=expand(outdir+"/homology/allvall/{code}_{seed}_fs.tsv", seed=config['seed'], code=codes)
    output:
        outdir+"/homology/{seed}_fs_brh.tsv"
    params: config['eval_brh']
    script: "../scripts/get_BRH.py"


# If foldseek, retrieve sequences from "translated version"
rule aln_3Di:
    input:
        ids=outdir+"/seeds/{seed}/{i}/{i}_{mode}.top",
        fa=rules.make_foldseekdb_seq.output
    output:
        seq=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.seqs",
        masked=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.masked.seqs",
        aln=outdir+"/seeds/{seed}/{i}/{i}_{mode}_3Di.alg"
    params: 
        structdir=config['structure_dir'],
        submat=config['subst_matrix'],
        min_lddt=config['min_lddt']
    log: outdir+"/log/mafft/{seed}_{i}_{mode}_3Di.log"
    benchmark: outdir+"/benchmarks/mafft/{seed}_{i}_{mode}_3Di.txt"
    threads: 4
    shell:'''
seqkit grep -f {input.ids} {input.fa} > {output.seq}
python workflow/scripts/mask_structures.py -i {output.seq} -o {output.masked} -m {params.min_lddt} -s {params.structdir}
mafft --auto --thread {threads} --aamatrix {params.submat} {output.masked} > {output.aln} 2> {log}
'''