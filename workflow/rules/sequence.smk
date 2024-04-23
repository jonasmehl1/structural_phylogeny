outdir=config['outdir']+config['dataset']

rule make_blastdb:
    input:
        fa=expand(config['sequence_dir']+'{code}.fa', code=codes),
        maps=expand(outdir+'/db/{code}.taxid', code=codes)
    output:
        fa=outdir+"/db/all_seqs.fa",
        blast=outdir+"/db/all_seqs_blastdb.pdb",
        mapid=outdir+"/db/taxidmap"
    shell:'''
cat {input.maps} > {output.mapid}
cat {input.fa} > {output.fa}
out_file=$(echo {output.blast} | rev | cut -f 2- -d '.' | rev)
makeblastdb -in {output.fa} -out $out_file -parse_seqids -taxid_map {output.mapid} -dbtype prot
'''

rule make_taxidmap_sp:
    input:
        taxid=rules.make_blastdb.output.mapid,
        groups=config['taxons_file']
    output: outdir+"/db/taxidmap_sps"
    shell:'''
awk 'NR>1' {input.groups} | cut -f2,11 | csvtk join -H -t -f"2;1" {input.taxid} - | \
awk '{{print $1"\\t"$3}}' > {output}
'''

# rule make_taxidmap_ale:
#     input: 
#         taxid=rules.make_blastdb.output.mapid,
#         groups=config['taxons_file']
#     output: outdir+"/db/taxidmap_ale"
#     shell:'''
# awk 'NR>1' {input.groups} | cut -f2,11 | csvtk join -H -t -f"2;1" {input.taxid} - | \
# awk '{{split($1,a,"-"); print $1"\\t"$3"_"a[2]}}' > {output}
# '''

rule make_blastdb_single:
    input:
        fa=rules.make_fastas.output,
        mapid=rules.make_blastdb.output.mapid
    output:
        outdir+"/db/single_dbs/{code}.pdb"
    shell:'''
out_file=$(echo {output} | rev | cut -f 2- -d '.' | rev)
makeblastdb -in {input.fa} -out $out_file -parse_seqids -taxid_map {input.mapid} -dbtype prot
'''

rule blast:
    input:
        q=config['sequence_dir']+'{seed}.fa',
        db=rules.make_blastdb.output.blast,
    output: outdir+"/homology/{seed}_blast.tsv"
    params: config['target_seqs']
    benchmark: outdir+"/benchmarks/homology/{seed}_blast.txt"
    threads: 24
    shell:'''
dbfile=$(echo {input.db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.q} -db $dbfile -out {output} -max_hsps 1 -max_target_seqs {params} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads}
'''


rule blast_allvall:
    input:
        q=config['sequence_dir']+'{seed}.fa',
        t=config['sequence_dir']+'{code}.fa',
        q_db=outdir+"/db/single_dbs/{seed}.pdb",
        t_db=outdir+"/db/single_dbs/{code}.pdb",
    output: 
        q_t=outdir+"/homology/allvall/{seed}_{code}_blast.tsv",
        t_q=outdir+"/homology/allvall/{code}_{seed}_blast.tsv"
    params: config['max_seqs_brh']
    benchmark: outdir+"/benchmarks/homology/{code}_{seed}_blast.txt"
    threads:
        4
    shell:'''
dbfile1=$(echo {input.t_db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.q} -db $dbfile1 -out {output.q_t} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads} -max_hsps 1 -max_target_seqs {params}

dbfile2=$(echo {input.q_db} | rev | cut -f 2- -d '.' | rev)
blastp -query {input.t} -db $dbfile2 -out {output.t_q} \
-outfmt "6 std qcovs qcovhsp qlen slen staxids" -num_threads {threads} -max_hsps 1 -max_target_seqs {params}
'''

rule blast_brh:
    input: 
        a=expand(outdir+"/homology/allvall/{seed}_{code}_blast.tsv", seed=config['seed'], code=codes),
        b=expand(outdir+"/homology/allvall/{code}_{seed}_blast.tsv", seed=config['seed'], code=codes)
    output: outdir+"/homology/{seed}_blast_brh.tsv"
    params: config['eval_brh']
    script: "../scripts/get_BRH.py"

rule aln_aa:
    input: 
        ids=outdir+"/seeds/{seed}/{i}/{i}_{mode}.ids",
        fa=rules.make_blastdb.output.fa
    output:
        seq=outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.seqs",
        aln=outdir+"/seeds/{seed}/{i}/{i}_{mode}_aa.alg"
    log: outdir+"/log/mafft/{seed}_{i}_{mode}_aa.log"
    benchmark: outdir+"/benchmarks/mafft/{seed}_{i}_{mode}_aa.txt"
    threads: 4
    shell:'''
seqkit grep -f {input.ids} {input.fa} > {output.seq}
mafft --auto --thread {threads} {output.seq} > {output.aln} 2> {log}
'''
