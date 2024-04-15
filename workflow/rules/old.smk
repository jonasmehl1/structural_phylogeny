# rule prepare_ALE:
#     input: config['species_tree']
#     output: outdir+dataset+"/ale/{seed}/sptree.nwk"
#     shell:'''
# Rscript scripts/add_bl.R -s {input} -o {output}    
# '''

# rule run_ALE:
#     input:
#         st=rules.prepare_ALE.output,
#         gt=rules.iqtree.output.ufboot,
#         names_map=rules.make_taxidmap_ale.output
#     output:
#         gt=outdir+dataset+"/ale/{seed}/{i}_{mode}_{alphabet}.nwk",
#         rec=outdir+dataset+"/ale/{seed}/sptree.nwk_{i}_{mode}_{alphabet}.nwk.ale.uml_rec"
#     shell:'''
# nw_rename {input.gt} {input.names_map} > {output.gt}
# ALEobserve {output.gt}
# ALEml_undated {input.st} {output.gt}.ale
# '''

#for seed in results/Hsap_opistho/seeds/UP000005640/*; do nm=$(basename $seed); num=$(comm -12 $seed/${nm}_fs.ids $seed/${nm}_blast.ids | wc -l); echo -e 
#"$nm\t$num\t$(wc -l < $seed/${nm}_blast.ids)\t$(wc -l < $seed/${nm}_fs.ids)"; done >  trees_sampling.txt

### SMK Foldtree

# rule get_distances_foldseek:
#     input: outdir+dataset+"/seeds/{seed}/{i}/{i}_fs.ids"
#     output: outdir+dataset+"/seeds/{seed}/{i}/{i}_fs_allvall.txt"
#     log: outdir+dataset+"/log/fs/{seed}_{i}_struct_fs.log"
#     benchmark: outdir+dataset+"/benchmarks/fs/{seed}_{i}_struct_fs.txt"
#     shell:'''
# seqsdir=$(dirname {input})/seqs
# mkdir -p $seqsdir

# for id in $(cat {input}); do
# ln -s -f -r {structure_dir}/*/cif/${{id}}-model_v4.cif.gz $seqsdir
# done

# foldseek easy-search $seqsdir $seqsdir {output} $TMPDIR/{wildcards.i} \
# --format-output 'query,target,fident,lddt,alntmscore' --exhaustive-search -e inf > {log}

# rm -r $seqsdir
# '''


# rule get_distmat_foldseek:
#     input: outdir+dataset+"/seeds/{seed}/{i}/{i}_fs_allvall.txt"
#     output:
#         outdir+dataset+"/seeds/{seed}/{i}/{i}_fident.txt",
#         outdir+dataset+"/seeds/{seed}/{i}/{i}_alntmscore.txt",
#         outdir+dataset+"/seeds/{seed}/{i}/{i}_lddt.txt"
#     shell:'''
# python ./scripts/foldseek2distmat.py -i {input} -o {output}
# '''

# rule get_foldtree:
#     input: outdir+dataset+"/seeds/{seed}/{i}/{i}_{mattype}.txt",
#     output: outdir+dataset+"/seeds/{seed}/{i}/{i}_{mattype}_fs.nwk"
#     # benchmark: outdir+dataset+"/benchmarks/quicktree/{seed}_{i}_{mattype}.txt"
#     shell:'''
# quicktree -i m {input} | tr -d '\n' | sed -e '$a\\' > {output}
# '''

# MCL if gene clusters

# rule foldseek_mcl:
#     input: 
#         expand(outdir+"foldseek/out/{comb}_fs_filtered.tsv", comb=combinations_codes_str)
#     output: 
#         multiext(outdir+"mcl/foldseek", ".abc", ".tab", ".mci", ".I"+mcl_inflation_str)
#     shell:'''
#     outdir=$(echo {output} | xargs dirname  | sort -u)
#     mkdir -p $outdir

#     cat {input} | cut -f 1,2,11 > ${{outdir}}/foldseek.abc
#     mcxload -abc ${{outdir}}/foldseek.abc --stream-mirror --stream-neg-log10 \
#     -stream-tf 'ceil(200)' -o ${{outdir}}/foldseek.mci -write-tab ${{outdir}}/foldseek.tab
#     mcl ${{outdir}}/foldseek.mci -I {mcl_inflation} -odir $outdir

#     mcxdump -icl ${{outdir}}/out.foldseek.mci.I{mcl_inflation_str} -tabr ${{outdir}}/foldseek.tab -o ${{outdir}}/foldseek.I{mcl_inflation_str}
# '''
