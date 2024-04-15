outdir=config['outdir']+config['dataset']

rule get_nums:
    input: config['input_table_file']
    output: outdir+'/db/db_num.tsv'
    shell:'''
echo -e "id\\tseqs\\tstruct\\tdiff" > {output}
for id in $(cut -f1 {input}); do
    seqs=$(wc -l < data/ids/${{id}}.txt)
    struct=$(find data/structures/$id/high_cif -name "*cif.gz" | wc -l)
    echo -e "$id\\t$seqs\\t$struct" | awk '{{print $0"\\t"$2-$3}}'
done >> {output}
'''

rule make_fastas:
    input: config['structure_dir']+'{code}/high_cif'
    output: config['sequence_dir']+'{code}.fa'
    params: config['min_len']
    threads: 8
    shell:'''
python workflow/scripts/cif2fasta.py -i {input} -o {output} -c {threads} -l {params}
'''

rule make_taxidmap:
    input: rules.make_fastas.output
    output: temp(outdir+'/db/{code}.taxid')
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell:'''
grep ">" {input} | sed 's/>//' | awk '{{print $0"\\t"{params.taxid}}}' > {output}
'''

