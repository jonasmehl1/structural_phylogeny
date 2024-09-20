outdir=config['outdir']+config['dataset']

input_table = pd.read_csv(config['taxids'], sep='\t')
input_table.columns = ['uniprot', 'taxid', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])

rule download_up_meta:
    output: outdir+'/meta/uniprot_genomes.tsv'
    shell: '''
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/STATS -O - | \
tail -n +16 | sed 's/ //' > {output}
'''

rule get_taxon_file:
    input:
        up=rules.download_up_meta.output,
        sps=config["taxids"]
    output: outdir+'/meta/'+config["dataset"]+'_uniprot_genomes.tsv'
    conda: "../envs/sp_utils.yaml"
    shell: '''
csvtk join -t -f "Tax_ID,Proteome_ID" {input.sps} {input.up} | cut -f1,2,3,7,9 | \
taxonkit reformat -I 2 -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" | awk 'NR>1' | \
csvtk add-header -t -n Proteome_ID,Tax_ID,mnemo,Assembly,Species_name,Kingdom,Phylum,Class,Order,Family,Genus,Species > {output}
'''

rule get_nums:
    input: 
        seqs=expand(config['data_dir']+'fastas/{code}.fa', code=codes),
        structs=expand(config['data_dir']+'structures/{code}/high_cif', code=codes)
    output:
        num=outdir+'/db/db_num.tsv',
        seqs=temp(outdir+'/db/db_num_seqs.tsv'),
        structs=temp(outdir+'/db/db_num_structs.tsv')
    conda: "../envs/sp_utils.yaml"
    shell:'''
find {input.structs} -name "*cif.gz" | rev | cut -d/ -f3 | rev | \
sort | uniq -c | awk -v OFS='\\t' '{{$1=$1;print $2,$1}}' > {output.structs}
echo {input.seqs} | tr ' ' '\\n' | seqkit stats -T -b -X - | sed 's/.fa//g' | \
cut -f1,4 | awk 'NR>1' > {output.seqs}

csvtk join -H -t -f1 {output.structs} {output.seqs} > {output.num}
'''

rule make_fastas:
    input: config['data_dir']+'structures/{code}/high_cif'
    output: config['data_dir']+'fastas/{code}.fa'
    params: config['min_len']
    threads: 8
    conda: "../envs/sp_python.yaml"
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

rule get_gff:
    output: config['data_dir']+'gffs/{code}.gff'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?format=gff&query=%28%28proteome%3A{wildcards.code}%29%29" \
-O /dev/stdout | awk 'NF' > {output}  
'''

rule get_CATH:
    output: config['data_dir']+'cath/{code}_cath.tsv'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_pfam%2Cxref_gene3d&format=tsv&query=%28%28proteome%3A{wildcards.code}%29%29" -O {output}  
'''

rule get_phygeno:
    output: config['data_dir']+'ids/{code}_ids.tsv'
    shell: '''
wget "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cxref_genetree%2Cxref_hogenom%2Cxref_inparanoid%2Cxref_phylomedb%2Cxref_orthodb%2Cxref_oma%2Cxref_treefam%2Cxref_eggnog&format=tsv&query=%28%28proteome%3A{wildcards.code}%29%29" \
-O {output}  
'''
