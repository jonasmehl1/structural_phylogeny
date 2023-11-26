import glob
import pandas as pd

uniprot_genomes = 'data/meta/uniprot_genomes.tsv'
sequence_dir = 'data/fastas/'
structure_dir = 'data/structures/'
outdir = 'results/'

low_confidence = 40

# dataset specific info
seed=config['seed']
files = ['cif', 'pae', 'confidence']

input_table = pd.read_csv(config['input_table_file'], header=None, sep='\t')
input_table.columns = ['uniprot', 'taxid', 'count1', 'count2', 'count3', 'genome', 'source', 'species', 'mnemo']
input_dict = input_table.set_index('uniprot').T.to_dict()

codes = list(input_table['uniprot'])

rule all:
    input:
        expand('data/ids/{code}.txt', code=codes),
        # expand(structure_dir+'{code}/cif', code=codes),
        expand(structure_dir+'{code}/pae', code=codes),
        expand(structure_dir+'{code}/confidence', code=codes),
        expand(structure_dir+'{code}/low_cif', code=codes),
        expand(structure_dir+'{code}/high_cif', code=codes)



rule download_pdbs:
    output:
        directory(structure_dir+'raw_files/{code}')
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell:'''
currdir=$PWD
mkdir -p {output}
cd {output}
gsutil -m cp -r gs://public-datasets-deepmind-alphafold-v4/proteomes/proteome-tax_id-{params.taxid}-* .
cd $currdir
'''

rule download_uniprot_ids:
    output:
        ids='data/ids/{code}.txt'
        # rev_ids='data/ids/{code}_rev.txt'
    params:
        taxid=lambda wcs: str(input_dict[wcs.code]['taxid'])
    shell: '''
set +o pipefail;

kingdom=$(grep {wildcards.code} {uniprot_genomes} | cut -f2 | taxonkit reformat -I 1 -f "{{k}}" | cut -f2)
url=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$kingdom/{wildcards.code}/{wildcards.code}_{params.taxid}.fasta.gz

if curl --head --silent --fail $url 2> /dev/null; then
    echo "Download from reference proteome"
    wget -O - $url | zgrep ">" | cut -f2 -d'|' > {output}
else
    echo "Download from fasta query"
    curl "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A{wildcards.code}%29%29" --output - | \
    zgrep ">" | cut -f2 -d'|' > {output}
fi
'''

rule untar_pdbs:
    input: 
        ids=rules.download_uniprot_ids.output,
        files=rules.download_pdbs.output
    output:
        cif=temp(directory(structure_dir+'{code}/cif')),
        pae=directory(structure_dir+'{code}/pae'),
        conf=directory(structure_dir+'{code}/confidence')
    shell:'''
mkdir -p {output.cif}
mkdir -p {output.pae}
mkdir -p {output.conf}

for tarfile in {input.files}/*tar; do

    ids=$(echo $tarfile | sed 's/.tar//')
    tar --list --file=$tarfile | awk -F- 'NR == FNR {{ a[$0]; next }} $2 in a' {input.ids} /dev/stdin > ${{ids}}.txt
    if [[ $(wc -l <${{ids}}.txt) -ge 1 ]]; then
        grep "model_v4.cif.gz" ${{ids}}.txt > ${{ids}}_cif.txt
        grep "predicted_aligned_error_v4.json.gz" ${{ids}}.txt > ${{ids}}_pae.txt
        grep "confidence_v4.json.gz" ${{ids}}.txt > ${{ids}}_conf.txt
        tar xf $tarfile -C {output.cif} -T ${{ids}}_cif.txt
        tar xf $tarfile -C {output.pae} -T ${{ids}}_pae.txt
        tar xf $tarfile -C {output.conf} -T ${{ids}}_conf.txt
    fi
done
'''


rule move_lowconf:
    input: 
        conf=rules.untar_pdbs.output.conf,
        cif=rules.untar_pdbs.output.cif
    output: 
        stats=structure_dir+'{code}/{code}_mean_plddt.tsv',
        low_dir=directory(structure_dir+'{code}/low_cif'),
        high_dir=directory(structure_dir+'{code}/high_cif')
    shell:'''
mkdir -p {output.low_dir}
mkdir -p {output.high_dir}

> {output.stats}
for struct in {input.conf}/*json.gz; do
    protein=$(basename $struct | cut -f2 -d'-')
    avg_lddt=$(zcat $struct | jq '.confidenceScore |  add/length*1000 | round/1000')
    echo -e "$protein\\t$avg_lddt" >> {output.stats}

    if (( $(echo "$avg_lddt < {low_confidence}" | bc -l) )); then
        cp {input.cif}/AF-${{protein}}-F1-model_v4.cif.gz {output.low_dir}
    else 
        cp {input.cif}/AF-${{protein}}-F1-model_v4.cif.gz {output.high_dir}
    fi
done
'''
