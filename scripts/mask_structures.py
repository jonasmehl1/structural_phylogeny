#!/usr/bin/env python
import argparse
from Bio import SeqIO
import pandas as pd
import glob


# Parse data
parser = argparse.ArgumentParser(
    description="Mask a fasta file based on structural confidence prediction. Residues with low lddt will be Xs"
)

parser.add_argument("-i", "--in", dest="infile", default="None", 
                    help="input file",required=True)
parser.add_argument("-o", "--out", dest="out", default="None", 
                    help="output file",required=True)
parser.add_argument("-s", "--struct", dest="struct", default="None", 
                    help="structure directory",required=True)
parser.add_argument('-m', '--min', dest='min',
                    help='Minimum lddt', type=int, default=50)

if __name__ == '__main__':
    args = parser.parse_args()
    min_lddt = args.min

    out_seqs = []

    for record in SeqIO.parse(args.infile, "fasta"):
        gene_id = record.name
        json_conf = glob.glob(args.struct+'/*/confidence/'+gene_id+'-confidence_v4.json.gz')[0]

        confidence = pd.read_json(json_conf)
        low_confidence_sites = list(confidence[confidence['confidenceScore']<min_lddt]['residueNumber'])

        for site in low_confidence_sites:
            record.seq = record.seq[:site-1] + 'X' + record.seq[site:]

        out_seqs.append(record)

    SeqIO.write(out_seqs, args.out, "fasta")
