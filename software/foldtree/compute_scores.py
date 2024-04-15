import treescore
import toytree
import pandas as pd
import subprocess
# import argparse
import numpy as np

from pathlib import Path
# Parse data
# parser = argparse.ArgumentParser(
#     description="Run Ranger DTL to compute most parsimonious number of Duplication and Losses"
# )

# parser.add_argument('--trees', dest='trees',required=True,
# 						help=('tree files for each method'))
# # parser.add_argument("-i", "--in", dest="infile", default="None", 
# #                     help="input file",required=True)
# parser.add_argument("-o", "--out", dest="out", default="None", 
#                     help="output file", required=True)
# parser.add_argument("-t", "--tax", dest="taxidmap", default="None", 
#                     help="taxidmap file", required=True)
# parser.add_argument("-l", "--lineage", dest="lineage", default="None", 
#                     help="lineage file", required=True)
# parser.add_argument("-s", "--sptree", dest="sptree", default="None", 
#                     help="species tree",required=True)
# parser.add_argument("--tmp", dest="tmp", default="/tmp/Ranger_input.txt", 
#                     help="tmp file for Ranger")


# def get_species_name(nodename):
#     taxid = taxidmap[taxidmap['query']==nodename].iloc[0]['species']
#     mnemo = meta[meta['taxid']==taxid].iloc[0]['mnemo']
#     return mnemo

T_score = 2000

if __name__ == '__main__':
    # args = parser.parse_args()
    uniprot_df = pd.read_csv(snakemake.input[1], sep="\t", names=['query', 'mnemo','tax'])
    taxidmap = pd.read_csv(snakemake.input[2], names=['query', 'species'], sep="\t")

    sptree = toytree.tree(snakemake.input[3])
    tmpfile = snakemake.output[0]

    df = []
    with open(snakemake.input[0]) as treefile:
        for line in treefile:
            line = line.strip().split()
            gene = line[0]
            targets = line[1]
            alphabet = line[2]
            tree = toytree.tree(line[3], format = 0)
            tree = treescore.label_leaves(tree, taxidmap, uniprot_df)
            overlap = treescore.getTaxOverlap(tree.treenode)

            leaves = tree.get_tip_labels()
            n_sps = len(set(leaves))
            n_tips = len(leaves)

            with open(tmpfile, 'w') as outfile:
                outfile.write(sptree.write(tree_format=9)+'\n')
                outfile.write(tree.write(tree_format=9)+'\n')

            cmd = './software/Ranger-DTL.linux -i %s -o /dev/stdout -s -T %s -q | grep reco | cut -f2 -d\'(\' | grep -E "[0-9]+" -o' % (tmpfile, T_score)
            result = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
            out, err = result.communicate()
            values = out.decode("utf-8").split('\n')
            dups = values[0]
            losses = values[2]

            # TCS score
            overlap = treescore.getTaxOverlap(tree.treenode)
            taxscore = tree.treenode.score
            # Branch lengths, normalized by total length
            lengths = np.array([node.dist for node in tree.treenode.traverse()])
            total_length = np.sum(lengths)
            lengths /= total_length

            # TCS score of root??
            treescore.getTaxOverlap_root(tree.treenode)
            root_score = treescore.sum_rootscore(tree.treenode)
            # Tip to root distance
            distances = np.array([ node.get_distance(tree.treenode) for node in tree.treenode.get_leaves() ])
            # distances_norm = distances / np.mean(distances)

            row = [gene, targets, alphabet, dups, losses, n_sps, n_tips, 
                    taxscore, total_length, np.mean(lengths), np.mean(distances), 
                    np.var(distances), root_score]
            df.append(row)

    outdf = pd.DataFrame(df, columns=['id', 'targets', 'alphabet', 'dups', 'losses', 'n_sps', 'n_tips',
                                      'score', 'tree_length', 'mean_normalized_length', 'mean_r2t',
                                      'variance_r2t', 'root_score'])

    outdf.to_csv(snakemake.output[1], index=False, sep="\t")