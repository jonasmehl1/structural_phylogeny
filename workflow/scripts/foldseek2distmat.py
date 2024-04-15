# import foldseek2tree
import numpy as np
import pandas as pd
import argparse

# Parse data
parser = argparse.ArgumentParser(
    description="Get a distance matrix from foldseek result"
)

parser.add_argument(
    "-i",
    "--input",
    dest="foldseek",
    action="store",
    default="None",
    help="Input foldseek result",
    required=True
)
parser.add_argument(
    "-o",
    "--output",
    dest="out",
    nargs='+',
    action="store",
    default="None",
    help="output files",
    required=True
)


def distmat_to_txt( identifiers , distmat, outfile):
	'''
	write out a distance matrix in fastme format

	Parameters
	----------
	identifiers : list
		list of identifiers for your proteins
	distmat : np.array  
		distance matrix
	outfile : str   
		path to output file

	'''

	#write out distmat in phylip compatible format
	outstr = str(len(identifiers)) + '\n'
	for i,pdb in enumerate(identifiers):
		outstr += pdb + ' ' + ' '.join( [ "{:.4f}".format(d) for d in list( distmat[i,: ] )  ]  ) + '\n'
	with open(outfile , 'w') as handle:
		handle.write(outstr)
		handle.close()
	return outfile

args = parser.parse_args()

if __name__=="__main__":
    res = pd.read_table(args.foldseek, header = None)
    # print(res.head())

    #get the folder of the input file
    infolder = args.foldseek.split('/')[:-1]
    infolder = ''.join( [i + '/' for i in infolder])+'/'
    # print(infolder)
    res[0] = res[0].map(lambda x :x.replace('-model_v4.cif.gz', ''))
    res[1] = res[1].map(lambda x :x.replace('-model_v4.cif.gz', ''))

    res.columns = 'query,target,fident,lddt,alntmscore'.split(',')
    # res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,lddt,alntmscore,rmsd,prob'.split(',')
    ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
    pos = { protid : i for i,protid in enumerate(ids)}
    kernels = ['fident', 'alntmscore', 'lddt']

    #set kernel columns to float
    for k in kernels:
        res[k] = res[k].astype(float)
    #change nan to 0
    res = res.fillna(0)
    matrices = { k:np.zeros((len(pos), len(pos))) for k in kernels }
    # print(res)

    #calc kernel for tm, aln score, lddt
    for idx,row in res.iterrows():
        for k in matrices:
            matrices[k][pos[row['query']] , pos[row['target']]] += row[k]
            matrices[k][pos[row['target']] , pos[row['query']]] += row[k]

    for i,k in enumerate(matrices):
        matrices[k] /= 2
        matrices[k] = 1-matrices[k]
        # print(matrices[k], np.amax(matrices[k]), np.amin(matrices[k]) )
        # np.save( outfolder + k + '_distmat.npy' , matrices[k])
        distmat_txt = distmat_to_txt(ids , matrices[k] , args.out[i])
