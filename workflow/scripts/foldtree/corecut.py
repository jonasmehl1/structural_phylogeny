import pandas as pd
import Bio.PDB
import numpy as np
import os
import tqdm
import sys

def extract_core(resdf, outfile,  hitthresh = .8, minthresh = .6, corefolder = 'core_structs/', structfolder = 'structs/', cterfolder = 'cter_structs/' , nterfolder = 'nter_structs/'):
	"""

	Extract a core of structures from a results file

	Parameters
		resdf: path to results file
		outfile: path to output file
		hitthresh: proportion of structures that need to map to a residue for it to be included in the core
		minthresh: if no residues meet the hitthresh, the minimum proportion of structures that need to map to a residue for it to be included in the core
		corefolder: name of folder to output core structures to
		structfolder: name of folder to find structures in
		cterfolder: name of folder to find cter structures in
		nterfolder: name of folder to find nter structures in


	"""
	
	#read all results
	# folder =''.join([ sub + '/' for sub in resdf.split('/')[:-1] ])
	# print(folder)
	resdf = pd.read_table(resdf, header = None)
	resdf.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')

	print('extracting core')
	print('hitthresh: ' + str(hitthresh))
	print('minthresh: ' + str(minthresh))

	queries = resdf['query'].unique()
	nqueries = len(queries)
	print(resdf.head(), nqueries)
	#map hits to each struc
	hits = {}
	#proportion of structures that need to map to a residue fo
	with tqdm.tqdm(total=nqueries) as pbar:
		for i,q in enumerate(queries):
			sub = resdf[resdf['query'] == q]
			hitvec = np.zeros((1, max(sub['qend'])) )
			for idx,r in sub.iterrows():
				hitvec[0,r['qstart']:r['qend']] = hitvec[0,r['qstart']:r['qend']]+1
			hitvec /= nqueries
			core = np.where(hitvec>hitthresh)[1]
			try:
				hits[q]= {'min': np.amin(core), 'max': np.amax(core)}
			except:
				#be more lenient...
				subthresh = np.amax(hitvec)
				print(hitvec, sub, 'be careful, non homologous sequences may have enterred the dataset!')
				if subthresh>=minthresh:
					print('new core threst set at ' + str(subthresh) )
					core = np.where(hitvec>=subthresh)[1]
					hits[q]= { 'min': np.amin(core), 'max': np.amax(core)}
					print(q , 'added')
				else:
					print(q , 'rejected')
			pbar.set_description('processed: %d' % (1 + i))
			pbar.update(1)
	#make core struct folder
	try:
		os.makedirs(corefolder)
	except:
		print(corefolder, 'folder already present')
	# try:
	# 	os.mkdir(cterfolder)
	# except:
	# 	print(cterfolder , 'folder already present')
	# try:
	# 	os.mkdir(folder+nterfolder)
	# except:
	# 	print(nterfolder , 'folder already present')
	
	#parse each struct and output core to folder
	cif_parser = Bio.PDB.MMCIFParser(QUIET=True)  # Parser for .cif files
	pdb_parser = Bio.PDB.PDBParser(QUIET=True)

	with tqdm.tqdm(total=len(hits)) as pbar:
		for i,q in enumerate(hits):

			# Determine file paths
			pdb_file = os.path.join(structfolder, q + ".pdb")
			cif_file = os.path.join(structfolder, q + ".cif")

			# Check which file exists
			if os.path.exists(pdb_file):
				parser = pdb_parser
				struct_file = pdb_file

			elif os.path.exists(cif_file):
				parser = cif_parser
				struct_file = cif_file

			struct = parser.get_structure(q, struct_file)
			outstruct = q.replace('.cif', '').replace('.pdb', '') + ".core.pdb"
			#zero based indexing...
			Bio.PDB.Dice.extract(struct, 'A', hits[q]['min']+1, hits[q]['max']+1, corefolder+outstruct)
			
			# struct_core = Bio.PDB.Dice.extract( struct ,'A' , 0, hits[q]['min']+1  ,folder+nterfolder+q  )
			# #select from max to end
			# struct_core = Bio.PDB.Dice.extract( struct ,'A' , hits[q]['max']+1 , len(struct[0]['A'])  ,folder+cterfolder+q  )

			hits[q]['len'] = [ len(chain) for model in struct for chain in model ][0]

			pbar.set_description('processed: %d' % (1 + i))
			pbar.update(1)


	hitsdf = pd.DataFrame.from_dict(hits, orient='index')
	hitsdf.to_csv(outfile, index_label='query')
	return outfile


