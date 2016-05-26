import subprocess
import textwrap
import os
import numpy as np
import cPickle as pkl
from multiprocessing import Pool

from featureExtraction import *
# get input sequence
# generate id
# Save sequence to temporary file id.fa 
# run mfold for that temp id
# delete id.fa
# read in id.plot
# delete id.plot
# marginalize id.plot table
# read in id.ann
# delete id.ann
#%start of base pair probability data

def predict_secondary_for_chunks(feature_extractor, outfile_base):
	for feature in feature_extractor.genefeatures:
		print feature
		gene_name = feature.name
		print gene_name
		for chunk_i,chunk_seq in enumerate(feature.chunks):
			_call_rnafold(chunk, chunk_i)
		print feature.name # gene name

def predict_secondary_for_chunks_parallel(fe,
										  outfile_base, 
										  n_processes=1):
	T = 68
	multi_params = []
	for feature in fe.genefeatures:
		gene_name = feature
		for chunk_i,chunk_seq in enumerate(fe.genefeatures[feature].chunks):
			multi_params.append([gene_name, chunk_seq, chunk_i, T])

	if n_processes > 1:
		pool = Pool(n_processes)
		results = pool.map(_call_rnafold, multi_params)
		pool.close()
		pool.join()
		multi_params = ""
	else:
		results = map(_call_rnafold, multi_params)

	# Create map from results
	features = {} # {gene name, [[chunk_id,mfe, contact_prob],...]}
	for result in results:
		gene_name = result[0]
		seqid = result[1]
		mfe = result[2]
		contact_probs = result[3]
		if features.get(gene_name):
			features[gene_name].append([seqid, mfe, contact_probs])
		else:
			features[gene_name] = [seqid, mfe, contact_probs]

	f = open(outfile_base+".npy", "w")
	pkl.dump(features, f)
	f.close()
	return features



def _call_rnafold(args): #T=68 is for ribozero # seq, seqid, T=68
	gene_name = args[0]
	seq = args[1]
	seqid = args[2]
	T = args[3]
	# Save the sequence to a temp file so it can be used by rnafold
	fname = gene_name + "_" + str(seqid)
	print fname + " started"
	with open(fname + ".fa", "w+") as f:
		f.write("> " + fname + "\n")
		for seqline in textwrap.wrap(seq, width=80):
			f.write(seqline+"\n")
	# Call rnafold
	subprocess.call("RNAFold -i %s -o %s -p -T %d" % \
				   (fname + ".fa", fname, float(T)), 
				   shell=True)
	# Remove temp file
	os.remove(fname + ".fa")
	seqlength = len(seq)
	pair_probs = _read_rnafold_dotplot(fname + "_dp.ps", seqlength)
	marg_pair_probs = _marginalize_dotplot(pair_probs)
	mfe = _read_rnafold_mfe(fname + "_" + fname + ".fold")
	os.remove(fname + "_dp.ps")
	os.remove(fname + "_ss.ps")
	os.remove(fname + "_" + fname + ".fold")
	# pnum = _read_mfold_pnum(str(seqid) + ".ann", seqlength)
	# hnum = _read_mfold_hnum(str(seqid) + ".hnum", seqlength)

	print fname + " finished"
	return gene_name, seqid, mfe, marg_pair_probs

def _read_rnafold_mfe(foldfile):
	with open(foldfile,"r") as f:
		lines = f.readlines()
		mfe = str(lines[1].split()[2])[:-1]
	return float(mfe)

def _read_rnafold_dotplot(plotfile, seqlength):
	pair_probs = np.zeros((seqlength, seqlength))
	with open(plotfile,"r") as f:
		start = False 
		for line in f:
			# Keep reading until we reach the start of the dot plot
			if line.startswith("%start of base pair probability data"):
				start = True
				continue
			# If we have not gotten to the start, proceed to next line
			elif start == False:
				continue
			# If we've reached the end, stop
			elif start == True and line.startswith("showpage"):
				break
			line = line.split()
			# Ignore "lbox" for dotplot according to doc
			if line[3] == "lbox":
				continue
			# Process MFE value
			else:
				assert int(line[0])
				i = int(line[0])-1
				j = int(line[1])-1
				p = float(line[2])
				pair_probs[i,j] = p
				pair_probs[j,i] = p
	return pair_probs

def _marginalize_dotplot(dotplot):
	# Return the (col) vector of probs that each base i touches anything else
	marg = np.sum(dotplot, axis=1)
	return marg

# def _read_mfold_pnum(pnumfile, seqlength):
# 	pnum = np.zeros(seqlength)
# 	with open(pnumfile, "r") as f:
# 		header = False #temp
# 		for line in f:
# 			if header == True:
# 				continue
# 			assert int(line[0])
# 			i = int(line[0])
# 			p = float(line[1])
# 			pnum[i] = p
# 	return pnum


# def _read_mfold_hnum(hnumfile, seqlength):
# 	hnum = np.zeros((seqlength, seqlength))
# 	with open(hnumfile,"r") as f:
# 		# There is probably a header
# 		header = False #temp
# 		for line in f:
# 			if header == True:
# 				continue
# 			assert int(line[0])
# 			i = int(line[2])
# 			j = int(line[3])
# 			h = float(line[4])
# 			print i,j,h
# 			hnum[i,j] = h
# 			hnum[j,i] = h
# 	return hnum

def main():
	print "ok I'm running........"
	seq = "GUCUACGGCCAUACCACCCUGAACGCGCCCGAUCUCGUCUGAUCUCGGAAGCUAAGCAGG\
GUCGGGCCUGGUUAGUACUUGGAUGGGAGACCGCCUGGGAAUACCGGGUGCUGUAGGCUU"
	i_seq = 1

	# _call_rnafold(seq, i_seq)

	base_to_rrna = "/Volumes/cycle/rrna/"
	print "Creating feature extractor"
	fe = FeatureExtractor()
	print "Loading gene features"
	fe.load_gene_features() # returns bool; returns true if precomputed
	print "Generating chunked seqs"
	fe.generate_chunked_sequences()
	print "Finished generating chunked seqs"
	# print fe.genefeatures

	print "Predicting secondary structure for all chunked seqs"
	results = predict_secondary_for_chunks_parallel(fe, "temp")
	print "Done!!"
	# generate chunked seqs








if __name__ == '__main__':
	main()

