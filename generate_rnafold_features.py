import subprocess
import textwrap
import os
import numpy as np
import cPickle as pkl
from multiprocessing import Pool
from operator import itemgetter
import sys
import signal

if sys.platform == "darwin":
  base = "/Volumes/cycle/"
else:
  base = "/homes/gws/" + os.getlogin() + "/" # make portable for katies and svens

source_home = base + "rrna/src/"
data_home = base + "rrna/data/"
# rnafold_home = base + "software/ViennaRNA/"

# TODO make everything save to a temporary folder instead of just the rundir

def predict_secondary_for_chunks(feature_extractor, out_folder, T=68):
	for feature in feature_extractor.genefeatures:
		print feature
		gene_name = feature
		print gene_name
		for chunk_i,chunk_seq in enumerate(feature_extractor.genefeatures[feature].chunks):
			_call_rnafold([gene_name, chunk_seq, chunk_i, T, out_folder])
		print feature.name # gene name

# def _write_struct_to_file(fname, q):
# 	f = open(fname, 'wb') 
# 	features = {}
# 	while 1:
# 		gene_name, seqid, mfe, marg_pair_probs = q.get()
# 		if features.get(gene_name):
# 			features[gene_name].append([seqid, mfe, marg_pair_probs])
# 		else:
# 			features[gene_name] = [seqid, mfe, marg_pair_probs]

# 		if m == 'kill':
# 			f.write('writer killed')
# 			break
#         pkl.dump(features, f)
#         # f.flush()
# 	f.close()

 #    	

	# print len(features.keys())
	# pkl.dump(features, outfile)


def predict_secondary_for_chunks_parallel(fe,
										  out_folder, 
										  gene_list=None,
										  n_processes=1,
										  T=68):
	
	# Homo_sapiens.GRCh38.84_gene_names_ribo.npy
	# genefeatures = np.load("/Volumes/cycle/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy")
	if gene_list is None:
		genefeatures = np.load(data_home + "featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy")
	else:
		genefeatures = np.load(gene_list)

	finished_genes = os.listdir(data_home + out_folder)
	finished_genes = [os.path.splitext(i)[0] for i in finished_genes]
	# f = open(outfile_base+".npy", "w")
	# manager = Manager()
	# q = manager.Queue()   
	print len(finished_genes)

	multi_params = []
	features = {}
	genes_not_in_db = []
	for feature in genefeatures:
		gene_name = feature
		if gene_name in finished_genes:
			continue

		try:
			fe.genefeatures[feature]
		except:
			print "Error: couldn't find gene " + gene_name + " in database."
			genes_not_in_db.append(gene_name)
			continue
		print "Gene " + gene_name + " has " + str(len(fe.genefeatures[feature].chunks)) + " chunks."
		for chunk_i,chunk_seq in enumerate(fe.genefeatures[feature].chunks):
			multi_params.append([gene_name, chunk_seq, chunk_i, T, out_folder])

	np.save("genes_not_in_db.npy", genes_not_in_db)
	if n_processes > 1:
		pool = Pool(n_processes, _init_pool)
		# gene_name, seqid, mfe, contact_probs = pool.map(_call_rnafold, multi_params)
		pool.map_async(_call_rnafold, multi_params)
		print "Closing pool"
		pool.close()
		print "Joining pool"
		pool.join()
		multi_params = ""
	else:
		results = map(_call_rnafold, multi_params)

	return features

def _init_pool():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def _call_rnafold(args): #T=68 is for ribozero # seq, seqid, T=68
	gene_name = args[0]
	seq = args[1]
	seqid = args[2] # chunk number
	T = args[3]
	out_folder = args[4]

	# Save the sequence to a temp file so it can be used by rnafold
	fname = gene_name + "_" + str(seqid)
	print fname + " started"
	with open(fname + ".fa", "w+") as f:
		f.write("> " + fname + "\n")
		for seqline in textwrap.wrap(seq, width=80):
			f.write(seqline+"\n")

	# Call rnafold
	print fname + "  calling rnafold"
	subprocess.call("RNAfold -i %s -o %s -p -T %d" % \
				   (fname + ".fa", fname, float(T)), 
				   shell=True)

	# Remove temp file
	os.remove(fname + ".fa")
	seqlength = len(seq)
	pair_probs = _read_rnafold_dotplot(fname + "_dp.ps", seqlength)
	contact_probs = _summarize_dotplot(pair_probs)
	mfe = _read_rnafold_mfe(fname + "_" + fname + ".fold")
	os.remove(fname + "_dp.ps")
	os.remove(fname + "_ss.ps")
	os.remove(fname + "_" + fname + ".fold")

	# Write results
	file_exists = os.path.isfile(data_home + out_folder + "/" + gene_name + ".npy")
	f = ""
	if file_exists:
		f = open(data_home + out_folder + "/" + gene_name + ".npy", "r")
		features = np.load(f)
		f.close()
		features[gene_name].append([seqid, mfe, contact_probs])
		f = open(data_home + out_folder + "/" + gene_name + ".npy", "w")
		pkl.dump(features, f)
		f.close()
	else:
		f = open(data_home + out_folder + "/" + gene_name + ".npy", "w")
		features = {}
		features[gene_name] = [seqid, mfe, contact_probs]
		pkl.dump(features, f)
		f.close()
	matrix_file = open(data_home + out_folder + "/matrices/" + gene_name + "_" + str(seqid) + "_matrix.npy", "w")
	pkl.dump(pair_probs, matrix_file)
	matrix_file.close()
	
	print fname + " finished"
	return gene_name, seqid, mfe, contact_probs

def _read_rnafold_mfe(foldfile):
	with open(foldfile,"r") as f:
		lines = f.readlines()
		mfe = str(lines[1].split()[-1]).replace("(", "").replace(")", "")
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

def _summarize_dotplot(dotplot):
	# Return the (col) vector of probs that each base i touches anything else
	marg = np.mean(dotplot, axis=1)
	return marg

def create_feature_vector(source_folder, out_folder, order):
	# source_folder = "/Volumes/cycle/rrna/data/rnafold/"
	# out_folder = "rnafold"
	# order = "/Volumes/cycle/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy"

	if sys.platform == "darwin":
	  base = "/Volumes/cycle/rrna/"
	else:
	  base = "/homes/gws/kdorosch/rrna/"

	finished_genes = [os.path.splitext(i)[0] for i in os.listdir(source_folder) if os.path.isfile(source_folder + i) and not i.startswith(".")]

	# genes_to_select = np.load(base+"data/featureExtractors/Homo_sapiens.GRCh38.84_gene_names_ribo.npy")
	# print genes_to_select
	# print finished_genes
	# finished_genes = set(finished_genes).intersection(set(genes_to_select))
	# print finished_genes
	
	gene_ordering = np.load(order)
	num_genes = len(gene_ordering)
	gene_lengths = np.zeros(num_genes, dtype=int)
	gene_features = [[] for i in range(num_genes)]
	chunk_len = 50

	for gene in finished_genes:
		# Figure out where this gene is in the list of all genes
		# Finished genes should be in the top N genes, but we read them from 
		#   the dir in abc order so need to re-get the index		
		i = int(list(gene_ordering).index(gene))

		# Read in all chunks for the gene. May or may not be in order
		try:
			f = open(source_folder + gene + ".npy", "r")
			chunks = np.load(f).values()
		except:
			print gene
			continue
		
		sorted(chunks, key=itemgetter(1))
		

		# Make a new list because the first item stupidly isn't a list (bug courtesy of KD)
		ordered_chunks = [chunks[0][:3]]
		chunk_len = max(chunk_len, len(chunks[0][2]))
		ordered_chunks.extend(chunks[0][3:])
		ordered_chunks = sorted(ordered_chunks, key=itemgetter(0))

		print "-----------", gene, "\t# chunks:\t", len(ordered_chunks), "\tindex:\t", i

		# Sanity check to make sure no chunks are missing
		chunk_numbers = [j[0] for j in ordered_chunks]
		mod = chunk_numbers[0]
		chunk_numbers = [j - mod for j in chunk_numbers]
		print range(len(ordered_chunks))[-1], chunk_numbers[-1], np.floor(len(ordered_chunks)-1 / 2.)

		##### TEMPORARY WORKAROUND
		
		if np.floor(len(ordered_chunks)-1 / 2.) != chunk_numbers[-1]:
			print "removing"
			f.close()
			os.remove(source_folder + gene + ".npy")
			continue
		
		assert range(len(ordered_chunks)) == chunk_numbers
		
		gene_lengths[i] = int(len(ordered_chunks))
		gene_features[i] = ordered_chunks

	# assert np.count_nonzero(gene_lengths) == len(finished_genes)
	# for i in range(len(gene_features)-1,-1,-1):
	# 	if gene_lengths[i] > 0:
	# 		break

	gene_lengths = gene_lengths[:i+1]
	print len(gene_lengths)
	gene_features = gene_features[:i+1]

	num_chunks = int(np.sum(gene_lengths))
	row_labels = [[] for j in range(num_chunks)]
	mfe = np.zeros(num_chunks)
	pair_probs = np.zeros((num_chunks, chunk_len))
	to_write = [[] for j in range(num_chunks)]
	index = 0
	for i,length in enumerate(gene_lengths):
		gene_name = gene_ordering[i]
		for j in range(length):
			# print gene_name + "_" + str(j)
			# print i+j
			chunk_number = gene_features[i][0]
			row_labels[index] = [gene_name, gene_name + "_" + str(j), j]
			mfe[index] = gene_features[i][j][1]
			pair_probs[index,:] = gene_features[i][j][2]
			index += 1
			# print row_labels[i+j]
			# print mfe[i+j]

	print i+j
	print index
	print num_chunks
	print len(row_labels)


	# np.save(data_home + out_folder + out_folder[:-1] + "_results/mfe.npy", mfe)
	# np.save(data_home + out_folder + out_folder[:-1] + "_results/pair_probs.npy", pair_probs)
	# np.save(data_home + out_folder + out_folder[:-1] + "_results/row_names.npy", row_labels)
	np.save(out_folder + source_folder.split("/")[-2] + "_mfe.npy", mfe)
	np.save(out_folder + source_folder.split("/")[-2] + "_pair_probs.npy", pair_probs)
	np.save(out_folder + source_folder.split("/")[-2] + "_row_names.npy", row_labels)

	print "Saving to " + out_folder + source_folder.split("/")[-2] + "_*.npy"




# def main():
# 	print "ok I'm running........"
# 	seq = "GUCUACGGCCAUACCACCCUGAACGCGCCCGAUCUCGUCUGAUCUCGGAAGCUAAGCAGG\
# GUCGGGCCUGGUUAGUACUUGGAUGGGAGACCGCCUGGGAAUACCGGGUGCUGUAGGCUU"
# 	i_seq = 1

# 	# _call_rnafold(seq, i_seq)

# 	base_to_rrna = "/Volumes/cycle/rrna/"
# 	print "Creating feature extractor"
# 	fe = FeatureExtractor()
# 	print "Loading gene features"
# 	fe.load_gene_features() # returns bool; returns true if precomputed
# 	print "Generating chunked seqs"
# 	fe.generate_chunked_sequences()
# 	print "Finished generating chunked seqs"
# 	# print fe.genefeatures

# 	print "Predicting secondary structure for all chunked seqs"
# 	results = predict_secondary_for_chunks_parallel(fe, "temp")
# 	print "Done!!"
# 	# generate chunked seqs




# from featureExtraction import *
# import generate_rnafold_features as rnafold
# fe = FeatureExtractor()
# fe.load_chunked_sequences()
# results = rnafold.predict_secondary_for_chunks_parallel(fe, "rnafold_mean", n_processes=40)





# if __name__ == '__main__':
# 	main()

