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
  # base = "/homes/gws/" + os.getlogin() + "/" # make portable for katies and svens
  base = "/projects/bio/rrna"

source_home = base + "/src/"
data_home = base + "/data/"
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
										  gene_list="/projects/bio/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy",
										  n_processes=1,
										  T=68):
	
	# Homo_sapiens.GRCh38.84_gene_names_ribo.npy
	genefeatures = np.load(gene_list)
	# print genefeatures
	if not os.path.exists(data_home + "/" + out_folder + "/raw"):
		os.makedirs(data_home + "/" + out_folder + "/raw")

	

	finished_genes = os.listdir(data_home + "/" + out_folder)
	finished_genes = [os.path.splitext(i)[0] for i in finished_genes]

	multi_params = []
	features = {}
	genes_not_in_db = []
	for feature in genefeatures:
		gene_name = feature
		# print gene_name
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
	contact_probs_mean = _summarize_dotplot(pair_probs, "mean")
	contact_probs_norm = _summarize_dotplot(pair_probs, "norm")
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
		features[gene_name].append([seqid, mfe, contact_probs_mean, contact_probs_norm])
		f = open(data_home + out_folder + "/" + gene_name + ".npy", "w")
		pkl.dump(features, f)
		f.close()
	else:
		f = open(data_home + out_folder + "/" + gene_name + ".npy", "w")
		features = {}
		features[gene_name] = [seqid, mfe, contact_probs_mean, contact_probs_norm]
		pkl.dump(features, f)
		f.close()

	raw_data_file = data_home + out_folder + "/raw/" + gene_name + "_rawprobs_matrix.npy"
	raw_data_file_exists = os.path.isfile(raw_data_file)
	f = ""
	if raw_data_file_exists:
		f = open(raw_data_file, "r")
		matrices = np.load(f)
		f.close()
		matrices[seqid] = pair_probs
		f = open(raw_data_file, "w")
		pkl.dump(matrices, f)
		f.close()
	else:
		f = open(raw_data_file, "w")
		matrices = {}
		matrices[seqid] = pair_probs
		pkl.dump(matrices, f)
		f.close()
	print fname + " finished"
	# matrix_file = open(data_home + out_folder + "/matrices/" + gene_name + "_" + str(seqid) + "_matrix.npy", "w")
	# pkl.dump(pair_probs, matrix_file)
	# matrix_file.close()
	
	
	return gene_name, seqid, mfe, contact_probs_mean, contact_probs_norm

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

def _summarize_dotplot(dotplot, method):
	# Return the (col) vector of probs that each base i touches anything else
	if method == "norm":
		marg = np.linalg.norm(dotplot, axis=1)
	elif method == "mean":
		marg = np.mean(dotplot, axis=1)
	return marg

def create_feature_vector(source_folder, out_folder, output_base_name=None, order=None):
	# source_folder = "/Volumes/cycle/rrna/data/rnafold/"
	# out_folder = "rnafold"
	# order = "/Volumes/cycle/rrna/data/featureExtractors/Homo_sapiens.GRCh38.84_ordered_gene_names.npy"

	finished_genes = [os.path.splitext(i)[0] for i in os.listdir(source_folder) if os.path.isfile(source_folder + "/" + i) and not i.startswith(".")]

	if order is not None:
		gene_ordering = np.load(order)
		sorted_genes = []
		for gene in gene_ordering:
			if gene in finished_genes:
				sorted_genes.append(gene)
		finished_genes = sorted_genes

	if output_base_name is None:
		output_base_name = source_folder.rstrip("/").split("/")[-1]

	num_genes = len(finished_genes)
	gene_lengths = np.zeros(num_genes, dtype=int)
	gene_features = [[] for i in range(num_genes)]
	detected_chunk_len = None

	error_file_name = "/".join(source_folder.rstrip("/").split("/")[:-1]) + "/create_feature_vector_errors.log"
	print error_file_name
	error_file = open(error_file_name, "w")

	for i, gene in enumerate(finished_genes):
		# Read in all chunks for the gene. May or may not be in order
		try:
			f = open(source_folder + "/" + gene + ".npy", "r")
			chunks = np.load(f).values()
		except:
			print "Error: Could not load chunks for " + gene
			error_file.write("Error: Could not load chunks for " + gene + "\n")
			continue
		
		sorted(chunks, key=itemgetter(1)) # Sort by chunk number

		# Make a new list because the first item stupidly isn't a list (bug courtesy of KD)
		ordered_chunks = [chunks[0][:3]]
		detected_chunk_len = max(detected_chunk_len, len(chunks[0][2]))
		ordered_chunks.extend(chunks[0][3:])
		ordered_chunks = sorted(ordered_chunks, key=itemgetter(0))

		print "-----------", gene, "\t# chunks:\t", len(ordered_chunks), "\tindex:\t", i

		# Sanity check to make sure no chunks are missing
		chunk_numbers = [j[0] for j in ordered_chunks]
		mod = chunk_numbers[0]
		chunk_numbers = [j - mod for j in chunk_numbers]
		# print range(len(ordered_chunks))[-1], chunk_numbers[-1], np.floor(len(ordered_chunks)-1 / 2.)

		##### TEMPORARY WORKAROUND
		
		if np.floor(len(ordered_chunks)-1 / 2.) != chunk_numbers[-1]:
			print "removing"
			error_file.write("Gene " + gene + " is missing chunks.\n")
			f.close()
			# os.remove(source_folder + gene + ".npy")
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
	pair_probs = np.zeros((num_chunks, detected_chunk_len))
	to_write = [[] for j in range(num_chunks)]
	index = 0
	for i,length in enumerate(gene_lengths):
		gene_name = finished_genes[i]
		for j in range(length):
			# print gene_name + "_" + str(j)
			# print i+j
			chunk_number = gene_features[i][0]
			row_labels[index] = [gene_name, gene_name + "_" + str(j), j]
			mfe[index] = gene_features[i][j][1]
			mean_pair_probs[index,:] = gene_features[i][j][2]
			norm_pair_probs[index,:] = gene_features[i][j][3]
			index += 1
			# print row_labels[i+j]
			# print mfe[i+j]

	# print i+j
	print index
	print num_chunks
	print len(row_labels)


	# np.save(out_folder + "/" + output_base_name + "_results/mfe.npy", mfe)
	# np.save(out_folder + "/" + output_base_name + "_results/pair_probs.npy", pair_probs)
	# np.save(out_folder + "/" + output_base_name + "_results/row_names.npy", row_labels)
	print source_folder.split("/")
	np.save(out_folder + "/" + output_base_name + "_mfe.npy", mfe)
	np.save(out_folder + "/" + output_base_name + "_pair_probs_mean.npy", mean_pair_probs)
	np.save(out_folder + "/" + output_base_name + "_pair_probs_norm.npy", norm_pair_probs)
	np.save(out_folder + "/" + output_base_name + "_row_names.npy", row_labels)

	print "Saving to " + out_folder + "/" + output_base_name + "_*.npy"




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

