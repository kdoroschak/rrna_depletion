import numpy as np
import os

def combine_matrices(indir, outdir):
	genes = set()
	gene_files = {}
	chunks = {} # key: gene, value: {chunk #, matrix}

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	print "Getting unique genes..."
	if os.path.exists(outdir + "/unique_genes.npy") \
	   and os.path.exists(outdir + "/unique_genes_matrix_filenames.npy"):
	   print "Yay! The files exist so we don't have to read everything in again from scratch."
	   genes = np.load(outdir + "/unique_genes.npy").item()
	   gene_files = np.load(outdir + "/unique_genes_matrix_filenames.npy").item()
	else:
		print "Reading in matrix files (this will take awhile)..."
		matrix_files = os.listdir(indir)
		matrix_files = [i for i in matrix_files if os.path.isfile(indir + "/" + i) and not i.startswith(".") and "unique" not in i]
		for matrix_file in matrix_files:
			print matrix_file
			gene = matrix_file.split("_")[-2]
			print gene
			genes.add(gene)
			try:
				gene_files[gene].append(matrix_file)
			except:
				gene_files[gene] = [indir + "/" + matrix_file]
		np.save(outdir + "/unique_genes.npy", list(genes))
		np.save(outdir + "/unique_genes_matrix_filenames.npy", gene_files)

	print "Processing each gene..."
	for gene in genes:
		files = gene_files[gene]
		process_gene(gene, files, outdir)
	print "Done..."
		

def process_gene(gene_name, gene_matrix_full_filenames, outdir):
	norm_chunks = {gene_name: []}
	matrices = {} 
	print "Reading in chunks for gene " + gene_name + " and computing norm."
	for matrix_file in gene_matrix_full_filenames:
		chunk_number = matrix_file.split("_")[-1]
		matrix = np.load(matrix_file)
		#I know this structure doesn't make sense. Following what was done before
		contact_probs = np.linalg.norm(matrix, axis=1)
		norm_chunks[gene_name].append([chunk_number, None, contact_probs])
		matrices[chunk_number] = matrix

	print "Saving normalized matrices for gene " + gene_name + "."
	np.save(outdir + "/" + gene_name + "_norm.npy", norm_chunks)

	print "Saving raw matrices for gene " + gene_name + "."
	np.save(outdir + "/" + gene_name + "_raw.npy", matrices)
	

