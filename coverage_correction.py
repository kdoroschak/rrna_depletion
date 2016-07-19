import numpy as np
import os
import sys
import pysam 			# pip install pysam
import gffutils 		# pip install gffutils


if sys.platform == "darwin":
  base = "/Volumes/cycle/"
else:
  base = "/homes/gws/" + os.getlogin() + "/rrna/" # make portable for katies and svens

source_home = base + "rrna/src/"
data_home = base + "rrna/data/"

# TODO
# Decide if we want this to output FASTQ or optionally run an aligner once it's all done
# Parallelize compute_expected_coverage and replete_coverage (must be run in serial but multiple instances)
# 

class CoverageCorrector(object):
	def __init__(self, depleted_genes=None, depleted_genes_file=None, 
					   gtf_db_file=None, gtf_file=None, 
					   bam_file=None, multiple_bam_file=None, 
					   out_fastq=True, out_bam=False):
		# depleted_genes: npy file containing a dict from a gene to region(s)
		# gtf_file: gtf containing at least the genes contained in depleted_genes
		# bam_file: bam containing read information for depleted_genes

		# In the init, we need to:
		#  - load the depleted genes from the file
		#  - either load or create the gtf database
		#  - make sure the bam file exists
		#  - 

		# Load depleted genes directly from a dict or from file
		self.depleted_genes = depleted_genes
		if self.depleted_genes is None:
			if depleted_genes_file is None:
				raise ValueError("Please supply a dictionary of depleted genes or a .npy file containing one.")
			else:
				self.depleted_genes = self.load_depleted_genes(depleted_genes_file)
		else:
			if not isinstance(depleted_genes, dict):
				raise TypeError("Loading the depleted genes file worked, but a dict was not found.")

		# Load gffutils database or generate it from a gtf file
		if gtf_db_file is None:
			if gtf_file is None:
				raise ValueError("Please supply a gffutils gtf db file or a gtf file to generate one.")
			else: # if db not already generated and gtf_file is given:
				self.gtf_db = self.generate_gtf_db(gtf_file)
		else: # if gtf_db_file is given:
			self.gtf_db = self.load_gtf_db(gtf_db_file)

		# Load a single or multiple bam files
		if bam_file is None:
			if multiple_bam_file is None:
				raise ValueError("Please supply a single bam file or a file that contains a list of bam files.")
			else:
				list_of_bam_files = np.loadtxt(multiple_bam_file)
				self.bam_files = self.load_bam_files(list_of_bam_files)
		else:
			self.bam_files = self.load_bam_files([bam_file])

		# Decide how to handle output
		if out_fastq:
			pass
		if out_bam:
			raise NotImplementedError("Got the message that you wanted to save a bam, but that's not set up yet.")
			print 

		pass

	def load_depleted_genes(self, depleted_genes_file):
		try:
			depleted_genes = np.load(depleted_genes_file)
		except:
			raise IOError("Could not open depleted genes file.")
		if not isinstance(depleted_genes, dict):
			raise TypeError("Loading the depleted genes file worked, but a dict was not found.")
		return depleted_genes

	def load_gtf_db(self, gtf_db_file):
		gtf_db = None
		if os.path.exists(gtf_db_file):
			gtf_db = gffutils.FeatureDB(gtf_db_file)
		return gtf_db

	def generate_gtf_db(self, gtf_file):
		gtf_db = None
		if os.path.exists(gtf_file):
			gtf_ext = gtf_file.split(".")[-1]
			gtf_db_file = ".".join(gtf_file.split(".")[:-1]) + ".db"
			print gtf_file, gtf_ext, gtf_db_file

			gtf_db = gffutils.create_db(gtf_file,
                                        dbfn=gtf_db_file,
                                        disable_infer_transcripts=True,
                                        disable_infer_genes=True,
                                        keep_order=True)
		return gtf_db

	def load_bam_files(self, list_of_bam_files):
		open_bam_files = []
		for bam_file in list_of_bam_files:
			bam = None
			try:
				bam = pysam.AlignmentFile(bam_file, "rb")
			except ValueError:
				print "File " + bam_file + " could not be read by pysam, and it will be ignored."
			if bam is not None:
				open_bam_files.append(bam)

		if len(open_bam_files) == 0:
			raise ValueError("Please provide at least one valid bam file.")

		return open_bam_files

	# def check_bam_file(self, bam_file):
	# 	try:
	# 		p = pysam.AlignmentFile(bam_file, "rb")
	# 		return True
	# 	except ValueError:
	# 		print "File " + bam_file + " could not be read by pysam."
	# 		return False


	def something_that_runs_this_whole_thing(self):
		pass

	def compute_expected_coverage(self, gene_name, depleted_region, gene_coverage):
		pass

	def replete_coverage(self):
		pass
	# 