import numpy as np
import pysam
import gffutils
import os

def get_coverage(bam, chromosome, start_pos, end_pos):
	coverage_by_position = []
	chromosome="chr"+chromosome
	# try:
	columns = bam.pileup(reference=chromosome, start=start_pos, end=end_pos)
	# except ValueError:
		# throw(ValueError)
		# return None
	for col in columns:
		coverage_by_position.append([col.pos, col.n])
	coverage_by_position = np.array(coverage_by_position)
	return coverage_by_position

# This function is here because these files shouldn't be reopened with every 
# call to compare_depleted_coverage but it's annoying to have to reload them
# in and remember the functions and proper files to do it. Will be unnecessary
# when this is part of a larger script, but for testing this is ok
def prep(bam_file="/projects/bio/rrna/data/rsubread_aligned_v2/SRR891244_s.bam", 
	     gtf_db_file="/projects/bio/rrna/data/annotations/Homo_sapiens.GRCh38.84.features_name.db"):
	# bam, gtf_db = cc.prep()
	bam = pysam.AlignmentFile(bam_file, "rb")
	gtf_db = gffutils.FeatureDB(gtf_db_file)
	return bam, gtf_db

def compare_depleted_coverage(gene_name, chunk_numbers, gtf_db, bam, chunk_length, exons_only=True): #exclude_zero_cov=False or exclude_cov_le=2 aka less than or equal to
	#cc.compare_depleted_coverage("PSMA2P2", [1,2], gtf_db, bam, 50, exons_only=False)
	gene = gtf_db[gene_name]
	gene_start = gene.start
	gene_end = gene.end
	gene_chrom = gene.chrom
	
	contiguous_chunks = determine_contiguous_chunks(chunk_numbers)
	avg_coverage = np.zeros(2) # (i, 0) = coverage in region; (i, 1) = coverage not in region
	coverage_in_region = []
	coverage_not_in_region = []

	if exons_only:
		# Make a list of ranges so we can check each exon for overlap
		depleted_regions = []
		for contig_chunk in contiguous_chunks:
			depleted_regions.append(range(calculate_pos_in_chr(gene_start, contig_chunk[0], contig_chunk[1], chunk_length)))

		# Tally up the coverage for exons in and not in the depleted regions
		# eventually surround this exon generator with try/catch but i want it to throw the error for now
		exon_generator = gtf_db.children(gene_name, featuretype="exon", order_by="start")
		for exon in exon_generator:
			print exon
			start_pos = exon.start
			end_pos = exon.end
			exon_range = range(start_pos, end_pos)
			in_depleted_region = check_overlap(exon_range, depleted_regions)
			coverage = get_coverage(bam, gene_chrom, start_pos, end_pos)
			if in_depleted_region:
				coverage_in_region.extend(list(coverage))
			else:
				coverage_not_in_region.extend(list(coverage))	

	else:
		for contig_i, contig_chunk in enumerate(contiguous_chunks):
			start_pos, end_pos = calculate_pos_in_chr(gene_start, contig_chunk[0], contig_chunk[1], chunk_length)
			coverage_in_region.extend(get_coverage(bam, gene_chrom, start_pos, end_pos))
			coverage_not_in_region.extend(list(get_coverage(bam, gene_chrom, gene_start, start_pos)))
			coverage_not_in_region.extend(list(get_coverage(bam, gene_chrom, end_pos, gene_end)))
	
	avg_coverage[0] = np.mean(np.array(coverage_in_region)[:,1])
	avg_coverage[1] = np.mean(np.array(coverage_not_in_region)[:,1])
	
	return avg_coverage

def check_overlap(range_to_test, target_ranges):
	range_to_test = set(range_to_test)
	for target_range in target_ranges:
		if range_to_test.intersection(target_range):
			return True
	return False

def determine_contiguous_chunks(chunk_numbers):
	chunk_numbers.sort()
	len_chunk_numbers = len(chunk_numbers)
	contig_chunks = []
	contig_chunk_start = chunk_numbers[0]
	last_chunk = chunk_numbers[0]
	
	if len_chunk_numbers == 0:
		return None
	elif len_chunk_numbers == 1:
		contig_chunks = [(chunk_numbers[0], chunk_numbers[0])]
	elif len_chunk_numbers == 2:
		if chunk_numbers[0] + 1 == chunk_numbers[1]:
			contig_chunks = [(chunk_numbers[0], chunk_numbers[1])]
		else:
			contig_chunks = [(chunk_numbers[0], chunk_numbers[0]), (chunk_numbers[1], chunk_numbers[1])]
	else:
		curr_chunk = chunk_numbers[1]
		pos = 2
		while pos < len_chunk_numbers:
			if last_chunk + 1 == curr_chunk: 	# continue contiguous region
				last_chunk = curr_chunk
				curr_chunk = chunk_numbers[pos]
				pos += 1
			else:								# start new region
				contig_chunks.append((contig_chunk_start, last_chunk))
				contig_chunk_start = curr_chunk
				last_chunk = curr_chunk
				curr_chunk = chunk_numbers[pos]
				pos += 1
		# terminate - figure out whether to create the last one on its own or add to the last region
		if last_chunk + 1 == curr_chunk:
			contig_chunks.append((contig_chunk_start, curr_chunk))
		else:
			contig_chunks.append((contig_chunk_start, last_chunk))
			contig_chunks.append((curr_chunk, curr_chunk))
			
	return contig_chunks

def calculate_pos_in_chr(gene_start, start_chunk, end_chunk, chunk_length):
	start_pos_within_gene = start_chunk * chunk_length
	end_pos_within_gene = (end_chunk + 1) * chunk_length
	start_pos = gene_start + start_pos_within_gene
	end_pos = gene_start + end_pos_within_gene
	return start_pos, end_pos
