# from matplotlib import pyplot as plt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
import gffutils
import os
from multiprocessing import Pool

# TODO expand this to plot coverage from multiple bam files vertically

# covplot.plot_from_list_parallel("../data/first5k.txt", "../data/rsubread_aligned_v2/SRR891244_s.bam", "../data/annotations/Homo_sapiens.GRCh38.84.gtf", "../data/coverage_plots", n_processes=40)

def plot_from_list_parallel(gene_names_file, bam_file, gtf_file_path, save_path, n_processes=1):
	# Determine if we can load or need to regenerate the gtf database
	gtf_db_path = gtf_file_path[:-4] + ".features_name.db"
	if not os.path.exists(gtf_db_path):
		print "Error: gtf db " + gtf_db_path + " does not exist."
		return False

	if not os.path.exists(bam_file):
		print "Error: bam file " + bam_file + " cannot be found."
		return False

	if os.path.exists(gene_names_file):
		with open(gene_names_file, "r") as f:
			list_of_gene_names = f.readlines()
	list_of_gene_names = [name.strip() for name in list_of_gene_names]

	multi_params = []

	for gene_name in list_of_gene_names:
		multi_params.append([gene_name, bam_file, save_path, gtf_db_path])

	if n_processes > 1:
		pool = Pool(n_processes)
		# pool.map(_plot_helper, multi_params)
		# results = pool.map_async(slowly_square, range(40)).get(9999999)
		pool.map_async(_plot_helper, multi_params).get(9999999)
		print "Closing pool"
		pool.close()
		pool.join()
	else:
		map(_plot_helper, multi_params)

def _plot_helper(args):
	gene_name = args[0]
	bam_file = args[1]
	save_path = args[2]
	gtf_db_path = args[3]

	gtf_db = gffutils.FeatureDB(gtf_db_path)
	try:
		gene = gtf_db[gene_name]
	except:
		print "Gene " + gene_name + " is not in the database."
		return

	print "Getting coverage for " + gene_name
	coverage = get_coverage(bam_file, gene.chrom, gene.start, gene.end)

	if coverage is None:
		print "Gene " + gene_name + " has an invalid chromosome id: " + gene.chrom
		return
	elif coverage is not None and coverage.size == 0:
		# TODO generate empty plot? message?
		print "Gene " + gene_name + " has no coverage."
		return

	plot_file = save_path + "/" + gene_name + "_coverage.png"

	print "Plotting " + gene_name
	bam_name = bam_file.split("/")[-1]
	plot_coverage(coverage, plot_file, start_pos=gene.start, end_pos=gene.end, title=gene_name+" Coverage Plot " + bam_name)


def get_coverage(bam_file, chromosome, start_pos, end_pos):
	coverage_data = []
	bam = pysam.AlignmentFile(bam_file, "rb")
	chromosome="chr"+chromosome
	try:
		columns = bam.pileup(reference=chromosome, start=start_pos, end=end_pos)
	except ValueError:
		return None
	for col in columns:
		coverage_data.append([col.pos, col.n])
	coverage_data = np.array(coverage_data)
	return coverage_data

# Much of this is copied from plotUtils histogram function
# TODO add exon as lower subplot
def plot_coverage(coverage_data, save_path, start_pos=None, end_pos=None, xlabel="Position", ylabel="Coverage", title=None, max_bar_height=None, width=None):
	# Coverage data is Nx2 with columns representing position, coverage
    fig, ax = plt.subplots()

    fig.patch.set_facecolor('white')

    if start_pos is None:
    	start_pos = coverage_data[0,0]
    if end_pos is None:
    	end_pos = coverage_data[-1,0]

    # Resize the window so if we have a super long range, it will widen
    plot_span = end_pos - start_pos
    if width is not None:
    	if width > 4:
    		w = width
    		fig.set_size_inches([w,6])
    elif plot_span > 400:
    	w = plot_span/50
    	w = min(w,20)
    	fig.set_size_inches([w, 6])

    if max_bar_height is None:
    	max_bar_height = max(coverage_data[:,1])

    ax.tick_params(axis='x', which='major', labelsize=20, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=18)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=18)
    if title is not None:
    	ax.set_title(title, fontsize=20)

    ax.bar(coverage_data[:,0], coverage_data[:,1], color="black", edgecolor="black")
    # plt.ylim(ymax=max_bar_height)
    plt.ylim((0,max_bar_height))
    # plt.yl

    # legend = ax.legend(loc="upper right", frameon=False, prop={'size': 18})

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()