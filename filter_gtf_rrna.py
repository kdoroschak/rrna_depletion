import numpy as np

def separate_gtf_rrna_protein(gtf_file="/homes/gws/kdorosch/rrna/data/annotations/Homo_sapiens.GRCh38.84.gtf", 
							  out_path="/homes/gws/kdorosch/rrna/data/annotations/gtf_subsets/"):
	gtf_fname = gtf_file.split("/")[-1][:-4]

	header_lines = []
	rrna_lines = []
	mt_rrna_lines = []
	ribo_protein_lines = []
	other_lines = []

	with open(gtf_file, "r") as f:
		for line_i,line in enumerate(f):
			# line = line.strip()
			if line.startswith("#!"):
				header_lines.append(line)
			else:
				if "rRNA" in line or "rrna" in line:
					if "MT-" in line or "Mt_" in line:
						mt_rrna_lines.append(line)
					else:
						rrna_lines.append(line)
				elif "RPS" in line or "RPL" in line:
					ribo_protein_lines.append(line)
				else:
					other_lines.append(line)


	# Create new GTF files
	#   mt_rrna and rrna
	with open(out_path + gtf_fname + ".rrna_mtrrna.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes only rRNA and mt_rRNA.\n")		
		f.writelines(rrna_lines)
		f.writelines(mt_rrna_lines)

	#   rrna
	with open(out_path + gtf_fname + ".rrna.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes only rRNA.\n")
		f.writelines(rrna_lines)
	
	# 	mt_rrna
	with open(out_path + gtf_fname + ".mtrrna.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes only mt_rRNA.\n")
		f.writelines(mt_rrna_lines)

	#	ribosomal protein
	with open(out_path + gtf_fname + ".ribo_proteins.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes only ribosomal protein genes.\n")
		f.writelines(ribo_protein_lines)

	#	original without rrna
	with open(out_path + gtf_fname + ".no_rrna.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes everything but strict rRNA (does have mt_rrna).\n")
		f.writelines(other_lines)
		f.writelines(mt_rrna_lines)
		f.writelines(ribo_protein_lines)

	#	original without rrna or mt_rrna
	with open(out_path + gtf_fname + ".no_rrna_or_mtrrna.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes everything but rRNA and mt_rrna.\n")
		f.writelines(other_lines)
		f.writelines(ribo_protein_lines)

	# 	original without ribosomal protein
	with open(out_path + gtf_fname + ".no_ribo_proteins.gtf", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Includes everything but ribosomal protein genes.\n")
		f.writelines(other_lines)
		f.writelines(mt_rrna_lines)
		f.writelines(rrna_lines)

	# 	original without either rrna or ribosomal protein
	with open(out_path + gtf_fname + ".exclude_all_ribo", "w") as f:
		f.writelines(header_lines)
		f.write("#! Note from KJD: Excludes anything having to do with the ribosome.\n")
		f.writelines(other_lines)





