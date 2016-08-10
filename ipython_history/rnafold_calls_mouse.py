# coding: utf-8
get_ipython().magic(u'cd ~/rrna/src')
import featureExtraction as featex
import generate_rnafold_features as rnafold

# rRNA
fe = featex.FeatureExtractor(fasta_file_path="/projects/bio/rrna/data/ref_genomes/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa",
					  gtf_file_path="/projects/bio/rrna/data/annotations/mouse/gtf_subsets/Mus_musculus.GRCm38.85.rrna.gtf",
					  save_path="/projects/bio/rrna/data/featureExtractors/mouse/Mus_musculus.GRCm38.85_gene_features.pkl",
					  chunk_save_path="/projects/bio/rrna/data/featureExtractors/mouse/Mus_musculus.GRCm38.85_gene_chunks.pkl")
fe.generate_gene_features()
fe.generate_chunked_sequences()
results = rnafold.predict_secondary_for_chunks_parallel(fe, "rnafold_mouse_rrna", gene_list="/projects/bio/rrna/data/annotations/mouse/gtf_subsets/Mus_musculus.GRCm38.85.rrna.names.npy", n_processes=40)

# all but rRNA
fe = featex.FeatureExtractor(fasta_file_path="/projects/bio/rrna/data/ref_genomes/mouse/Mus_musculus.GRCm38.dna.primary_assembly.fa",
					  gtf_file_path="/projects/bio/rrna/data/annotations/mouse/gtf_subsets/Mus_musculus.GRCm38.85.no_rrna.gtf",
					  save_path="/projects/bio/rrna/data/featureExtractors/mouse/Mus_musculus.GRCm38.85_gene_features.pkl",
					  chunk_save_path="/projects/bio/rrna/data/featureExtractors/mouse/Mus_musculus.GRCm38.85_gene_chunks.pkl")
fe.generate_gene_features()
fe.generate_chunked_sequences()
results = rnafold.predict_secondary_for_chunks_parallel(fe, "rnafold_mouse_no_rrna", gene_list="/projects/bio/rrna/data/annotations/mouse/gtf_subsets/Mus_musculus.GRCm38.85.no_rrna.names.npy", n_processes=40)