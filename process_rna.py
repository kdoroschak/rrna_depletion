import glob
from multiprocessing import Pool
import numpy as np
import os
import psutil
import re
import shutil
import subprocess
import sys
import time

# reference_genome_path = "/homes/gws/sdorkenw/rrna/data/ref_genomes/GRCh38.p3_genomic.fa"
# reference_genome_path = "/homes/gws/sdorkenw/reference_genome_38/genome_GRCh38.fa"
# reference_genome_path = "/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa"
reference_genome_path = "/homes/gws/sdorkenw/reference_genome_38/GRCh38_o.p3.genome.fa"
# reference_genome_path = "/homes/gws/sdorkenw/rrna/data/ref_genomes/m10_genome.fa"
reference_genome_gt_path = "/homes/gws/sdorkenw/reference_genome_38/genes_GRCh38.gtf"
# reference_genome_gt_path = "/homes/gws/sdorkenw/reference_genome_38/data/ref_genomes/m10_genes.gtf"
# reference_genome_gt_path = "/homes/gws/sdorkenw/rrna/data/ref_genomes/rrna_hg38.gtf"
# reference_genome_exon_gt_path = "/homes/gws/sdorkenw/reference_genome_38/genes_exons.gtf"


def create_sampled_accession_file(accession_folder, n_accessions_per_file=2):
    accessions = []
    accession_paths = glob.glob(accession_folder + "/*.txt")
    for path in accession_paths:
        if not "sample" in path:
            with open(path, "r") as f:
                lines = f.readlines()
                line_ids = []
                while len(line_ids) < n_accessions_per_file:
                    r = np.random.randint(0, len(lines))
                    if r not in line_ids:
                        line_ids.append(r)
                for line_id in line_ids:
                    acc = lines[line_id]
                    accessions.append(acc)

    with open(accession_folder + "sample.txt", "w") as f:
        # for accession in accessions:
        f.writelines(accessions)


def _download_from_accession_file(args):
    name = args[0]
    save_path = args[1]

    call = "/homes/gws/sdorkenw/sratoolkit.2.5.7-ubuntu64/bin/prefetch %s" % \
           name
    subprocess.call(call, shell=True)

    shutil.move("/homes/gws/sdorkenw/ncbi/public/sra/%s.sra" % name,
                save_path)


def download_from_accession_file(save_path, accession_file, n_processes=1):
    multi_params = []
    with open(accession_file, "r") as f:
        for line in f.readlines():
           multi_params.append([line.rstrip('\n'), save_path])

    if n_processes > 1:
        pool = Pool(n_processes)
        results = pool.map(_download_from_accession_file, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_download_from_accession_file, multi_params)


def _convert_sra_to_fastq_worker(args):
    path = args[0]

    print path
    dir_name = os.path.dirname(path)
    subprocess.call("/homes/gws/sdorkenw/sratoolkit.2.5.7-ubuntu64/bin/fastq-dump -I --split-files --gzip --outdir %s %s" %
                    (dir_name, path), shell=True)


def convert_sra_to_fastq(paths, n_processes=1):
    multi_params = []
    for path in paths:
        multi_params.append([path])

    if n_processes > 1:
        pool = Pool(n_processes)
        results = pool.map(_convert_sra_to_fastq_worker, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_convert_sra_to_fastq_worker, multi_params)


def _run_rsubread_worker(args):
    paths = args[0]
    n_threads = args[1]
    n_gb = args[2]

    dir_name = os.path.dirname(paths[0]) + "/"
    base_name = os.path.basename(paths[0])

    prefix = ""
    for string in re.findall("[\w]+", base_name)[::-1]:
        if not string in ["fastq", "gz", "bam", "fq"]:
            prefix = string[:-2]
            break

    call = "Rscript --vanilla /homes/gws/sdorkenw/src/rna/run_rsubread_on_fastq.R %s %s %s %s %s %s %d %d" % \
                    (prefix, dir_name, reference_genome_path, paths[0], paths[1],
                     reference_genome_gt_path, n_gb*1000, n_threads)
    print call
    subprocess.call(call, shell=True)


def run_rsubread(fastq_files, n_processes=2, n_gb=10, n_threads=2):
    if n_processes % n_threads > 0:
        n_processes += n_threads - n_processes % n_threads

    mem = (psutil.virtual_memory().available)/1024**3
    if n_processes / 2. * n_gb > mem*.9:
        raise Exception("RAM usage too high for requested task: %f.1 GB / %f.1 GB" %
                        (n_processes / 2. * n_gb, mem))

    paths = []
    for path_1 in fastq_files:
        if "_1.fastq" in path_1:
            this_id = int(re.findall("[\d]+", path_1)[-2])
            this_tuple = ()
            for path_2 in fastq_files:
                if "_2.fastq" in path_2:
                    if this_id == int(re.findall("[\d]+", path_2)[-2]):
                        this_tuple = (path_1, path_2)
                        break

            if len(this_tuple) == 0:
                this_tuple = (path_1, "NULL")

            paths.append(this_tuple)


    multi_params = []
    for i in range(0, len(paths)):
        multi_params.append([paths[i], n_threads, n_gb])

    if n_processes > 2:
        pool = Pool(n_processes / 2)
        results = pool.map(_run_rsubread_worker, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_run_rsubread_worker, multi_params)


def _run_featureCounts_worker(args):
    path = args[0]
    n_threads = args[1]
    n_gb = args[2]
    paired_end = args[3]

    dir_name = os.path.dirname(path) + "/"
    base_name = os.path.basename(path)

    prefix = ""
    for string in re.findall("[\w]+", base_name):
        if not string in ["fastq", "gz", "bam", "fq"]:
            prefix = string
            break

    call = "Rscript --vanilla /homes/gws/sdorkenw/src/rna/run_featureCounts.R %s %s %s %s %s %s %d %d" % \
                    (prefix, dir_name, path, reference_genome_path, paired_end,
                     reference_genome_gt_path, n_gb*1000, n_threads)
    print call
    subprocess.call(call, shell=True)


def run_featureCounts(bam_files, paired_end="TRUE",
                      n_processes=2, n_gb=10, n_threads=2):
    if n_processes % n_threads > 0:
        n_processes += n_threads - n_processes % n_threads

    # mem = psutil.virtual_memory().free/1024**3
    # if n_processes / 2. * n_gb > mem*.7:
    #     raise Exception("RAM usage too high for requested task")

    multi_params = []
    for i in range(0, len(bam_files)):
        multi_params.append([bam_files[i], n_threads, n_gb, paired_end])

    if n_processes > n_threads:
        pool = Pool(n_processes / n_threads)
        results = pool.map(_run_featureCounts_worker, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_run_featureCounts_worker, multi_params)


def _sort_and_index_bams_worker(args):
    path = args[0]
    suffix = args[1]

    print "Start"
    call = "Rscript --vanilla /homes/gws/sdorkenw/src/rna/sort_index_bam.R %s %s" % \
            (path, path[:-4] + suffix)
    print call
    subprocess.call(call, shell=True)


def sort_and_index_bams(bam_files, out_suffix="_s", n_processes=2):
    multi_params = []
    for i in range(0, len(bam_files)):
        multi_params.append([bam_files[i], out_suffix])

    if n_processes > 1:
        pool = Pool(n_processes)
        results = pool.map(_sort_and_index_bams_worker, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_sort_and_index_bams_worker, multi_params)


def _run_seqbias_worker(args):
    path = args[0]
    n_threads = args[1]
    n_gb = args[2]
    paired_end = args[3]

    dir_name = os.path.dirname(path) + "/"
    base_name = os.path.basename(path)

    prefix = ""
    for string in re.findall("[\w]+", base_name):
        if not string in ["fastq", "gz", "bam", "fq"]:
            prefix = string
            break

    call = "Rscript --vanilla /homes/gws/sdorkenw/src/rna/run_seqbias.R %s %s %s %s %s %s %d %d" % \
                    (prefix, dir_name, path, reference_genome_path, paired_end,
                     reference_genome_exon_gt_path, n_gb*1000, n_threads)
    print call
    subprocess.call(call, shell=True)


def run_seqbias(bam_files, paired_end="TRUE",
                n_processes=1, n_gb=10, n_threads=4):
    if n_processes % n_threads > 0:
        n_processes += n_threads - n_processes % n_threads

    multi_params = []
    for i in range(0, len(bam_files)):
        multi_params.append([bam_files[i], n_threads, n_gb, paired_end])

    if n_processes > n_threads:
        pool = Pool(n_processes / n_threads)
        results = pool.map(_run_seqbias_worker, multi_params)
        pool.close()
        pool.join()
    else:
        results = map(_run_seqbias_worker, multi_params)


def main(input_dir, accession_file=None, count_only=False,
         n_processes=1):
    if accession_file is not None:
        download_from_accession_file(input_dir, accession_file, n_processes=n_processes)

    sra_files = glob.glob(input_dir + "/*.sra")

    # Convert sra files to fastq.gz files
    convert_sra_to_fastq(sra_files, n_processes=n_processes)

    fastq_files = glob.glob(input_dir + "/*.fastq.gz")

    if count_only:
        bam_files = glob.glob(input_dir + "/*/*out.bam")
        run_featureCounts(bam_files, paired_end="TRUE",
                          n_processes=n_processes)
    else:
        fastq_files = glob.glob(input_dir + "/*.fastq.gz")
        run_rsubread(fastq_files, n_processes)

    bam_files = glob.glob(input_dir + "/*.bam")

    sort_and_index_bams(bam_files, n_processes=n_processes)
    # sorted_bam_files = glob.glob(input_dir + "/*/*_s.bam")
    # time_s = time.time()
    # run_seqbias(sorted_bam_files, n_processes=n_processes, n_threads=1)
    # print "Time: %.3fs" % (time.time() - time_s)

if __name__ == "__main__":
    if not len(sys.argv) >= 3:
        print "Usage: python2 process_rna.py <input_dir> <n_processes>"
        raise Exception("Wrong arguments")

    if len(sys.argv) == 4:
        accession_file = sys.argv[3]
    else:
        accession_file = None

    main(sys.argv[1], accession_file, count_only=False,
         n_processes=int(sys.argv[2]))
