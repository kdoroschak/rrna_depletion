import glob
from matplotlib import pyplot as plt
import numpy as np
import re


def mean_std(data, names, save_path, xlabel=None, ylabel=None, y_max=20):
    # order = np.argsort(names)
    # data = np.array(data)[order]
    # names = np.array(names)[order]

    # data = [data[i][data[i] > 0] for i in range(len(data))]

    lengths = [len(data[i]) for i in range(len(data))]
    print lengths
    print [np.mean(data[i]) for i in range(len(data))]
    print [np.median(data[i]) for i in range(len(data))]

    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')

    y = [np.mean(data[i]) for i in range(len(data))]
    yerr = [np.std(data[i]) / (lengths[i])**(1./2) for i in range(len(lengths))]

    ax.errorbar(range(1, 1+len(data)), y, yerr=yerr, fmt='s', color=".3", capthick=2,
                lw=5, ms=7, mec="None", ecolor="k", elinewidth=3)

    ax.tick_params(axis='x', which='major', labelsize=16, direction='out',
                   length=0, width=0, right="off", top="off", pad=10)

    ax.tick_params(axis='y', which='major', labelsize=20, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.tick_params(axis='x', which='minor', labelsize=12, direction='out',
                   length=0, width=0, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.ylim(0, y_max)
    # ax.set_yscale("log")
    # ax.set_yticks(np.arange(0, 1.55, .5))
    # ax.set_yticks(np.arange(0, 1.55, .25), minor=True)
    # ax.set_yticks(np.arange(0, 1.1, .1), minor=True)
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=20)

    plt.xticks(np.arange(1, len(names) + 1), names, rotation=90)
    ax.set_xlim(0, len(data)+2)
    # legend = plt.legend(loc=loc, frameon=False, prop={'size': 23})
    # legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    # plt.show(block=True)
    plt.close()



def boxplot(data, names, save_path, xlabel=None, ylabel=None, y_max=20):
    # data = [data[i][data[i] > 0] for i in range(len(data))]
    print [len(data[i]) for i in range(len(data))]
    print [np.mean(data[i]) for i in range(len(data))]
    print [np.median(data[i]) for i in range(len(data))]

    fig, ax = plt.subplots()
    fig.patch.set_facecolor('white')

    bp = plt.boxplot(data,  patch_artist=True, notch=0, sym='-', vert=1, whis=1.5, showmeans=True)

    colors = [".7"] * len(data)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    for box, color in zip(bp['boxes'], colors):
        # change outline color
        box.set(color='.7', linewidth=2)
        # change fill color
        box.set(facecolor="w")

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='0.3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='0.3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='0.', linewidth=2)

    # for mean in bp['means']:
    #     mean.set(color='0.', linewidth=2)

    ax.tick_params(axis='x', which='major', labelsize=16, direction='out',
                   length=0, width=0, right="off", top="off", pad=10)

    ax.tick_params(axis='y', which='major', labelsize=20, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.tick_params(axis='x', which='minor', labelsize=12, direction='out',
                   length=0, width=0, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.ylim(0, y_max)
    # ax.set_yscale("log")
    # ax.set_yticks(np.arange(0, 1.55, .5))
    # ax.set_yticks(np.arange(0, 1.55, .25), minor=True)
    # ax.set_yticks(np.arange(0, 1.1, .1), minor=True)
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=20)

    plt.xticks(np.arange(1, len(names) + 1), names, rotation=30)
    # legend = plt.legend(loc=loc, frameon=False, prop={'size': 23})
    # legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    # plt.show(block=True)
    plt.close()


def histogram(data, names, save_path, n_bars=40, xlabel=None, ylabel=None):
    data = [data[i][data[i] > 0] for i in range(len(data))]
    fig, ax = plt.subplots()

    fig.patch.set_facecolor('white')

    ax.tick_params(axis='x', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.tick_params(axis='x', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=18)
    # ax.set_yticks(np.arange(0, 0.5, 0.1))
    # ax.set_xticks(np.arange(0, 1.1, 0.5))
    # ax.set_xticks(np.arange(0, 1.1, 0.25), minor=True)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=18)

    weights0 = np.ones_like(data[0]) / len(data[0])
    weights1 = np.ones_like(data[1]) / len(data[1])
    # raise()
    ax.hist(data[0], bins=np.arange(n_bars) / float(n_bars - 1)*20, alpha=0.7,
            label=names[0], weights=weights0, lw=0,
            color=[39 / 255., 184 / 255., 148 / 255.], histtype="stepfilled")

    ax.hist(data[1], bins=np.arange(n_bars) / float(n_bars - 1)*20, alpha=0.7,
            color=[239 / 255., 65 / 255., 50 / 255.], label=names[1], lw=0,
            histtype="stepfilled", weights=weights1)

    legend = ax.legend(loc="upper right", frameon=False, prop={'size': 18})

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
 

def plot_all_basepair(k, analysis_folder="/homes/gws/sdorkenw/rrna/data/kmer_analysis/",
                      save_folder="/homes/gws/sdorkenw/rrna/", y_max=100,
                      cleaned=True, contiguous=False, paths=None):
    coverage_dict = {}
    if paths is None or not isinstance(paths, list) or len(paths) == 0:
        paths = glob.glob(analysis_folder + "/*k%d*" % k)

    for path in paths:
        if not "/read_" in path and (cleaned == bool("cleaned" in path)):
            print path
            if contiguous and "cont." in path:
                seq_id = re.findall("[\d]+", path)[-1]
                if not seq_id in coverage_dict:
                    coverage_dict[seq_id] = {}

                if "nonribo" in path:
                    coverage_dict[seq_id][0] = this_data
                else:
                    coverage_dict[seq_id][1] = this_data
            elif not contiguous and not "cont." in path:
                seq_id = re.findall("[\d]+", path)[-1]

                this_data = np.load(path)
                if isinstance(this_data[0], list):
                    this_data = [item for sublist in this_data for item in sublist]

                if not seq_id in coverage_dict:
                    coverage_dict[seq_id] = {}
                if "nonribo" in path:
                    coverage_dict[seq_id][0] = this_data
                else:
                    coverage_dict[seq_id][1] = this_data

    data = []
    names = []
    for key in coverage_dict.keys():
        data.append(np.array(coverage_dict[key][0]))
        data.append(np.array(coverage_dict[key][1]))
        print "nonribo", np.median(coverage_dict[key][0]), np.mean(coverage_dict[key][0])
        print "ribo", np.median(coverage_dict[key][1]), np.mean(coverage_dict[key][1])
        names.append("nonribo %s" % key)
        names.append("ribo %s" % key)

        # if contiguous:
        #     histogram(data[-2:], names[-2:], save_folder+"/hist_%d_%s_cont.png" % (k, key),
        #               xlabel="kmer av expresion", ylabel="normed ratio")
        # else:
        #     histogram(data[-2:], names[-2:], save_folder+"/hist_%d_%s.png" % (k, key),
        #               xlabel="kmer av expresion", ylabel="normed ratio")

    if contiguous:
        mean_std(data, names, save_path=save_folder + "mean_std_summary_%d_cont.png" % k,
                ylabel="kmer average expression", y_max=y_max)
    else:
        mean_std(data, names, save_path=save_folder + "mean_std_summary_%d.png" % k,
                ylabel="kmer average expression", y_max=y_max)


def plot_all_reads(k, analysis_folder="/homes/gws/sdorkenw/rrna/data/kmer_analysis/",
                   save_folder="/homes/gws/sdorkenw/rrna/",
                   y_max=100, start_only=True, paths=None):
    coverage_dict = {}
    if paths is None or not isinstance(paths, list) or len(paths) == 0:
        paths = glob.glob(analysis_folder + "/*k%d*" % k)

    for path in paths:
        if "/read_" in path:

            this_data = np.load(path).item().values()

            if start_only and "start_only" in path:
                print path
                seq_id = re.findall("[\d]+", path)[-1]
                if not seq_id in coverage_dict:
                    coverage_dict[seq_id] = {}

                if "nonribo" in path:
                    coverage_dict[seq_id][0] = this_data
                else:
                    coverage_dict[seq_id][1] = this_data

            elif not start_only and not "start_only" in path:
                seq_id = re.findall("[\d]+", path)[-1]

                this_data = np.load(path)
                if isinstance(this_data[0], list):
                    this_data = [item for sublist in this_data for item in sublist]

                if not seq_id in coverage_dict:
                    coverage_dict[seq_id] = {}
                if "nonribo" in path:
                    coverage_dict[seq_id][0] = this_data
                else:
                    coverage_dict[seq_id][1] = this_data

    data = []
    names = []
    for key in coverage_dict.keys():
        data.append(coverage_dict[key][0])
        data.append(coverage_dict[key][1])
        print np.median(coverage_dict[key][0]), np.mean(coverage_dict[key][0])
        print np.median(coverage_dict[key][1]), np.mean(coverage_dict[key][1])
        names.append("nonribo %s" % key)
        names.append("ribo %s" % key)

        # if start_only:
        #     histogram(data[-2:], names[-2:],
        #               save_folder + "/hist_reads_%d_%s_start_only.png" % (k, key),
        #               xlabel="kmer av expresion", ylabel="normed ratio")
        # else:
        #     histogram(data[-2:], names[-2:],
        #               save_folder + "/hist_reads_%d_%s.png" % (k, key),
        #               xlabel="kmer av expresion", ylabel="normed ratio")

    if start_only:
        boxplot(data, names,
                save_path=save_folder + "box_plot_reads_summary_%d_start_only.png" % k,
                ylabel="kmer av expression", y_max=y_max)
    else:
        boxplot(data, names, save_path=save_folder + "box_plot_reads_summary_%d.png" % k,
                ylabel="kmer av expression", y_max=y_max)


def plot_relative_coverage(distances, rel_cov, nbins=40):
    fig, ax = plt.subplots()

    fig.patch.set_facecolor('white')

    ax.tick_params(axis='x', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.tick_params(axis='x', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    bins = np.linspace(0, 1, nbins)

    idx = np.digitize(distances, bins)
    running_median = [np.mean(rel_cov[idx == k]) for k in range(nbins)]

    plt.scatter(distances, rel_cov, color='k', alpha=.005, s=2)
    plt.plot(bins - 1 / 2, running_median, 'r--', lw=2, alpha=.8)

    ax.set_xlim(0, 0.45)
    ax.set_ylim(0.0, 2)
    # ax.set_yscale("log")


    plt.tight_layout()
    # plt.savefig(save_path, dpi=300)
    plt.show(block=True)
    plt.close()


def plot_pca(X, Y):
    fig, ax = plt.subplots()

    fig.patch.set_facecolor('white')

    ax.tick_params(axis='x', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='major', labelsize=22, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.tick_params(axis='x', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)
    ax.tick_params(axis='y', which='minor', labelsize=12, direction='out',
                   length=8, width=3, right="off", top="off", pad=10)

    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.scatter(X, Y, color='k', alpha=.2, s=2)
       
    plt.tight_layout()
    # plt.savefig(save_path, dpi=300)
    plt.show(block=True)
    plt.close()
