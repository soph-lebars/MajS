import argparse
import json
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm


def cli_parser() -> argparse.ArgumentParser:
    """Define the CLI parser

    Returns:
        argparse.ArgumentParser: Argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-path', '-d', type=str, required=True,
                        help='Path to the folder contains experimentd data (str)')
    parser.add_argument('--out-path', '-o', type=str, required=True,
                        help='Output path (str')
    parser.add_argument('--benchmark', '-b', type=int, required=True,
                        help='Benchmark number (1 or 2)')
    return parser


def parse_args(args: iter = None) -> dict:
    """Parse arguments

    Args:
        args (iter, optional): Arguments to parse. Defaults to None.

    Returns:
        dict: Arguments
    """
    return cli_parser().parse_args(args)


def get_scores(hc: float, lc: float, path: str, benchmark: str) -> dict:
    """Get scores of genes for a specific hc and lc

    Args:
        hc (float): High confidence value
        path (str): Path to folder contains all experiments data
        benchmark (str): Banchmark number

    Returns:
        dict: Current scores of genes for epsilon
    """
    with open(f'{path}/hc_{hc}_lc_{lc}/{benchmark}/gene_scores.json', 'r') as file_:
        current_scores = json.load(file_)
    for gene in current_scores.keys():
        for to_del in ['Sign-1', 'Sign1', 'Sign0', 'SignMaj', 'iggy', 'logFC']:
            del current_scores[gene][to_del]
    return current_scores


def compare_scores(majs_score: float, iggy_score: float) -> str:
    """Compare scores value of both method

    Args:
        majs_score (float): MajS score
        iggy_score (float): Iggy score

    Returns:
        str: Comparison of the two scores
    """
    if majs_score > iggy_score:
        return 'majs'
    elif iggy_score > majs_score:
        return 'iggy'
    else:
        return 'equal'


def get_genes_list(path: str, benchmark: str) -> list:
    """Get gene list from a folder

    Args:
        path (str): Path of the folder
        benchmark (str): Benchmark number

    Returns:
        list: List of genes present in the folder for the specific benchmark
    """
    expe_list = os.listdir(path)
    genes_list = [f for f in os.listdir(
        f'{path}/{expe_list[0]}/{benchmark}') if f.endswith('.pdf')]
    genes_list.remove(f'score_plot_{benchmark}.pdf')
    genes_list = [x.split('.')[0] for x in genes_list]
    return genes_list


if __name__ == '__main__':
    # Parse arguments
    args = parse_args()

    data_path = args.data_path
    out_path = args.out_path
    benchmark = f'benchmark{args.benchmark}'

    out_path = f'{out_path}/{benchmark}'

    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    if not os.path.isdir(f'{out_path}/surface_plots'):
        os.makedirs(f'{out_path}/surface_plots')

    folders = os.listdir(data_path)

    # Get statistics of the comparison of both methods
    results = dict()
    for expe in folders:
        _, hc, _, lc = expe.split('_')
        with open(f'{data_path}/{expe}/{benchmark}/out.log', 'r') as f:
            lines = f.readlines()
            idx = lines.index(
                '### Statistics of the comparison of MajS and Iggy ###\n')
            idx += 1  # Get the index of stats
            stats = lines[idx]
            stats = stats.replace("\'", "\"")
            stats = json.loads(stats)
            to_keep = ['majs>iggy', 'iggy>majs', 'equal']
            reduced_stats = {key: stats[key] for key in to_keep}
            if hc in results.keys():
                results[hc][lc] = reduced_stats
            else:
                results[hc] = {lc: reduced_stats}

    # Export statistics
    json.dump(results, open(
        f'{out_path}/methods_comparison.json', 'w'), sort_keys=True)

    # Compute surface plot
    hc_values = ['0.01', '0.02', '0.03', '0.04',
                 '0.05', '0.06', '0.07', '0.08', '0.09', '0.10']
    lc_values = ['0.1', '0.2', '0.3', '0.4',
                 '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']

    fig, axes = plt.subplots(1, 3, subplot_kw=dict(projection='3d'),
                             figsize=plt.figaspect(1/3))
    for feature, ax in zip(['majs>iggy', 'iggy>majs', 'equal'], axes):
        # Make data
        X = np.arange(10)
        Y = np.arange(10)
        X, Y = np.meshgrid(X, Y)
        Z = list()
        for i in range(10):
            Z.append([])
            for j in range(10):
                Z[i].append(results[hc_values[i]][lc_values[j]][feature])
        Z = np.array(Z)

        # Plot the surface
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        # Legend configuration
        ax.set_xticklabels(hc_values)
        ax.set_yticklabels(lc_values)
        ax.view_init(40, 20)
        ax.set_xlabel(r'$\sigma_{hc}$')
        ax.set_ylabel(r'$\sigma_{lc}$')
        ax.set_zlabel('# of genes')
        ax.locator_params(axis="x", nbins=10)
        ax.locator_params(axis="y", nbins=10)
        if feature == 'majs>iggy':
            subtitle = 'score(MajS) > score(Iggy)'
        elif feature == 'iggy>majs':
            subtitle = 'score(Iggy) > score(MajS)'
        else:
            subtitle = 'score(MajS) = score(Iggy)'
        ax.set_title(subtitle)
    plt.savefig(f'{out_path}/surface_plots/plot.pdf')
    plt.close()
