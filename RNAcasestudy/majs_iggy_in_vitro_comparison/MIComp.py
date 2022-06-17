import argparse
import json
import os
import matplotlib.pyplot as plt
import pandas as pd
from typing import Callable
import numpy as np
from math import exp, pi, sqrt, pow
from scipy.integrate import quad


def cli_parser() -> argparse.ArgumentParser:
    """Define the CLI parser

    Returns:
        argparse.ArgumentParser: Argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', '-o', type=str, required=True,
                        help='out directory name')
    parser.add_argument('--majs-file', '-mf', type=str, required=True,
                        help='majS data file (csv)')
    parser.add_argument('--iggy-file', '-if', type=str, required=True,
                        help='iggy data file (csv)')
    parser.add_argument('--rna-file', '-rf', type=str, required=True,
                        help='RNA data file (csv)')
    parser.add_argument('--obs-file', '-of', type=str, required=True,
                        help='observations data file (csv)')

    parser.add_argument('--high_confidence_coefficient', '-hc', type=float, required=False,
                        help='high confidence coefficient (float)')
    parser.add_argument('--low_confidence_coefficient', '-lc', type=float, required=False,
                        help='low confidence coefficient (float)')
    parser.add_argument('--epsilon', '-eps', type=float, required=False,
                        help='epsilon value (float)')
    parser.add_argument('--export', '-e', action='store_true',
                        help='export all predictions in a JSON file')
    return parser


def parse_args(args: iter = None) -> dict:
    """Parse arguments

    Args:
        args (iter, optional): Arguments to parse. Defaults to None.

    Returns:
        dict: Arguments
    """
    return cli_parser().parse_args(args)


def load_majs(filename: str, all_predictions: dict) -> dict:
    """Load MajS' prediction

    Args:
        filename (str): MajS' predictions filename
        all_predictions (dict): Gene predictons

    Returns:
        dict: Gene predictons
    """
    print("Load majS' file...")
    majs_data = pd.read_csv(filename)
    features = majs_data.columns
    for _, row in majs_data.iterrows():
        gene_name = row[0].upper().strip()
        if gene_name not in all_predictions.keys():
            all_predictions[gene_name] = {}
        for i in range(len(row[1:])):
            pred = row[i+1].upper().strip()
            if pred[0] == '[':
                pred = list(map(lambda x: int(x), pred.replace(
                    "[", "").replace("]", "").split(",")))
            else:
                pred = pred.split()
            all_predictions[gene_name][features[i+1]] = pred

    return all_predictions


def load_rna(filename: str, all_predictions: dict) -> dict:
    """Load RNA data

    Args:
        filename (str): RNA fold changes filename
        all_predictions (dict): Gene predictons

    Returns:
        dict: Gene predictons
    """
    print("Load RNA's file...")
    rna_data = pd.read_csv(filename)
    features = rna_data.columns
    for _, row in rna_data.iterrows():
        gene_name = row[0].upper().strip()
        if gene_name not in all_predictions.keys():
            all_predictions[gene_name] = {}
        all_predictions[gene_name][features[1]] = row[1]
    return all_predictions


def load_iggy(filename: str, all_predictions: dict) -> dict:
    """Load Iggy's prediction

    Args:
        filename (str): Iggy's predictions filename
        all_predictions (dict): Gene predictons

    Returns:
        dict: Gene predictons
    """
    print("Load Iggy's file...")
    iggy_data = pd.read_csv(filename, sep='=', header=None)
    for _, row in iggy_data.iterrows():
        gene_name = row[0].upper().strip()
        if gene_name not in all_predictions.keys():
            all_predictions[gene_name] = {}
        all_predictions[gene_name]['iggy'] = row[1].strip()
    return all_predictions


def load_observations(filename: str) -> list:
    """Get observed genes

    Args:
        filename (str): Observed gene filename

    Returns:
        list: Observed genes
    """
    print("Load observations' file...")
    obs_data = pd.read_csv(filename, sep='=', header=None)
    return [x.upper().strip() for x in obs_data.iloc[:, 0]]


def export_to_json(dict_: str, filename: str = 'out.json'):
    """Export dictionnary to JSON file

    Args:
        dict_ (str): Dictionnary to export
        filename (str, optional): Name of the output JSON file. Defaults to 'out.json'.
    """
    with open(filename, 'w') as outfile:
        json.dump(dict_, outfile, sort_keys=True, indent=4)


def normal(x: float, m: float, sd: float) -> float:
    """Compute a normal distribution

    Args:
        x (float): X value
        m (float): Mean
        sd (float): Standard deviation

    Returns:
        list: Normal value of X
    """
    return exp(-((x-m)/sd)**2/2)/(sd*sqrt(2*pi))


def fnormal(m: float, sd: float) -> Callable:
    """Get Normal distribution

    Args:
        m (float): Mean
        sd (float): Standard deviation

    Returns:
        Callable: Normal distribution
    """
    return lambda x: normal(x, m, sd)


def scale(l: float, f: Callable) -> Callable:
    """Scale the normal distribution into the mixture

    Args:
        l (float): Coefficient
        f (Callable): Normal distribution function

    Returns:
        Callable: Scaled normal distribution function
    """
    return lambda x: l*f(x)


def add(f1: Callable, f2: Callable) -> Callable:
    """Add a distribution to another

    Args:
        f1 (Callable): Normal distribution function
        f2 (Callable): Normal distribution function

    Returns:
        Callable: Added normal distribution function
    """
    return lambda x: f1(x)+f2(x)


def mixture(listdist: list) -> add:
    """Compute the mixture of a list of distribution

    Args:
        listdist (list): List of distributions

    Returns:
        add: Computation of the adding of a distribution to the mixture
    """
    if len(listdist) == 0:
        return lambda x: 0
    else:
        elem = listdist[0]
        return add(scale(elem["coefficient"], elem["distribution"]), mixture(listdist[1:]))


def weight_to_sd(w: int, high_confidence_coefficient: float, low_confidence_coefficient: float) -> float:
    """Transform a the weight into standard deviation value for the mixture computation

    Args:
        w (int): Weight to be transformed
        high_confidence_coefficient (float): High confidence coefficient for the mixture computation.
        low_confidence_coefficient (float): Low confidence coefficient for the mixture computation.

    Returns:
        float: Weight transformed
    """
    return (w/100)*high_confidence_coefficient+(1-w/100)*low_confidence_coefficient


def get_mixture_gene(g: dict, means: dict, high_confidence_coefficient: float, low_confidence_coefficient: float) -> mixture:
    """Get the distribution mixture of a gene

    Args:
        g (dict): Gene predictions and fold change value
        means (dict): Calculated means from fold change dataset
        high_confidence_coefficient (float): Minimum boundary marker.
        low_confidence_coefficient (float): Maximum boundary marker.

    Returns:
        mixture: Mixture of the gene
    """
    sum_ = g["Sign-1"][0]+g["Sign1"][0]+g["Sign0"][0]
    Lmixture = list()
    for sign in ["Sign-1", "Sign1", "Sign0"]:
        mean_ = means[sign]
        sd = weight_to_sd(
            g[sign][1], high_confidence_coefficient, low_confidence_coefficient)
        Lmixture.append(
            {"coefficient": g[sign][0]/sum_, "distribution": fnormal(mean_, sd)})
    return mixture(Lmixture)


def iggy_results_transformation(sign: str) -> dict:
    """Transformation of Iggy's sign into MajS prediction result

    Args:
        sign (str): Iggy's prediction

    Returns:
        dict: MajS prediction format
    """
    if sign == "0":
        return {
            "Sign0": [1, 100, 0],
            "Sign1": [0, 100, 0],
            "Sign-1": [0, 100, 0]
        }
    elif sign == "-":
        return {
            "Sign0": [0, 100, 0],
            "Sign1": [0, 100, 0],
            "Sign-1": [1, 100, 0]
        }
    elif sign == "+":
        return {
            "Sign0": [0, 100, 0],
            "Sign1": [1, 100, 0],
            "Sign-1": [0, 100, 0]
        }
    elif sign == "notPlus":
        return {
            "Sign0": [1, 100, 0],
            "Sign1": [0, 100, 0],
            "Sign-1": [1, 100, 0]
        }
    elif sign == "notMinus":
        return {
            "Sign0": [1, 100, 0],
            "Sign1": [1, 100, 0],
            "Sign-1": [0, 100, 0]
        }
    else:
      # no
        return {
            "Sign0": [0, 0, 0],
            "Sign1": [0, 0, 0],
            "Sign-1": [0, 0, 0]
        }


def compute_probability(dist_mixt: list, fc: float, epsilon: float) -> float:
    """Compute the probalities of a fold change in the a mixture

    Args:
        dist_mixt (list): Distribution mixture
        fc (float): Fold change value
        epsilon (float): Epsilon value for the probability computation

    Returns:
        float: Probability value
    """
    return (2*epsilon)*(dist_mixt(fc-epsilon)+dist_mixt(fc+epsilon))/2


def compute_mixtures_and_scores(g: dict, gene: str, means: dict, epsilon: float, high_confidence_coefficient: float, low_confidence_coefficient: float, outfile: str = None, X: list = np.arange(-1, 1, 0.01)) -> dict:
    """Plot the mixture of a gene for the paper and return significance scores of the two methods

    Args:
        g (dict): Gene predictions and fold change value
        gene (str): Gene name
        means (dict): Sign1 and Sign-1 means for the mixture computation
        epsilon (float): Epsilon to compute probabilites
        high_confidence_coefficient (float): High confidence coefficient for the mixture computation
        low_confidence_coefficient (float): Low confidence coefficient for the mixture computation
        outfile (str, optional): Name of the output graph file. Defaults to None.
        X (list, optional): X-axis values. Defaults to np.arange(-1, 1, 0.01).

    Returns:
        dict: MajS' and Iggy's scores
    """
    if not outfile:
        outfile = gene+'.pdf'

    if 'iggy' not in g.keys():
        return None

    fc = g['logFC']

    # Compute mixtures
    majs_dist_mix = get_mixture_gene(
        g, means, high_confidence_coefficient, low_confidence_coefficient)
    iggy_dist_mix = get_mixture_gene(
        iggy_results_transformation(g['iggy']), means, high_confidence_coefficient, low_confidence_coefficient)

    # Compute probabilities and scores for the fc
    majs_proba = compute_probability(majs_dist_mix, fc, epsilon)
    iggy_proba = compute_probability(iggy_dist_mix, fc, epsilon)

    majs_score = majs_proba/MAX_PROBA
    iggy_score = iggy_proba/MAX_PROBA

    # Compute significance score distributions
    majs_score_distribution = [(compute_probability(majs_dist_mix, x, epsilon)/MAX_PROBA)
                               for x in X]
    iggy_score_distribution = [(compute_probability(iggy_dist_mix, x, epsilon)/MAX_PROBA)
                               for x in X]

    # Plot the mixtures
    plt.plot(X, iggy_score_distribution, color='#E66101', label='Iggy')
    plt.plot(X, majs_score_distribution, color='#5E3C99', label='MajS')
    # Plot logFC
    plt.axvline(x=fc, color='black',
                label=f'logFC ({round(fc,3)})', linestyle=':')

    # Setting of Y axis range
    plt.ylim([0, 1.05])
    # Reordering the labels
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [1, 0, 2]
    plt.legend([handles[i] for i in order], [labels[i] for i in order])
    plt.xlabel('logFC')
    plt.ylabel('Significance score')
    plt.title(gene)
    plt.savefig(outfile)
    plt.close()

    return {'majs_score': majs_score, 'iggy_score': iggy_score}


def plot_scores_for_all_genes(all_predictions: dict, outfile: str, threshold: float = None):
    """Plot MajS' and Iggy's scores

    Args:
        all_predictions (dict): Predictions of both method
        outfile (str): Name of the output graph file
        threshold (float, optional): Threshold value to plot the boundary into the graph. Defaults to None.
    """
    X = all_predictions.keys()
    majs_score = [v['majs_score'] for v in all_predictions.values()]
    iggy_score = list()
    for v in all_predictions.values():
        if 'iggy_score' in v.keys():
            iggy_score.append(v['iggy_score'])
        else:
            iggy_score.append(-1)

    df = pd.DataFrame({
        'gene': X,
        'majs': majs_score,
        'iggy': iggy_score
    })

    df = df.sort_values(by=['majs'])

    fig, ax = plt.subplots(figsize=(60, 40), constrained_layout=True)
    ax.scatter(df['gene'], df['iggy'], color='#E66101',
               label='Iggy', s=1500, alpha=0.75)
    ax.scatter(df['gene'], df['majs'], color='#5E3C99',
               label='MajS', s=1500, alpha=0.75)

    if threshold:
        ax.axhline(y=threshold, color='black',
                   label='threshold', linestyle=':')
        labels_order = [1, 0, 2]
    else:
        labels_order = [1, 0]

    plt.xticks(rotation=90, ha='center', size=50)
    plt.yticks(size=50)
    ax.xaxis.grid(linestyle=':', color='black')
    # Reordering the labels
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[i] for i in labels_order], [labels[i] for i in labels_order],
               borderaxespad=2, fontsize=50, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(outfile)
    plt.close()


def compute_means_sd_fc(filename: str) -> dict:
    """Calculate means and standard deviation of gene values between +/- 'zero_threshold' and +/- 'plus_threshold'

    Args:
        filename (str): Name of the file contains the dataset

    Returns:
        dict: Means and standard deviation of the two subsets
    """
    zero_threshold = 0.15
    plus_threshold = 1.5
    fc_data = pd.read_csv(filename)

    plus_data = fc_data.loc[(fc_data['logFC'] > zero_threshold) & (
        fc_data['logFC'] < plus_threshold)]
    minus_data = fc_data.loc[(fc_data['logFC'] < -zero_threshold)
                             & (fc_data['logFC'] > -plus_threshold)]

    return {
        'plus': {
            'mean': plus_data['logFC'].mean(),
            'sd': plus_data['logFC'].std(),
            'len': len(plus_data)
        },
        'minus': {
            'mean': minus_data['logFC'].mean(),
            'sd': minus_data['logFC'].std(),
            'len': len(minus_data)
        }
    }


def calculate_real_proba(sigma: float, epsilon: float) -> float:
    """Calculate the real probability for specific sigma and espsilon

    Args:
        sigma (float): Sigma value
        epsilon (float): Epsilon value

    Returns:
        float: Calculated probality
    """
    def f(x):
        return 1/(sigma*sqrt(2*pi))*exp(-1/2*pow((x/sigma), 2))

    res, _ = quad(f, -epsilon, epsilon)
    return res


if __name__ == '__main__':
    all_predictions = dict()

    # Parse arguments
    args = parse_args()

    # Fix global variables
    HIGH_CONFIDENCE_COEFFICIENT = 0.05 if not args.high_confidence_coefficient else args.high_confidence_coefficient
    LOW_CONFIDENCE_COEFFICIENT = 0.5 if not args.low_confidence_coefficient else args.low_confidence_coefficient
    EPSILON = 0.005 if not args.epsilon else args.epsilon

    MIN_PROBA = calculate_real_proba(
        LOW_CONFIDENCE_COEFFICIENT, EPSILON)  # Not used
    MAX_PROBA = calculate_real_proba(
        HIGH_CONFIDENCE_COEFFICIENT, EPSILON)  # Pmax in the paper

    # Load data
    all_predictions = load_majs(args.majs_file, all_predictions)
    all_predictions = load_rna(args.rna_file, all_predictions)
    all_predictions = load_iggy(args.iggy_file, all_predictions)
    observed_genes = load_observations(args.obs_file)

    # Check the benchmark number
    benchmark = 'benchmark1' if 'benchmark1' in args.iggy_file else 'benchmark2'

    # Calculate mean and standard deviation of genes from all fold change dataset
    all_fc = compute_means_sd_fc('./All-logFC_RNA.csv')
    print('\n\n ------------\n\n')
    print("### Mean and Std - All fold change dataset ###")
    print(all_fc)

    # Define means for mixture computation
    means = {
        'Sign1': all_fc['plus']['mean'],
        'Sign-1': all_fc['minus']['mean'],
        'Sign0': 0
    }
    # Check if out dir if already created
    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if not os.path.isdir(outdir+'/'+benchmark):
        os.makedirs(outdir+'/'+benchmark)

    outdir = outdir+'/'+benchmark

    # Select only non-observed genes
    non_observed = dict()
    to_del = list()
    for gene, values in all_predictions.items():
        if gene not in observed_genes:
            non_observed[gene] = values

    # Compute mixtures and scores for non-observed genes (Fig. 4 of the article)
    for gene in non_observed.keys():
        gene_score = compute_mixtures_and_scores(
            non_observed[gene], gene, means, EPSILON, HIGH_CONFIDENCE_COEFFICIENT, LOW_CONFIDENCE_COEFFICIENT, outfile=f'{outdir}/{gene}.pdf')
        # No Iggy's prediction
        if not gene_score:
            to_del.append(gene)
        else:
            non_observed[gene].update(gene_score)
    for gene in to_del:
        del non_observed[gene]
    export_to_json(non_observed, f'{outdir}/gene_scores.json')

    # Plot Iggy's and MajS' scores for comparison
    plot_scores_for_all_genes(
        non_observed, f'{outdir}/score_plot_{benchmark}.pdf')

    # Calculate some statistics about both methods
    stats = {'majs>iggy': 0, 'iggy>majs': 0, 'equal': 0, 'no_iggy': 0, 'dist_majs': {
        'sum': 0.0, 'mean': 0.0}, 'dist_iggy': {'sum': 0.0, 'mean': 0.0}}
    for v in non_observed.values():
        if 'iggy_score' in v.keys():
            if v['majs_score'] > v['iggy_score']:
                stats['majs>iggy'] += 1
                stats['dist_majs']['sum'] += (v['majs_score'] -
                                              v['iggy_score'])
            else:
                if v['majs_score'] == v['iggy_score']:
                    stats['equal'] += 1
                else:
                    stats['iggy>majs'] += 1
                    stats['dist_iggy']['sum'] += (v['iggy_score'] -
                                                  v['majs_score'])
        else:
            stats['no_iggy'] += 1
    stats['dist_majs']['mean'] = stats['dist_majs']['sum']/stats['majs>iggy']
    stats['dist_iggy']['mean'] = stats['dist_iggy']['sum']/stats['iggy>majs']

    print('\n\n ------------\n\n')
    print("### Statistics of the comparison of MajS and Iggy ###")
    print(stats)

    # Calculate the number of gene score higher than a threeshold
    threshold = 0.01
    sup_threshold = {'majs': 0, 'iggy': 0, 'no_iggy': 0}
    for v in non_observed.values():
        if 'iggy_score' in v.keys():
            if v['majs_score'] >= threshold:
                sup_threshold['majs'] += 1
            if v['iggy_score'] >= threshold:
                sup_threshold['iggy'] += 1
        else:
            sup_threshold['no_iggy'] += 1

    print('\n\n ------------\n\n')
    print(f'### Gene score higher than the threshold \'{threshold}\' ###')
    print(sup_threshold)
    print('\n\n ------------\n\n')

    # Export non-observed data and calculated scores into JSON file
    if args.export:
        print('Export non-observed data...')
        export_to_json(
            non_observed, f'{outdir}/non_observed_export.json')

    print("DONE.")
