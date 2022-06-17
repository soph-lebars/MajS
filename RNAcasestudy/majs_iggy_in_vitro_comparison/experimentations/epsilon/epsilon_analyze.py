import argparse
import json
import os
import matplotlib.pyplot as plt


def cli_parser() -> argparse.ArgumentParser:
    """Define the CLI parser

    Returns:
        argparse.ArgumentParser: Argument parser
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-path', '-d', type=str, required=True,
                        help='Path to the folder contains experimentd data')
    parser.add_argument('--out-path', '-o', type=str, required=True,
                        help='Output path')
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


def get_scores(epsilon_value: float, path: str, benchmark: str) -> dict:
    """Get scores of genes for a specific epsilon

    Args:
        epsilon_value (float): Espilon value
        path (str): Path to folder contains all experiments data
        benchmark (str): Banchmark number

    Returns:
        dict: Current scores of genes for epsilon
    """
    with open(f'{path}/{epsilon_value}/{benchmark}/gene_scores.json', 'r') as file_:
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


if __name__ == '__main__':
    # Parse arguments
    args = parse_args()

    data_path = args.data_path
    out_path = args.out_path
    benchmark = f'benchmark{args.benchmark}'

    out_path = f'{out_path}/{benchmark}'

    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    if not os.path.isdir(f'{out_path}/score_traces'):
        os.makedirs(f'{out_path}/score_traces')

    # Get scores for each genes for each epsilon
    epsilon_values = os.listdir(data_path)
    epsilon_values.sort()
    scores = dict()
    for epsilon_value in epsilon_values:

        current_scores = get_scores(epsilon_value, data_path, benchmark)
        for gene in current_scores.keys():
            if gene in scores:
                scores[gene]['majs'][epsilon_value] = current_scores[gene]['majs_score']
                scores[gene]['iggy'][epsilon_value] = current_scores[gene]['iggy_score']
            else:
                scores[gene] = {
                    'majs': {epsilon_value: current_scores[gene]['majs_score']},
                    'iggy': {epsilon_value: current_scores[gene]['iggy_score']}
                }
    json.dump(scores, open('scores.json', 'w'))

    epsilon_breaks = dict()
    trace_score = dict()

    for gene in scores:
        current_comp = compare_scores(
            scores[gene]['majs'][epsilon_values[0]], scores[gene]['iggy'][epsilon_values[0]])
        cpt = 1
        # While no change, compare next score
        while current_comp == compare_scores(scores[gene]['majs'][epsilon_values[cpt]], scores[gene]['iggy'][epsilon_values[cpt]]):
            current_comp = compare_scores(
                scores[gene]['majs'][epsilon_values[cpt]], scores[gene]['iggy'][epsilon_values[cpt]])
            cpt += 1
            # If it is the end of scores
            if cpt == len(epsilon_values):
                cpt = "BREAK"
                break
        # If a change is present
        if not cpt == 'BREAK':
            before = compare_scores(
                scores[gene]['majs'][epsilon_values[cpt-1]], scores[gene]['iggy'][epsilon_values[cpt-1]])
            after = compare_scores(
                scores[gene]['majs'][epsilon_values[cpt]], scores[gene]['iggy'][epsilon_values[cpt]])
            epsilon_break = f'{epsilon_values[cpt]} -- FROM {before} TO {after}'
        else:
            epsilon_break = 'NO CHANGE'
        epsilon_breaks[gene] = epsilon_break

        diff_list = list()
        # Compute the difference between the two scores
        for i in range(len(epsilon_values)):
            diff = scores[gene]['majs'][epsilon_values[i]] - \
                scores[gene]['iggy'][epsilon_values[i]]
            diff_list.append(diff)
        trace_score[gene] = diff_list
        json.dump(trace_score, open(f'{out_path}/trace_score.json', 'w'))

        outfile = f'{out_path}/score_traces/{gene}.pdf'

        # Plot the difference in function of the espilon values
        X = epsilon_values
        Y = diff_list
        plt.plot(X, Y)
        eps_break = None
        if 'NO CHANGE' not in epsilon_break:
            eps_break = epsilon_break.split(' -- ')[0]

        if eps_break != None:
            plt.axvline(x=eps_break, color='gray',
                        label=f'epsilon ({str(eps_break)})', linestyle=':')
            plt.legend()
        plt.axhline(y=0, color='black',  linewidth=1)
        plt.xlabel('epsilon')
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='major', labelsize=3)
        plt.ylabel('score(majs) - score(iggy)')
        plt.ylim([-.4, .4])
        plt.title(gene)
        plt.savefig(outfile)
        plt.close()

    json.dump(epsilon_breaks, open(f'{out_path}/epsilon_breaks.json', 'w'))
