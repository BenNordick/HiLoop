import argparse
from countmotifs import hasrepression
from identityholder import IdentityHolder
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition
import networkx as nx
import numpy as np
import pandas as pd
import permutenetwork
import random
from scipy.special import comb
import scipy.stats

# PyPy virtual environment recommended for performance

def isfused(cycle_sets):
    """Determine whether all cycles (represented as sets of node IDs) share at least one node."""
    intersection = cycle_sets[0]
    for cycle in cycle_sets[1:]:
        intersection = intersection.intersection(cycle)
    return len(intersection) > 0

def findconnector(cycle_sets):
    """Determine the index of the connector cycle for a Type-II-like motif, or None if the cycles do not comprise such a motif."""
    for connector_c, connector in enumerate(cycle_sets):
        other1 = cycle_sets[(connector_c + 1) % 3]
        other2 = cycle_sets[(connector_c + 2) % 3]
        if connector.isdisjoint(other1) or connector.isdisjoint(other2):
            continue
        if other1.isdisjoint(other2):
            return connector_c
    return None

def summarize(graph, samples, max_motif_size=None, max_cycle_length=None):
    """
    Gather high-feedback-related metrics on one network using a sampling approach.

    Arguments:
    - graph: NetworkX DiGraph
    - samples: number of cycle tuples to sample
    - max_motif_size: maximum number of nodes in a motif
    - max_cycle_length: maximum length of cycle to enumerate

    Returns a tuple of metrics: PFLs, PFL/FL, Type1, Type2, MISA, MixHF, T1/F3, F2/Brdg, Excite, MISSA, MISSA/F, uMISSA.
    PFL count and feedback loop positivity ratio are computed exactly; others are based on sampling.
    """
    def safediv(a, b):
        return a / b if b != 0 else 0.0
    cycles = [(frozenset(cycle), ispositive(graph, cycle), cycle, hasrepression(graph, cycle)) for cycle in liuwangcycles.cyclesgenerator(graph, max_cycle_length)]
    pfls = len([0 for cycle in cycles if cycle[1]])
    pfl_ratio = safediv(pfls, len(cycles))
    if len(cycles) < 2:
        return pfls, pfl_ratio, 0, 0, 0, 0, 0, 0, 0, 0
    fused3, bridged2, fusedpfls, type1, type2, misa, mixed, missa, minimissa, excitable = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    sample_size = min(3, len(cycles))
    for _ in range(samples):
        chosen = random.sample(cycles, sample_size)
        chosen_cycles = [entry[0] for entry in chosen]
        if max_motif_size:
            cycle_pair_nodes = chosen_cycles[0].union(chosen_cycles[1])
            if len(cycle_pair_nodes) > max_motif_size:
                continue
        if isfused(chosen_cycles[:2]):
            if chosen[0][1] and chosen[1][1]:
                fusedpfls += 1
                if chosen[0][3] != chosen[1][3]:
                    missa += 1
                    if len(chosen_cycles[0]) == 1 or len(chosen_cycles[1]) == 1:
                        minimissa += 1
            elif chosen[0][1] != chosen[1][1]:
                excitable += 1
        if sample_size < 3:
            continue
        if max_motif_size:
            all_nodes = cycle_pair_nodes.union(chosen_cycles[2])
            if len(all_nodes) > max_motif_size:
                continue
        all_positive = all(entry[1] for entry in chosen)
        if isfused(chosen_cycles):
            fused3 += 1
            if all_positive:
                type1 += 1
            elif any(entry[1] for entry in chosen):
                mixed += 1
        else:
            connector_index = findconnector(chosen_cycles)
            if connector_index is not None:
                bridged2 += 1
                if all_positive:
                    type2 += 1
                    if ismutualinhibition(graph, chosen[connector_index][2], chosen_cycles[(connector_index + 1) % 3], chosen_cycles[(connector_index + 2) % 3]):
                        misa += 1
    sample3_adjustment = comb(len(cycles), 3) / samples
    type1_est = type1 * sample3_adjustment
    type2_est = type2 * sample3_adjustment
    sample2_adjustment = comb(len(cycles), 2) / samples
    return (pfls, pfl_ratio, type1_est, type2_est, misa * sample3_adjustment, mixed * sample3_adjustment, safediv(type1, fused3), safediv(type2, bridged2),
            excitable * sample2_adjustment, missa * sample2_adjustment, safediv(missa, fusedpfls), minimissa * sample2_adjustment)

def evaluate(graph, permutations, samples, base_trials=10, use_full_permutation=True, max_nodes_for_sample=None, max_motif_size=None, max_cycle_length=None, 
             fixed_sign_sources=None, try_scc=False, base_callback=None):
    """
    Evaluate the enrichment of a network in high-feedback-related metrics relative to permutations.

    Arguments:
    - graph: original network as a NetworkX DiGraph
    - permutations: number of permutations to compare against
    - samples: number of samples to take per permutation
    - base_trials: how many sampling (summarize function) runs of the original network to average together
    - use_full_permutation: whether to use rejection sampling of directed configuration models (only suitable for small networks)
    - max_nodes_for_sample: use subgraphs of permuted networks with limited size (i.e. nested sampling, very much not recommended)
    - max_motif_size: maximum motif size in nodes to consider
    - max_cycle_length: maximum length of cycle to enumerate (must be limited for moderately large networks)
    - fixed_sign_sources: regular expression to match names of nodes whose outgoing edge signs should remain constant (e.g. miR)
    - try_scc: whether to require strongly connected permutations of strongly connected input (not recommended)
    - base_callback: callable to receive the list of metrics for the base network early

    Returns a tuple:
    - list of enrichment empirical p-values (proportion of permutations meeting or exceeding the base network in each metric)
    - 1-sample t-test result for the hypothesis that the base result is the mean of the permutation population's for each metric
    - raw results as list of lists, first index specifying metric, second index specifying permutation
    - list of base network's metrics

    Metrics are ordered as listed in the summarize docstring and explained in the README.
    """
    require_scc = try_scc and nx.algorithms.is_strongly_connected(graph)
    base_results = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # PFLs, PFL/FL, Type1, Type2, MISA, MixHF, T1/F3, F2/Brdg, Excite, MISSA, MISSA/F, uMISSA
    for _ in range(base_trials):
        feasible_base_subgraph = graph if max_nodes_for_sample is None else permutenetwork.randomsubgraph(graph, max_nodes_for_sample)
        for component, value in enumerate(summarize(feasible_base_subgraph, samples, max_motif_size, max_cycle_length)):
            base_results[component] += value
    for component in range(len(base_results)):
        base_results[component] /= base_trials
    if base_callback is not None:
        base_callback(base_results)
    permutation_results = [[] for _ in base_results]
    for checked_permutations, permutation in permutenetwork.generatepermutations(graph, require_scc, use_full_permutation, max_nodes_for_sample, fixed_sign_sources):
        for component, value in enumerate(summarize(permutation, samples, max_motif_size, max_cycle_length)):
            permutation_results[component].append(value)
        if checked_permutations == permutations:
            break
    empirical_cdfs = []
    p_values = []
    for component, value in enumerate(base_results):
        extreme_counts = 0
        for permutation_value in permutation_results[component]:
            if permutation_value >= value:
                extreme_counts += 1
        empirical_cdfs.append(extreme_counts / permutations)
        t, p = scipy.stats.ttest_1samp(permutation_results[component], value)
        p_values.append(p)
    return empirical_cdfs, p_values, permutation_results, base_results

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='input GraphML file to process')
    parser.add_argument('permutations', type=int, help='number of network permutations to generate')
    parser.add_argument('samples', type=int, help='number of cycle triplets to sample from each permutation')
    parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
    parser.add_argument('--basetrials', type=int, default=10, help='number of trials to average for the input network')
    parser.add_argument('--dcm', action='store_true', help='use directed configuration model for permutation')
    parser.add_argument('--maxcycle', type=int, help='maximum number of nodes in a cycle')
    parser.add_argument('--maxmotifsize', type=int, help='maximum number of nodes in a motif')
    parser.add_argument('--fixedsign', type=str, help='regex matching nodes whose source regulations to preserve signs of')
    parser.add_argument('--saveraw', type=str, help='CSV file to save raw sample results in')
    parser.add_argument('--scc', action='store_true', help='only consider strongly connected permutations')
    parser.add_argument('--showbase', action='store_true', help='display the raw base results list')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    def callback(base_results):
        if args.showbase:
            print('Base results:', base_results)
    empirical_cdfs, p_values, raw_results, base_results = evaluate(graph, args.permutations, args.samples, args.basetrials, args.dcm,
                                                                   args.maxnodes, args.maxmotifsize, args.maxcycle, args.fixedsign, args.scc, callback)
    column_names = ['PFLs', 'PFL/FL', 'Type 1', 'Type 2', 'MISA', 'MixHF', 'T1/Fus3', 'T2/Brdg', 'Excite', 'MISSA', 'MISSA/F', 'uMISSA']
    print('\tFracLE')
    for i, column in enumerate(column_names):
        print(column, empirical_cdfs[i], sep='\t')
    if args.saveraw:
        arr = np.transpose(np.array(raw_results))
        arr = np.vstack((np.array(base_results), arr))
        df = pd.DataFrame(arr, columns=column_names)
        df['IsBase'] = False
        df.loc[0, 'IsBase'] = True
        df.to_csv(args.saveraw, index=False)
