import argparse
from countmotifs import hasrepression
from identityholder import IdentityHolder
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition
import networkx as nx
import numpy as np
import permutenetwork
import random
from scipy.special import comb
import scipy.stats

def isfused(cycle_sets):
    intersection = cycle_sets[0]
    for cycle in cycle_sets[1:]:
        intersection = intersection.intersection(cycle)
    return len(intersection) > 0

def findconnector(cycle_sets):
    for connector_c, connector in enumerate(cycle_sets):
        other1 = cycle_sets[(connector_c + 1) % 3]
        other2 = cycle_sets[(connector_c + 2) % 3]
        if connector.isdisjoint(other1) or connector.isdisjoint(other2):
            continue
        if other1.isdisjoint(other2):
            return connector_c
    return None

def summarize(graph, samples, max_motif_size=None, max_cycle_length=None):
    cycles = [(frozenset(cycle), ispositive(graph, cycle), cycle, hasrepression(graph, cycle)) for cycle in liuwangcycles.cyclesgenerator(graph, max_cycle_length)]
    pfls = len([0 for cycle in cycles if cycle[1]])
    pfl_ratio = pfls / len(cycles) if len(cycles) > 0 else 0.0
    if len(cycles) < 2:
        return pfls, pfl_ratio, 0, 0, 0, 0, 0, 0, 0
    fused3, bridged2, fusedpfls, type1, type2, misa, excitable, minimisa, fpnp = 0, 0, 0, 0, 0, 0, 0, 0, 0
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
                if chosen[0][3] or chosen[0][3]:
                    minimisa += 1
            elif chosen[0][1] != chosen[1][1]:
                fpnp += 1
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
                excitable += 1
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
    return (pfls, pfl_ratio, type1_est, type2_est, misa * sample3_adjustment, excitable * sample3_adjustment, type1 / fused3, type2 / bridged2,
            fpnp * sample2_adjustment, minimisa * sample2_adjustment, minimisa / fusedpfls)

def evaluate(graph, permutations, samples, base_trials=10, use_full_permutation=True, max_nodes_for_sample=None, max_motif_size=None, max_cycle_length=None, fixed_sign_sources=None):
    base_connected = nx.algorithms.is_strongly_connected(graph)
    base_results = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # PFLs, PFL/FL, Type1, Type2, MISA, Excite, T1/F3, F2/Brdg, FPNP, uMISA, uMISA/F
    for _ in range(base_trials):
        feasible_base_subgraph = graph if max_nodes_for_sample is None else permutenetwork.randomsubgraph(graph, max_nodes_for_sample)
        for component, value in enumerate(summarize(feasible_base_subgraph, samples, max_motif_size, max_cycle_length)):
            base_results[component] += value
    for component in range(len(base_results)):
        base_results[component] /= base_trials
    print('Base results:', base_results)
    permutation_results = [[] for _ in base_results]
    for checked_permutations, permutation in permutenetwork.generatepermutations(graph, base_connected, use_full_permutation, max_nodes_for_sample, fixed_sign_sources):
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
    return empirical_cdfs, p_values

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='input GraphML file to process')
    parser.add_argument('permutations', type=int, help='number of network permutations to generate')
    parser.add_argument('samples', type=int, help='number of cycle triplets to sample from each permutation')
    parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
    parser.add_argument('--basetrials', type=int, default=10, help='number of trials to average for the input network')
    parser.add_argument('--desonly', action='store_true', help='use only double-edge swaps for permutation')
    parser.add_argument('--maxcycle', type=int, help='maximum number of nodes in a cycle')
    parser.add_argument('--maxmotifsize', type=int, help='maximum number of nodes in a motif')
    parser.add_argument('--fixedsign', type=str, help='regex matching nodes whose source regulations to preserve signs of')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    empirical_cdfs, p_values = evaluate(graph, args.permutations, args.samples, args.basetrials, not args.desonly, args.maxnodes, args.maxmotifsize, args.maxcycle, args.fixedsign)
    column_names = ['PFLs', 'PFL/FL', 'Type 1', 'Type 2', 'MISA', 'Excite', 'T1/Fus3', 'T2/Brdg', 'FPNP', 'uMISA', 'uMISA/F']
    print('\tFracLE')
    for i, column in enumerate(column_names):
        print(column, empirical_cdfs[i], sep='\t')
