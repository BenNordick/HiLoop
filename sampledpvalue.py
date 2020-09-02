import argparse
from collections import deque
from identityholder import IdentityHolder
from minimumtopologies import ispositive
import networkx as nx
import permutenetwork
import random

def istype1(cycle_sets):
    intersection = cycle_sets[0]
    for cycle in cycle_sets[1:]:
        intersection = intersection.intersection(cycle)
    return len(intersection) > 0

def istype2(cycle_sets):
    for connector_c, connector in enumerate(cycle_sets):
        other1 = cycle_sets[(connector_c + 1) % 3]
        other2 = cycle_sets[(connector_c + 2) % 3]
        if connector.isdisjoint(other1) or connector.isdisjoint(other2):
            continue
        if other1.isdisjoint(other2):
           return True
    return False

def summarize(graph, samples):
    cycles = [(IdentityHolder(frozenset(cycle)), ispositive(graph, cycle)) for cycle in nx.algorithms.simple_cycles(graph)]
    pfls = len([0 for cycle in cycles if cycle[1]])
    if len(cycles) < 3:
        return pfls, 0, 0
    type1, type2 = 0, 0
    checked = set()
    for _ in range(samples):
        chosen = random.sample(cycles, 3)
        if frozenset(chosen) in checked:
            continue
        checked.add(frozenset(chosen))
        if [entry[1] for entry in chosen] != [True] * 3:
            continue
        chosen_cycles = [entry[0].value for entry in chosen]
        if istype1(chosen_cycles):
            type1 += 1
        elif istype2(chosen_cycles):
            type2 += 1
    return pfls, type1, type2

def randomsubgraph(graph, max_nodes):
    queue = deque(maxlen=len(graph.nodes))
    selected = set()
    queue.append(random.sample(graph.nodes, 1)[0])
    while len(selected) < max_nodes and len(queue) > 0:
        head = queue.popleft()
        if head in selected:
            continue
        selected.add(head)
        neighbors = list(graph.successors(head))
        random.shuffle(neighbors)
        for n in neighbors:
            queue.append(n)
    return graph.subgraph(selected)

def evaluate(graph, permutations, samples, base_trials=10, use_full_permutation=True, max_nodes_for_sample=None):
    if max_nodes_for_sample is None:
        base_samples = summarize(graph, samples * base_trials)
        base_results = (base_samples[0], base_samples[1] / base_trials, base_samples[2] / base_trials)
    else:
        base_results = [0, 0, 0]
        for _ in range(base_trials):
            feasible_base_subgraph = randomsubgraph(graph, max_nodes_for_sample)
            for component, value in enumerate(summarize(feasible_base_subgraph, samples)):
                base_results[component] += value / base_trials
    print('Base results:', base_results)
    permutation_results = [[], [], []]
    last_permutation = graph
    for n in range(permutations):
        if n % 20 == 0 and use_full_permutation:
            last_permutation = permutenetwork.permutenetwork(graph)
        else:
            last_permutation = permutenetwork.permuteedgeswaps(permutenetwork.permuteregulations(last_permutation))
        feasible_subgraph = last_permutation if max_nodes_for_sample is None else randomsubgraph(last_permutation, max_nodes_for_sample)
        for component, value in enumerate(summarize(feasible_subgraph, samples)):
            permutation_results[component].append(value)
    p_values = []
    for component, value in enumerate(base_results):
        extreme_counts = 0
        for permutation_value in permutation_results[component]:
            if permutation_value >= value:
                extreme_counts += 1
        p_values.append(extreme_counts / permutations)
    return p_values

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='input GraphML file to process')
    parser.add_argument('permutations', type=int, help='number of network permutations to generate')
    parser.add_argument('samples', type=int, help='number of cycle triplets to sample from each permutation')
    parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
    parser.add_argument('--basetrials', type=int, default=10, help='number of trials to average for the input network')
    parser.add_argument('--desonly', action='store_true', help='use only double-edge swaps for permutation')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    p_values = evaluate(graph, args.permutations, args.samples, args.basetrials, not args.desonly, args.maxnodes)
    print('PFLs:', p_values[0], sep='\t')
    print('Type 1:', p_values[1], sep='\t')
    print('Type 2:', p_values[2], sep='\t')
