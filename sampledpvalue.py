from minimumtopologies import ispositive
import networkx as nx
from permutenetwork import permutenetwork
import random
import sys

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
    cycles = [(frozenset(cycle), ispositive(graph, cycle)) for cycle in nx.algorithms.simple_cycles(graph)]
    type1, type2 = 0, 0
    for _ in range(samples):
        chosen = random.sample(cycles, 3)
        if [entry[1] for entry in chosen] != [True] * 3:
            continue
        chosen_cycles = [entry[0] for entry in chosen]
        if istype1(chosen_cycles):
            type1 += 1
        elif istype2(chosen_cycles):
            type2 += 1
    pfls = len([0 for cycle in cycles if cycle[1]])
    return pfls, type1, type2

def evaluate(graph, permutations, samples):
    base_results = summarize(graph, samples)
    print('Base results:', base_results)
    permutation_results = [[], [], []]
    for _ in range(permutations):
        permutation = permutenetwork(graph)
        for component, value in enumerate(summarize(permutation, samples)):
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
    graph = nx.read_graphml(sys.argv[1])
    p_values = evaluate(graph, int(sys.argv[2]), int(sys.argv[3]))
    print('PFLs:', p_values[0], sep='\t')
    print('Type 1:', p_values[1], sep='\t')
    print('Type 2:', p_values[2], sep='\t')
