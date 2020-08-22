import networkx as nx
import permutenetwork
import scipy.stats
import sys

graph = nx.read_graphml(sys.argv[1])
possible_edges = len(graph.nodes)**2 # NetworkX's density uses n(n - 1) in the denominator instead
edge_existence_chance = len(graph.edges) / possible_edges
samples = 150
expected_values = [edge_existence_chance * samples] * possible_edges
permutations = [graph.copy() for _ in range(samples)]
swaps = 0
for n in range(30):
    for _ in range(5):
        permutenetwork.edgeswap(permutations)
    swaps += 5
    observed_values = []
    for n1 in graph.nodes:
        for n2 in graph.nodes:
            present = 0
            for permutation in permutations:
                if permutation.has_edge(n1, n2):
                    present += 1
            observed_values.append(present)
    chi2, pval = scipy.stats.chisquare(observed_values, expected_values)
    print(swaps, chi2, pval, sep='\t')
