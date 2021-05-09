from collections import defaultdict
import networkx as nx
import sys

def loadcitedtsv(reader):
    """Creates a NetworkX DiGraph from a reader of a TRRUST-style TSV."""
    evidence = defaultdict(int)
    directions = {'Activation': 1, 'Repression': -1, 'Unknown': 0}
    graph = nx.DiGraph()
    for line in reader:
        tf, target, direction, refs = line.rstrip().split('\t')
        evidence[tf, target] += directions[direction] * (refs.count(';') + 1)
        graph.add_node(tf)
        graph.add_node(target)
    for (tf, target), balance in evidence.items():
        if balance == 0:
            continue
        graph.add_edge(tf, target, repress=(balance < 0))
    for node in graph.nodes:
        graph.nodes[node]['name'] = node
    return nx.convert_node_labels_to_integers(graph)

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        network = loadcitedtsv(f)
        nx.write_graphml(network, sys.argv[2])
