from identityholder import IdentityHolder
from minimumtopologies import ispositive
import networkx as nx
import rendergraph
import sys

def countmotifs(network):
    print('Finding cycles')
    cycle_sets = [IdentityHolder(frozenset(cycle), cycle) for cycle in nx.algorithms.simple_cycles(network) if ispositive(network, cycle)]
    cycle_edge_sets = dict()
    for holder in cycle_sets:
        cycle = holder.tag
        edge_set = frozenset((cycle[n], cycle[(n + 1) % len(cycle)]) for n in range(len(cycle)))
        cycle_edge_sets[holder] = edge_set
    print('Creating cycle intersection graphs')
    cycle_graph = nx.Graph()
    cycle_graph.add_nodes_from(cycle_sets)
    cycle_edge_graph = nx.Graph()
    cycle_edge_graph.add_nodes_from(cycle_sets)
    for i1, holder1 in enumerate(cycle_sets):
        for i2 in range(i1 + 1, len(cycle_sets)):
            holder2 = cycle_sets[i2]
            shared_nodes = holder1.value.intersection(holder2.value)
            if len(shared_nodes) > 0:
                cycle_graph.add_edge(holder1, holder2, shared=shared_nodes)
                shared_edges = cycle_edge_sets[holder1].intersection(cycle_edge_sets[holder2])
                if len(shared_edges) > 0:
                    cycle_edge_graph.add_edge(holder1, holder2, shared=shared_edges)
    def coverextracycle(holder1, holder2, ignoring=None):
        for common_neighbor in set(cycle_edge_graph[holder1]).intersection(set(cycle_edge_graph[holder2])):
            if common_neighbor is ignoring:
                continue
            if cycle_edge_sets[common_neighbor] < cycle_edge_sets[holder1].union(cycle_edge_sets[holder2]):
                return True
        return False
    print(len(cycle_sets), 'cycles,', len(cycle_graph.edges), 'node sharings')
    print('Searching for Type I motifs')
    checked = 0
    type1 = 0
    for a, b, c in findtriangles(cycle_graph):
        shared_nodes = cycle_graph.edges[a, b]['shared'].intersection(cycle_graph.edges[b, c]['shared'])
        if len(shared_nodes) == 0:
            continue
        if not shared_nodes.isdisjoint(cycle_graph.edges[a, c]['shared']):
            if coverextracycle(a, b, c) or coverextracycle(a, c, b) or coverextracycle(b, c, a):
                continue # TODO: Avoid rejecting emergent cycles that are necessary to the motif
            type1 += 1
        checked += 1
    print('Searching for Type II motifs')
    checked = 0
    type2 = 0
    for holder in cycle_sets:
        neighbors = list(cycle_graph.neighbors(holder))
        for i1, neigh1 in enumerate(neighbors):
            if coverextracycle(holder, neigh1):
                continue
            for i2 in range(i1 + 1, len(neighbors)):
                neigh2 = neighbors[i2]
                if cycle_graph.has_edge(neigh1, neigh2):
                    continue
                if coverextracycle(holder, neigh2):
                    continue
                type2 += 1
        checked += 1
        print(f'{checked}\r', end='')
    return (len(cycle_sets), type1, type2)

def findtriangles(graph):
    for a, b in graph.edges:
        for c in frozenset(graph[a]).intersection(frozenset(graph[b])):
            if c.isbefore(a) and c.isbefore(b):
                yield (a, b, c)

if __name__ == "__main__":
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(sys.argv[1]))
    result = countmotifs(graph)
    print('PFLs', result[0], '\nType1', result[1], '\nType2', result[2], sep='\t')
    