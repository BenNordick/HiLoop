import argparse
from collections import Counter
from identityholder import IdentityHolder
import liuwangcycles
from minimumtopologies import ispositive
import networkx as nx
import rendergraph
import sys

launched_specifically = False

def countmotifs(network, max_cycle_length=None, max_motif_size=None):
    if launched_specifically:
        print('Finding cycles')
    cycle_sets = [IdentityHolder(frozenset(cycle), cycle) for cycle in liuwangcycles.generatecycles(network, max_cycle_length) if ispositive(network, cycle)]
    cycle_edge_sets = dict()
    for holder in cycle_sets:
        cycle = holder.tag
        edge_set = frozenset((cycle[n], cycle[(n + 1) % len(cycle)]) for n in range(len(cycle)))
        cycle_edge_sets[holder] = edge_set
    if launched_specifically:
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
    def findinducedcycles(holders, ignoring=None):
        common_neighbors = set(cycle_edge_graph[holders[0]])
        for cn in range(1, len(holders)):
            common_neighbors.intersection_update(set(cycle_edge_graph[holders[cn]]))
        for common_neighbor in common_neighbors:
            if common_neighbor is ignoring:
                continue
            all_edges = cycle_edge_sets[holders[0]]
            for cn in range(1, len(holders)):
                all_edges = all_edges.union(cycle_edge_sets[holders[cn]])
            if cycle_edge_sets[common_neighbor] < all_edges:
                yield common_neighbor
    def coverextracycle(holder1, holder2):
        for common_neighbor in set(cycle_edge_graph[holder1]).intersection(set(cycle_edge_graph[holder2])):
            if cycle_edge_sets[common_neighbor] < cycle_edge_sets[holder1].union(cycle_edge_sets[holder2]):
                return True
        return False
    if launched_specifically:
        print(len(cycle_sets), 'cycles,', len(cycle_graph.edges), 'node sharings')
        print('Searching for Type I motifs')
    type1 = 0
    for a, b, c in findtriangles(cycle_graph):
        shared_nodes = cycle_graph.edges[a, b]['shared'].intersection(cycle_graph.edges[b, c]['shared'])
        if len(shared_nodes) == 0:
            continue
        if not shared_nodes.isdisjoint(cycle_graph.edges[a, c]['shared']):
            if max_motif_size:
                all_nodes = a.value.union(b.value).union(c.value)
                if len(all_nodes) > max_motif_size:
                    continue
            extra_cycles = [*findinducedcycles([a, b], c), *findinducedcycles([a, c], b), *findinducedcycles([b, c], a), *findinducedcycles([a, b, c])]
            if len(extra_cycles) > 0:
                double_counting = False
                for extra in extra_cycles:
                    if extra.isbefore(c):
                        double_counting = True
                        break
                if double_counting:
                    continue
                all_cycles = [a, b, c] + extra_cycles
                all_edges = cycle_edge_sets[a].union(cycle_edge_sets[b]).union(cycle_edge_sets[c])
                found_extra = False
                for edge in all_edges:
                    node_uses_after_elimination = Counter()
                    for cycle in all_cycles:
                        if edge not in cycle_edge_sets[cycle]:
                            node_uses_after_elimination.update(cycle.value)
                    if len(node_uses_after_elimination) > 0 and max(node_uses_after_elimination.values()) >= 3:
                        found_extra = True
                        break
                if found_extra:
                    continue
            type1 += 1
    if launched_specifically:
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
                if max_motif_size:
                    all_nodes = holder.value.union(neigh1.value).union(neigh2.value)
                    if len(all_nodes) > max_motif_size:
                        continue
                if cycle_graph.has_edge(neigh1, neigh2):
                    continue
                if coverextracycle(holder, neigh2):
                    continue
                type2 += 1
        if launched_specifically:
            checked += 1
            print(f'{checked}\r', end='')
    return (len(cycle_sets), type1, type2)

def findtriangles(graph):
    checked = 0
    for a, b in graph.edges:
        for c in frozenset(graph[a]).intersection(frozenset(graph[b])):
            if a.isbefore(c) and b.isbefore(c):
                yield (a, b, c)
        if launched_specifically:
            checked += 1
            if checked % 10 == 0:
                print(f'{checked}\r', end='')

if __name__ == "__main__":
    launched_specifically = True
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='GraphML file to process')
    parser.add_argument('--maxcycle', type=int, help='maximum cycle length')
    parser.add_argument('--maxnodes', type=int, help='maximum number of nodes in a motif')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    result = countmotifs(graph, args.maxcycle, args.maxnodes)
    print('PFLs', result[0], '\nType1', result[1], '\nType2', result[2], sep='\t')
    