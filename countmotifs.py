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
    #print('Finding triangles')
    #cycle_triangles = findtriangles(cycle_graph)
    print('Searching for Type II motifs')
    checked = 0
    type2 = 0
    for holder in cycle_sets:
        neighbors = list(cycle_graph.neighbors(holder))
        edge_sharing_neighbors = set(cycle_edge_graph[holder])
        for i1, neigh1 in enumerate(neighbors):
            esns1 = set(cycle_edge_graph[neigh1])
            extra_cycle = False
            for esn in esns1.intersection(edge_sharing_neighbors):
                if cycle_edge_sets[esn] < cycle_edge_sets[holder].union(cycle_edge_sets[neigh1]):
                    extra_cycle = True
                    break
            if extra_cycle:
                continue
            for i2 in range(i1 + 1, len(neighbors)):
                neigh2 = neighbors[i2]
                if cycle_graph.has_edge(neigh1, neigh2):
                    continue
                esns2 = set(cycle_edge_graph[neigh2])
                extra_cycle = False
                for esn in esns2.intersection(edge_sharing_neighbors):
                    if cycle_edge_sets[esn] < cycle_edge_sets[holder].union(cycle_edge_sets[neigh2]):
                        extra_cycle = True
                        break
                if extra_cycle:
                    continue
                type2 += 1
        checked += 1
        print(f'{checked}\r', end='')
    # TODO: Type 1         
    return (len(cycle_sets), 0, type2)

def findtriangles(graph):
    triangles = set()
    for a, b in graph.edges:
        for c in frozenset(graph[a]).intersection(frozenset(graph[b])):
            triangles.add(frozenset([a, b, c]))
    return triangles

if __name__ == "__main__":
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(sys.argv[1]))
    result = countmotifs(graph)
    print('PFLs', result[0], '\nType1', result[1], '\nType2', result[2], sep='\t')
    