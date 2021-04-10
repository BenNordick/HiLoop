import argparse
from collections import Counter
from identityholder import IdentityHolder
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition
import networkx as nx
import rendergraph
import sys

launched_specifically = False

def countmotifs(network, max_cycle_length=None, max_motif_size=None, check_nfl=False):
    if launched_specifically:
        print('Finding cycles')
    cycle_sets = [IdentityHolder(frozenset(cycle), (cycle, ispositive(network, cycle), hasrepression(network, cycle))) for
                  cycle in liuwangcycles.generatecycles(network, max_cycle_length) if check_nfl or ispositive(network, cycle)]
    cycle_edge_sets = dict()
    for holder in cycle_sets:
        cycle = holder.tag[0]
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
    def coverextrapfl(holder1, holder2):
        for common_neighbor in set(cycle_edge_graph[holder1]).intersection(set(cycle_edge_graph[holder2])):
            if common_neighbor.tag[1] and cycle_edge_sets[common_neighbor] < cycle_edge_sets[holder1].union(cycle_edge_sets[holder2]):
                return True
        return False
    if launched_specifically:
        print(len(cycle_sets), 'cycles,', len(cycle_graph.edges), 'node sharings')
        print('Searching for Type I and excitable motifs')
    type1 = 0
    excitable = 0
    for a, b, c in findtriangles(cycle_graph):
        possible_type1 = (not check_nfl) or (a.tag[1] and b.tag[1] and c.tag[1])
        if check_nfl and not (a.tag[1] or b.tag[1] or c.tag[1]):
            continue
        shared_nodes = cycle_graph.edges[a, b]['shared'].intersection(cycle_graph.edges[b, c]['shared'])
        if len(shared_nodes) == 0:
            continue
        if not shared_nodes.isdisjoint(cycle_graph.edges[a, c]['shared']):
            if max_motif_size:
                all_nodes = a.value.union(b.value).union(c.value)
                if len(all_nodes) > max_motif_size:
                    continue
            extra_cycles = [*findinducedcycles([a, b], c), *findinducedcycles([a, c], b), *findinducedcycles([b, c], a), *findinducedcycles([a, b, c])]
            relevant_extras = [c for c in extra_cycles if c.tag[1]] if possible_type1 else extra_cycles
            if len(relevant_extras) > 0:
                double_counting = False
                for extra in relevant_extras:
                    if extra.isbefore(c):
                        double_counting = True
                        break
                if double_counting:
                    continue
                all_cycles = [a, b, c] + relevant_extras
                all_edges = cycle_edge_sets[a].union(cycle_edge_sets[b]).union(cycle_edge_sets[c])
                found_extra = False
                for edge in all_edges:
                    node_uses_after_elimination = Counter()
                    pfls_after_elimination = 0
                    nfls_after_elimination = 0
                    for cycle in all_cycles:
                        if edge not in cycle_edge_sets[cycle]:
                            node_uses_after_elimination.update(cycle.value)
                            if cycle.tag[1]:
                                pfls_after_elimination += 1
                            else:
                                nfls_after_elimination += 1
                    still_excitable = pfls_after_elimination > 0 and nfls_after_elimination > 0
                    if (possible_type1 or not still_excitable) and len(node_uses_after_elimination) > 0 and max(node_uses_after_elimination.values()) >= 3:
                        found_extra = True
                        break
                if found_extra:
                    continue
            if possible_type1:
                type1 += 1
            else:
                excitable += 1
    if launched_specifically:
        print('Checking fused pairs')
    minimisa = 0
    fpnp = 0
    for holder1, holder2 in cycle_graph.edges:
        if max_motif_size:
            all_nodes = holder1.value.union(holder2.value)
            if len(all_nodes) > max_motif_size:
                continue
        if holder1.tag[1] and holder2.tag[1] and (holder1.tag[2] or holder2.tag[2]):
            minimisa += 1
        elif check_nfl and holder1.tag[1] != holder2.tag[1]:
            fpnp += 1
    if launched_specifically:
        print('Searching for Type II and MISA motifs')
    checked = 0
    type2 = 0
    misa = 0
    for holder in cycle_sets:
        if not holder.tag[1]:
            continue
        neighbors = [h for h in cycle_graph.neighbors(holder) if h.tag[1]]
        for i1, neigh1 in enumerate(neighbors):
            if coverextrapfl(holder, neigh1):
                continue
            for i2 in range(i1 + 1, len(neighbors)):
                neigh2 = neighbors[i2]
                if max_motif_size:
                    all_nodes = holder.value.union(neigh1.value).union(neigh2.value)
                    if len(all_nodes) > max_motif_size:
                        continue
                if cycle_graph.has_edge(neigh1, neigh2):
                    continue
                if coverextrapfl(holder, neigh2):
                    continue
                type2 += 1
                if ismutualinhibition(network, holder.tag[0], neigh1.value, neigh2.value):
                    misa += 1
        if launched_specifically:
            checked += 1
            print(f'{checked}\r', end='')
    pfls = sum(1 for holder in cycle_sets if holder.tag[1])
    return (pfls, type1, type2, misa, minimisa, excitable, fpnp)

def hasrepression(graph, cycle):
    return any(graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]['repress'] for i in range(len(cycle)))

def findtriangles(graph):
    checked = 0
    for a, b in graph.edges:
        for c in frozenset(graph[a]).intersection(frozenset(graph[b])):
            if a.isbefore(c) and b.isbefore(c):
                yield (a, b, c)
        if launched_specifically:
            checked += 1
            if checked % 100 == 0:
                print(f'{checked}\r', end='')

if __name__ == "__main__":
    launched_specifically = True
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='GraphML file to process')
    parser.add_argument('--maxcycle', type=int, help='maximum cycle length')
    parser.add_argument('--maxnodes', type=int, help='maximum number of nodes in a motif')
    parser.add_argument('--checknfl', action='store_true', help='also search for motifs involving negative feedback')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    result = countmotifs(graph, args.maxcycle, args.maxnodes, args.checknfl)
    print('PFLs', result[0], '\nType1', result[1], '\nType2', result[2], '\nMISA', result[3], '\nuMISA', result[4], sep='\t')
    if args.checknfl:
        print('Excit', result[5], '\nFuse+-', result[6], sep='\t')
    