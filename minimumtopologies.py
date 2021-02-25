from disjointcache import DisjointCache
from identityholder import IdentityHolder
import networkx as nx

def reducetopologies(graph):
    '''Produces a DAG of reduced topologies containing at least one Type I or II motif.'''
    graph = nx.convert_node_labels_to_integers(graph)
    cycle_sets = preparecyclesets(graph)
    atlas = nx.DiGraph()
    disjoint_cache = DisjointCache(len(graph.nodes))
    root_motifs = getmotifs(cycle_sets, len(graph.nodes), disjoint_cache)
    if (root_motifs == (False, False)):
        return atlas
    no_nodes_removed = frozenset()
    atlas.add_node(no_nodes_removed, topology=graph, motifs=root_motifs)
    tryreduce(atlas, graph, set(), disjoint_cache, no_nodes_removed, cycle_sets, graph)
    return atlas

# All below functions are implementation details

def tryreduce(atlas, original_graph, checked_removals, disjoint_cache, removed_edges, cycle_sets, graph):
    for edge in graph.edges:
        new_removed_edges = removed_edges.union(frozenset([edge]))
        if (new_removed_edges in checked_removals):
            if (new_removed_edges in atlas):
                atlas.add_edge(removed_edges, new_removed_edges, removed_edge=edge)
            continue
        checked_removals.add(new_removed_edges)
        reduced_graph = nx.classes.function.restricted_view(original_graph, [], new_removed_edges)
        remaining_cycles = cycle_sets.difference(graph.edges[edge]['cycles'])
        motifs = getmotifs(remaining_cycles, len(graph.nodes), disjoint_cache)
        if (motifs == (False, False)):
            continue
        atlas.add_node(new_removed_edges, topology=reduced_graph, motifs=motifs)
        atlas.add_edge(removed_edges, new_removed_edges, removed_edge=edge)
        tryreduce(atlas, original_graph, checked_removals, disjoint_cache, new_removed_edges, remaining_cycles, reduced_graph)

def preparecyclesets(graph):
    for edge in graph.edges:
        graph.edges[edge]['cycles'] = set()
    all_pfls = [cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if ispositive(graph, cycle)]
    cycle_sets = [IdentityHolder(frozenset(cycle)) for cycle in all_pfls]
    for index, cycle in enumerate(all_pfls):
        for n in range(len(cycle)):
            graph.edges[cycle[n], cycle[(n + 1) % len(cycle)]]['cycles'].add(cycle_sets[index])
    return frozenset(cycle_sets)

def getmotifs(cycle_sets, node_count, disjoint_cache):
    return (hastype1(cycle_sets, node_count), hastype2(cycle_sets, disjoint_cache))

def hastype1(cycle_sets, node_count):
    node_uses = [0] * node_count
    for cycle_holder in cycle_sets:
        for n in cycle_holder.value:
            node_uses[n] += 1
            if (node_uses[n] == 3):
                return True
    return False

def hastype2(cycle_sets, disjoint_cache):
    cycle_sets_list = list(cycle_sets)
    for connector_c in range(len(cycle_sets_list)):
        connector = cycle_sets_list[connector_c]
        for c1 in range(len(cycle_sets_list)):
            if (c1 == connector_c):
                continue
            cycle1 = cycle_sets_list[c1]
            if (disjoint_cache.isdisjoint(connector, cycle1)):
                continue
            for c2 in range(c1 + 1, len(cycle_sets_list)):
                if (c2 == connector_c):
                    continue
                cycle2 = cycle_sets_list[c2]
                if ((not disjoint_cache.isdisjoint(connector, cycle2)) and disjoint_cache.isdisjoint(cycle1, cycle2)):
                    return True
    return False

def ispositive(graph, cycle):
    positive = True
    for i in range(len(cycle)):
        positive ^= graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]['repress']
    return positive
    
def ismutualinhibition(graph, connector_cycle, cycle1, cycle2):
    cycle1 = cycle1 if isinstance(cycle1, frozenset) else frozenset(cycle1)
    cycle2 = cycle2 if isinstance(cycle2, frozenset) else frozenset(cycle2)
    cc_len = len(connector_cycle)
    i = 0
    while not (connector_cycle[i] in cycle1 and connector_cycle[(i + 1) % cc_len] not in cycle1):
        i += 1
    repress = False
    while connector_cycle[i % cc_len] not in cycle2:
        repress ^= graph.edges[connector_cycle[i % cc_len], connector_cycle[(i + 1) % cc_len]]['repress']
        i += 1
    if not repress:
        return False
    while connector_cycle[(i + 1) % cc_len] in cycle2:
        i += 1
    repress = False
    while connector_cycle[i % cc_len] not in cycle1:
        repress ^= graph.edges[connector_cycle[i % cc_len], connector_cycle[(i + 1) % cc_len]]['repress']
        i += 1
    return repress

def ishalfinhibition(graph, connector_cycle, source_cycle, target_cycle):
    cc_len = len(connector_cycle)
    i = 0
    while not (connector_cycle[i] in source_cycle and connector_cycle[(i + 1) % cc_len] not in source_cycle):
        i += 1
    repress = False
    while connector_cycle[i % cc_len] not in target_cycle:
        repress ^= graph.edges[connector_cycle[i % cc_len], connector_cycle[(i + 1) % cc_len]]['repress']
        i += 1
    return repress
