from identityholder import IdentityHolder
import networkx as nx

# PyPy virtual environment recommended for performance

def reducetopologies(graph):
    """
    Create a DAG of reduced topologies (complexity atlas, like Ye et al. 2019) containing at least one Type I or II motif.

    Each node in the resulting DiGraph is a frozenset of the edges that were eliminated from the input DiGraph to produce the reduced topology.
    The "topology" attribute stores the reduced DiGraph; the "motifs" attribute is a tuple of Booleans: has Type I, has Type II.
    Each edge in the atlas has a "removed_edge" attribute indicating which one edge was removed from the input graph.
    This process has exponential complexity in the number of edges, so it is only suitable for very small networks.
    """
    graph = nx.convert_node_labels_to_integers(graph)
    cycle_sets = preparecyclesets(graph)
    atlas = nx.DiGraph()
    root_motifs = getmotifs(cycle_sets, len(graph.nodes))
    if (root_motifs == (False, False)):
        return atlas
    no_removals = frozenset()
    atlas.add_node(no_removals, topology=graph, motifs=root_motifs)
    tryreduce(atlas, graph, set(), no_removals, cycle_sets, graph)
    return atlas

def summarizenetwork(graph):
    """Count minimal Type I and Type II subnetworks of the input DiGraph using topology reduction, like Ye et al. 2019."""
    simple_cycles = len([cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if ispositive(graph, cycle)])
    atlas = reducetopologies(graph)
    type1 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (True, False) and len(atlas.out_edges(n)) == 0])
    type2 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (False, True) and len(atlas.out_edges(n)) == 0])
    return {'pfls': simple_cycles, 'type1': type1, 'type2': type2, 'sum12': type1 + type2}

# All below functions are implementation details

def tryreduce(atlas, original_graph, checked_removals, removed_edges, cycle_sets, graph):
    """Recursively reduce the graph to smaller graphs that still have Type I or Type II motifs, expanding the atlas."""
    for edge in graph.edges:
        new_removed_edges = removed_edges.union(frozenset([edge]))
        if (new_removed_edges in checked_removals):
            if (new_removed_edges in atlas):
                atlas.add_edge(removed_edges, new_removed_edges, removed_edge=edge)
            continue
        checked_removals.add(new_removed_edges)
        reduced_graph = nx.classes.function.restricted_view(original_graph, [], new_removed_edges)
        remaining_cycles = cycle_sets.difference(graph.edges[edge]['cycles'])
        motifs = getmotifs(remaining_cycles, len(graph.nodes))
        if (motifs == (False, False)):
            continue
        atlas.add_node(new_removed_edges, topology=reduced_graph, motifs=motifs)
        atlas.add_edge(removed_edges, new_removed_edges, removed_edge=edge)
        tryreduce(atlas, original_graph, checked_removals, new_removed_edges, remaining_cycles, reduced_graph)

def preparecyclesets(graph):
    """Tag each edge in the input graph with the set of cycles it belongs to, removing the need to recalculate cycles after each reduction."""
    for edge in graph.edges:
        graph.edges[edge]['cycles'] = set()
    all_pfls = [cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if ispositive(graph, cycle)]
    cycle_sets = [IdentityHolder(frozenset(cycle)) for cycle in all_pfls]
    for index, cycle in enumerate(all_pfls):
        for n in range(len(cycle)):
            graph.edges[cycle[n], cycle[(n + 1) % len(cycle)]]['cycles'].add(cycle_sets[index])
    return frozenset(cycle_sets)

def getmotifs(cycle_sets, node_count):
    """Return a tuple of Booleans indicating whether a graph with the given cycle sets has any Type I or Type II (respectively) motifs."""
    return (hastype1(cycle_sets, node_count), hastype2(cycle_sets))

def hastype1(cycle_sets, node_count):
    """Determine whether any three cycles share at least one node (Type I motif)."""
    node_uses = [0] * node_count
    for cycle_holder in cycle_sets:
        for n in cycle_holder.value:
            node_uses[n] += 1
            if (node_uses[n] == 3):
                return True
    return False

def hastype2(cycle_sets):
    """Determine whether there are any pairs of independent cycles bridged by a third (Type II motif)."""
    cycle_sets_list = list(cycle_sets)
    for connector_c in range(len(cycle_sets_list)):
        connector = cycle_sets_list[connector_c]
        for c1 in range(len(cycle_sets_list)):
            if (c1 == connector_c):
                continue
            cycle1 = cycle_sets_list[c1]
            if (connector.value.isdisjoint(cycle1.value)):
                continue
            for c2 in range(c1 + 1, len(cycle_sets_list)):
                if (c2 == connector_c):
                    continue
                cycle2 = cycle_sets_list[c2]
                if ((not connector.value.isdisjoint(cycle2.value)) and cycle1.value.isdisjoint(cycle2.value)):
                    return True
    return False

# Utility functions used in other scripts

def ispositive(graph, cycle):
    """Determine whether the cycle (list of node IDs) is net-positive in the given DiGraph."""
    positive = True
    for i in range(len(cycle)):
        positive ^= graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]['repress']
    return positive
    
def ismutualinhibition(graph, connector_cycle, cycle1, cycle2):
    """
    Determine whether the Type II motif represented by the connection of cycle1 and cycle2 by connector_cycle is a mutual inhibition.

    Arguments:
    - graph: NetworkX DiGraph the cycles were obtained from
    - connector_cycle: list of node IDs for the cycle that connects the two independent cycles
    - cycle1 and cycle2: list or set of node IDs for the independent cycles
    """
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
