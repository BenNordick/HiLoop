from identityholder import IdentityHolder
import networkx as nx

def reducetopologies(graph):
    '''Produces a DAG of reduced topologies containing at least one Type I or II motif.'''
    graph = nx.convert_node_labels_to_integers(graph)
    cycle_sets = preparecyclesets(graph)
    atlas = nx.DiGraph()
    root_motifs = getmotifs(cycle_sets, len(graph.nodes))
    if (root_motifs == (False, False)):
        return atlas
    no_nodes_removed = frozenset()
    atlas.add_node(no_nodes_removed, topology=graph, motifs=root_motifs)
    tryreduce(atlas, graph, set(), no_nodes_removed, cycle_sets, graph)
    return atlas

def summarizenetwork(graph):
    simple_cycles = len([cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if ispositive(graph, cycle)])
    atlas = reducetopologies(graph)
    type1 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (True, False) and len(atlas.out_edges(n)) == 0])
    type2 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (False, True) and len(atlas.out_edges(n)) == 0])
    return {'pfls': simple_cycles, 'type1': type1, 'type2': type2, 'sum12': type1 + type2}

# All below functions are implementation details

def tryreduce(atlas, original_graph, checked_removals, removed_edges, cycle_sets, graph):
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
    for edge in graph.edges:
        graph.edges[edge]['cycles'] = set()
    all_pfls = [cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if ispositive(graph, cycle)]
    cycle_sets = [IdentityHolder(frozenset(cycle)) for cycle in all_pfls]
    for index, cycle in enumerate(all_pfls):
        for n in range(len(cycle)):
            graph.edges[cycle[n], cycle[(n + 1) % len(cycle)]]['cycles'].add(cycle_sets[index])
    return frozenset(cycle_sets)

def getmotifs(cycle_sets, node_count):
    return (hastype1(cycle_sets, node_count), hastype2(cycle_sets))

def hastype1(cycle_sets, node_count):
    node_uses = [0] * node_count
    for cycle_holder in cycle_sets:
        for n in cycle_holder.value:
            node_uses[n] += 1
            if (node_uses[n] == 3):
                return True
    return False

def hastype2(cycle_sets):
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
