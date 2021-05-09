import argparse
from collections import Counter, defaultdict
from identityholder import IdentityHolder
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition
import networkx as nx
import pandas as pd

# PyPy virtual environment recommended for performance

def countmotifs(network, max_cycle_length=None, max_motif_size=None, check_nfl=False, callback=None):
    """
    Systematically count/enumerate instances of high-feedback topologies in a network.
    
    Each instance is a way to choose a set of edges such that no edge can be removed without abolishing the high-feedback nature of the subnetwork.
    This is similar to the count of "minimum topologies" from Ye et al. 2019, but some topologies with more than the required number of cycles can
    still be considered minimal by this function if removing any edge would leave fewer cycles than necessary. 

    Arguments:
    - network: NetworkX DiGraph representing the network
    - max_cycle_length: maximum length of cycle to enumerate (needed for moderately large networks since all cycles must be held in memory)
    - max_motif_size: maximum size in nodes of counted motif instances
    - check_nfl: whether to count excitable and mixed-sign high-feedback motifs, which require enumerating negative feedback loops as well (much slower)
    - callback: callable to receive progress information and motif instances

    The callback, if provided, must take two positional arguments: notification type and data. These notifications will be given:
    - stage: the counting process moved into a new stage, name specified in data as string
    - instance: a motif instance was found; data is a tuple of the motif name, frozenset of nodes involved, and list of cycle holders
    - cycle_count: the number of cycles was determined, specified in data as int
    - overlap_count: the number of cycle intersections was determined, specified in data as int
    - overlap_progress: progress in enumerating triangles in the cycle intersection graph (e.g. Type 1), data is number of cycle intersections checked
    - cycle_progress: progress in enumerating 3-paths in the cycle intersection graph (e.g. Type 2), data is number of cycles checked

    Returns a tuple: counts of PFLs, Type 1, Type 2, MISA, MISSA, mini-MISSA, mixed-sign high-feedback, excitable.
    """
    if callback is not None:
        callback('stage', 'Finding cycles')
    cycle_sets = [IdentityHolder(frozenset(cycle), (cycle, ispositive(network, cycle), hasrepression(network, cycle))) for
                  cycle in liuwangcycles.generatecycles(network, max_cycle_length) if check_nfl or ispositive(network, cycle)]
    cycle_edge_sets = dict()
    for holder in cycle_sets:
        cycle = holder.tag[0]
        edge_set = frozenset((cycle[n], cycle[(n + 1) % len(cycle)]) for n in range(len(cycle)))
        cycle_edge_sets[holder] = edge_set
        if callback is not None:
            callback('instance', ('Cycle', holder.value, [holder]))
            if holder.tag[1]:
                callback('instance', ('PFL', holder.value, [holder]))
    if callback is not None:
        callback('stage', 'Creating cycle intersection graphs')
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
    if callback is not None:
        callback('cycle_count', len(cycle_sets))
        callback('overlap_count', len(cycle_graph.edges))
        callback('stage', 'Searching for Type I and mixed-sign high-feedback motifs')
    type1 = 0
    mixed = 0
    for triplet in findtriangles(cycle_graph, callback):
        a, b, c = triplet
        possible_type1 = (not check_nfl) or (a.tag[1] and b.tag[1] and c.tag[1])
        if check_nfl and not (a.tag[1] or b.tag[1] or c.tag[1]):
            continue
        shared_nodes = cycle_graph.edges[a, b]['shared'].intersection(cycle_graph.edges[b, c]['shared'])
        if len(shared_nodes) == 0:
            continue
        if not shared_nodes.isdisjoint(cycle_graph.edges[a, c]['shared']):
            all_nodes = a.value.union(b.value).union(c.value) # Not a performance problem - big networks will have a max motif size set anyway
            if max_motif_size and len(all_nodes) > max_motif_size:
                continue
            extra_cycles = list({*findinducedcycles([a, b], c), *findinducedcycles([a, c], b), *findinducedcycles([b, c], a), *findinducedcycles([a, b, c])})
            relevant_extras = [c for c in extra_cycles if c.tag[1]] if possible_type1 else extra_cycles
            if len(relevant_extras) > 0:
                double_counting = False
                if possible_type1:
                    double_counting = any(extra.isbefore(c) for extra in relevant_extras)
                else:
                    rare_sign = a.tag[1] ^ b.tag[1] ^ c.tag[1] # True if there's only one PFL
                    for holder in triplet:
                        if any(extra.isbefore(holder) and (extra.tag[1] == rare_sign or holder.tag[1] != rare_sign) for extra in extra_cycles):
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
                    still_mixed = pfls_after_elimination > 0 and nfls_after_elimination > 0
                    if (possible_type1 or still_mixed) and len(node_uses_after_elimination) > 0 and max(node_uses_after_elimination.values()) >= 3:
                        found_extra = True
                        break
                if found_extra:
                    continue
            if possible_type1:
                type1 += 1
                if callback is not None:
                    callback('instance', ('Type1', all_nodes, triplet))
            else:
                mixed += 1
                if callback is not None:
                    callback('instance', ('MixHF', all_nodes, triplet))
    if callback is not None:
        callback('stage', 'Checking fused pairs')
    missa = 0
    minimissa = 0
    excitable = 0
    for pair in cycle_graph.edges:
        holder1, holder2 = pair
        all_nodes = holder1.value.union(holder2.value)
        if max_motif_size and len(all_nodes) > max_motif_size:
                continue
        if holder1.tag[1] and holder2.tag[1] and (holder1.tag[2] != holder2.tag[2]):
            missa += 1
            if callback is not None:
                callback('instance', ('MISSA', all_nodes, pair))
            if len(holder1.value) == 1 or len(holder2.value) == 1:
                minimissa += 1
                if callback is not None:
                    callback('instance', ('uMISSA', all_nodes, pair))
        elif check_nfl and holder1.tag[1] != holder2.tag[1]:
            excitable += 1
            if callback is not None:
                callback('instance', ('Excite', all_nodes, pair))
    if callback is not None:
        callback('stage', 'Searching for Type II and MISA motifs')
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
                all_nodes = holder.value.union(neigh1.value).union(neigh2.value)
                if max_motif_size and len(all_nodes) > max_motif_size:
                    continue
                if cycle_graph.has_edge(neigh1, neigh2):
                    continue
                if coverextrapfl(holder, neigh2):
                    continue
                type2 += 1
                triplet = (neigh1, holder, neigh2)
                if callback is not None:
                    callback('instance', ('Type2', all_nodes, triplet))
                if ismutualinhibition(network, holder.tag[0], neigh1.value, neigh2.value):
                    misa += 1
                    if callback is not None:
                        callback('instance', ('MISA', all_nodes, triplet))
        if callback is not None:
            checked += 1
            callback('cycle_progress', checked)
    pfls = sum(1 for holder in cycle_sets if holder.tag[1])
    return (pfls, type1, type2, misa, missa, minimissa, mixed, excitable)

def hasrepression(graph, cycle):
    """Return whether the cycle (list of node IDs) in the graph (NetworkX DiGraph) includes any repressions."""
    return any(graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]['repress'] for i in range(len(cycle)))

def findtriangles(graph, callback):
    """
    Generate triangles (tuples of cycle holders) of the cycle intersection graph for countmotifs.
    To avoid duplicates, the third node C is always after A and B in the IdentityHolder total order.
    """
    checked = 0
    for a, b in graph.edges:
        for c in frozenset(graph[a]).intersection(frozenset(graph[b])):
            if a.isbefore(c) and b.isbefore(c):
                yield (a, b, c)
        if callback is not None:
            checked += 1
            if checked % 100 == 0:
                callback('overlap_progress', checked)

def countmotifspernode(callback, *args, **kwargs):
    """
    Systematically count motif instances and how many of each motif each node is involved in.

    Wraps countmotifs. The callback, if specified, will not receive "instance" notifications.
    Returns a tuple: countmotifs results tuple, dict of dicts {motif: {node: instances}}.
    """
    motif_involvement = defaultdict(Counter)
    def counting_callback(notification, data):
        if notification == 'instance':
            motif, node_ids, _ = data
            motif_involvement[motif].update(node_ids)
        elif callback is not None:
            callback(notification, data)
    counts = countmotifs(*args, **kwargs, callback=counting_callback)
    named_motif_involvement = {motif: {graph.nodes[node]['name']: counts for node, counts in motif_counts.items()} for motif, motif_counts in motif_involvement.items()}
    return counts, named_motif_involvement

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help='GraphML file to process')
    parser.add_argument('--maxcycle', type=int, help='maximum cycle length')
    parser.add_argument('--maxnodes', type=int, help='maximum number of nodes in a motif')
    parser.add_argument('--checknfl', action='store_true', help='also search for motifs involving negative feedback')
    parser.add_argument('--nodecounts', nargs='?', const=1, type=str, help='count how many motifs each node is in (optional CSV output)')
    parser.add_argument('-q', '--quiet', action='store_true', help='do not print progress')
    args = parser.parse_args()
    graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
    cycle_count, node_sharing_count = 0, 0
    def progress_callback(notification, data):
        global cycle_count, node_sharing_count
        if notification == 'stage':
            print(data)
        elif notification == 'cycle_count':
            cycle_count = data
            print(f'{cycle_count} cycles, ', end='')
        elif notification == 'overlap_count':
            node_sharing_count = data
            print(f'{node_sharing_count} node sharings')
        elif notification == 'overlap_progress':
            print(f'{data}/{node_sharing_count}\r', end='')
        elif notification == 'cycle_progress':
            print(f'{data}/{cycle_count}\r', end='')
    callback = None if args.quiet else progress_callback
    if args.nodecounts is not None:
        result, motif_involvement = countmotifspernode(callback, graph, args.maxcycle, args.maxnodes, args.checknfl)
    else:
        result = countmotifs(graph, args.maxcycle, args.maxnodes, args.checknfl, callback)
    print(''.ljust(20, ' '), '\r', sep='', end='')
    print('PFL', result[0], '\nType1', result[1], '\nType2', result[2], '\nMISA', result[3], '\nMISSA', result[4], '\nuMISSA', result[5], sep='\t')
    if args.checknfl:
        print('MixHF', result[6], '\nExcite', result[7], sep='\t')
    if args.nodecounts is not None:
        df = pd.DataFrame.from_dict(motif_involvement).fillna(0).astype(int)
        if isinstance(args.nodecounts, str):
            df.to_csv(args.nodecounts)
        else:
            print('\n', df, sep='')
    