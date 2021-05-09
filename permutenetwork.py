from collections import deque
import networkx as nx
import random
import re

# Support functions for making random permutations to regulatory networks
# PyPy virtual environment recommended for performance

def permutenetwork(graph):
    """
    Randomly reposition regulations by rejection sampling of directed configuration models.

    Because directed_configuration_model often produces multigraphs which must be rejected, this process is only suitable for small networks.
    It does not change the original graph.
    """
    repression_count = len([edge for edge in graph.edges if graph.edges[edge]['repress']])
    in_degrees = [graph.in_degree(n) for n in graph.nodes]
    out_degrees = [graph.out_degree(n) for n in graph.nodes]
    while True:
        multigraph = nx.generators.degree_seq.directed_configuration_model(in_degrees, out_degrees)
        flattened = nx.DiGraph(multigraph)
        if len(flattened.edges) != len(graph.edges):
            continue
        for edge in flattened.edges:
            flattened.edges[edge]['repress'] = False
        for edge in random.sample(list(flattened.edges), repression_count):
            flattened.edges[edge]['repress'] = True
        for index, node in enumerate(graph.nodes):
            flattened.nodes[index]['name'] = graph.nodes[node]['name']
        return flattened

def permuteedgeswaps(graph):
    """Attempt a random number of double edge swaps and directed triangle reversals in place."""
    for _ in range(int(len(graph.edges) * (1 + 2 * random.random()))):
        edges = list(graph.edges)
        ab, cd = random.sample(edges, 2)
        a, b = ab
        c, d = cd
        if graph.has_edge(a, d) or graph.has_edge(c, b):
            a_pred = list(graph.predecessors(a))
            if len(a_pred) == 0:
                continue
            e = random.choice(a_pred)
            if e != b and graph.has_edge(b, e):
                reversedirectedtriangle(graph, a, b, e)
            continue
        ab_data = graph.edges[ab]
        cd_data = graph.edges[cd]
        graph.remove_edge(a, b)
        graph.remove_edge(c, d)
        graph.add_edge(a, d, **ab_data)
        graph.add_edge(c, b, **cd_data)
    return graph

def reversedirectedtriangle(graph, a, b, c):
    """Try to reverse a directed triangle a -> b -> c -> a in place (if no other existing edges conflict)."""
    legs = [(a, b), (b, c), (c, a)]
    for src, dest in legs:
        if graph.has_edge(dest, src):
            return
    for src, dest in legs:
        data = graph.edges[src, dest]
        graph.remove_edge(src, dest)
        graph.add_edge(dest, src, **data)
        
def permuteregulations(graph):
    """Randomly change which regulations are repressions, maintaining activation and repression counts and directions."""
    edges = list(graph.edges)
    copy = graph.copy()
    repressions = 0
    for edge in edges:
        edge_data = copy.edges[edge]
        if edge_data['repress']:
            repressions += 1
            edge_data['repress'] = False
    for new_repression in random.sample(edges, repressions):
        copy.edges[new_repression]['repress'] = True
    return copy

def restorefixedsigns(original, graph, pattern):
    """
    Swap edge signs in place to ensure that the specified nodes' outgoing edge signs are the same as in the original graph.

    Arguments:
    - original: the original network as a NetworkX DiGraph
    - graph: the permuted network (to be corrected) as a NetworkX DiGraph
    - pattern: regular expression to match node names whose outgoing edge signs should not have been changed (e.g. "miR")
    """
    fixed_sign_source_nodes = set(n for n in graph.nodes if re.match(pattern, graph.nodes[n]['name']) is not None)
    available_repressions = [e for e in graph.edges if graph.edges[e]['repress'] and e[0] not in fixed_sign_source_nodes]
    available_activations = [e for e in graph.edges if (not graph.edges[e]['repress']) and e[0] not in fixed_sign_source_nodes]
    for src in fixed_sign_source_nodes:
        if graph.nodes[src]['name'] != original.nodes[src]['name']:
            raise AssertionError('ID mismatch')
        current_repressions = [(src, n) for n in graph[src] if graph.edges[src, n]['repress']]
        current_activations = [(src, n) for n in graph[src] if not graph.edges[src, n]['repress']]
        n_original_repressions = len([n for n in original[src] if original.edges[src, n]['repress']])
        while len(current_repressions) > n_original_repressions:
            extra_repression = random.choice(current_repressions)
            old_activation = random.choice(available_activations)
            graph.edges[extra_repression]['repress'] = False
            graph.edges[old_activation]['repress'] = True
            current_repressions.remove(extra_repression)
            current_activations.append(extra_repression)
            available_activations.remove(old_activation)
            available_repressions.append(old_activation)
        while len(current_repressions) < n_original_repressions:
            extra_activation = random.choice(current_activations)
            old_repression = random.choice(available_repressions)
            graph.edges[extra_activation]['repress'] = True
            graph.edges[old_repression]['repress'] = False
            current_activations.remove(extra_activation)
            current_repressions.append(extra_activation)
            available_repressions.remove(old_repression)
            available_activations.append(old_repression)

def randomsubgraph(graph, max_nodes):
    """Select a random connected subgraph with the specified number of nodes."""
    queue = deque(maxlen=len(graph.nodes))
    selected = set()
    queue.append(random.sample(graph.nodes, 1)[0])
    while len(selected) < max_nodes and len(queue) > 0:
        head = queue.popleft()
        if head in selected:
            continue
        selected.add(head)
        neighbors = list(graph.successors(head))
        random.shuffle(neighbors)
        for n in neighbors:
            queue.append(n)
    return graph.subgraph(selected)

def neighborhoodsubgraph(graph, start_node, depth):
    """
    Select a subgraph containing all nodes at most the specified depth away from a starting node ID.

    Likely only useful for debugging.
    """
    current_ring = set([start_node])
    next_ring = set()
    selected = set(current_ring)
    level = 0
    while level < depth:
        for n in current_ring:
            for neighbor in graph.successors(n):
                if neighbor not in selected:
                    selected.add(neighbor)
                    next_ring.add(neighbor)
        current_ring = next_ring
        next_ring = set()
        level += 1
    return graph.subgraph(selected)

def generatepermutations(graph, require_connected, use_full_permutation=True, max_nodes_for_sample=None, fixed_sign_sources=None):
    """
    Generate permutations of an input graph.

    Arguments:
    - graph: NetworkX DiGraph to generate permutations of
    - require_connected: whether to produce only strongly connected permutations
    - use_full_permutation: whether to use rejection sampling of directed configuration models
    - max_nodes_for_sample: produce subgraphs of this size (not recommended)
    - fixed_sign_sources: regex to match names of nodes whose outgoing edge signs should remain constant

    Yields a tuple of permutation count and permuted DiGraph.
    """
    last_permutation = graph
    checked_permutations = 0
    while True:
        if checked_permutations % 20 == 0 and use_full_permutation:
            last_permutation = permutenetwork(graph)
        else:
            last_permutation = permuteedgeswaps(permuteregulations(last_permutation))
        if require_connected and not nx.algorithms.is_strongly_connected(last_permutation):
            continue
        if fixed_sign_sources is not None:
            restorefixedsigns(graph, last_permutation, fixed_sign_sources)
        feasible_subgraph = last_permutation if max_nodes_for_sample is None else randomsubgraph(last_permutation, max_nodes_for_sample)
        checked_permutations += 1
        yield checked_permutations, feasible_subgraph
