import networkx as nx
import random

def permutenetwork(graph):
    '''Randomly moves regulations, preserving degree sequence.'''
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
    '''Attempts a random number of double edge swaps in place.'''
    for _ in range(int(len(graph.edges) * (1 + 2 * random.random()))):
        edges = list(graph.edges)
        ab, cd = random.sample(edges, 2)
        a, b = ab
        c, d = cd
        if graph.has_edge(a, d) or graph.has_edge(c, b):
            continue
        ab_data = graph.edges[ab]
        cd_data = graph.edges[cd]
        graph.remove_edge(a, b)
        graph.remove_edge(c, d)
        graph.add_edge(a, d, **ab_data)
        graph.add_edge(c, b, **cd_data)
    return graph
        
def permuteregulations(graph):
    '''Randomly changes which regulations are repressions, maintaining activation and repression counts and directions.'''
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
