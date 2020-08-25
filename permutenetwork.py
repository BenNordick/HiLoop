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
        