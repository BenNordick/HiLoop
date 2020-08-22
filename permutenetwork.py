import networkx as nx
import random

def permutenetwork(graph, swaps):
    if len(graph.edges) < 2:
        raise ValueError("graph must have at least two edges")
    if not nx.algorithms.components.is_strongly_connected(graph):
        raise ValueError("graph must be strongly connected")
    while True:
        copy = graph.copy()
        for _ in range(swaps):
            edgeswap(graph)
        if nx.algorithms.components.is_strongly_connected(copy):
            return copy

def edgeswap(graph):
    if isinstance(graph, list):
        for g in graph:
            edgeswap(g)
        return
    while True:
        if tryedgeswap(graph):
            break

def tryedgeswap(graph):
    edges = list(graph.edges)
    ab, cd = random.choices(edges, k=2)
    if (ab == cd):
        return False
    a, b = ab
    c, d = cd
    if graph.has_edge(a, d) or graph.has_edge(c, b):
        return False
    ab_data = graph.edges[ab]
    cd_data = graph.edges[cd]
    graph.remove_edge(a, b)
    graph.remove_edge(c, d)
    createedge(graph, a, d, ab_data)
    createedge(graph, c, b, cd_data)
    return True
    
def createedge(graph, src, dst, data):
    graph.add_edge(src, dst)
    for k, v in data.items():
        graph.edges[src, dst][k] = v
