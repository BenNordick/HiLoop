from collections import deque
import networkx as nx

def generatecycles(graph, max_length=None):
    '''Generates simple cycles in the directed graph using Liu & Wang's algorithm.'''
    queue = deque([n] for n in graph.nodes)
    while queue:
        path = queue.popleft()
        if graph.has_edge(path[-1], path[0]):
            yield path
        if max_length is not None and len(path) >= max_length:
            continue
        for n in set(graph.successors(path[-1])).difference(path):
            if n > path[0]:
                queue.append(path + [n])

def cyclesgenerator(graph, max_length=None):
    '''Returns the best generator for simple cycles based on the desired maximum length.'''
    if max_length:
        return generatecycles(graph, max_length)
    else:
        return nx.simple_cycles(graph)
