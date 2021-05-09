from collections import deque
import networkx as nx

# PyPy virtual environment recommended for performance

def generatecycles(graph, max_length=None):
    """
    Generate simple cycles, in increasing order of length, from the directed graph using Liu & Wang's algorithm.

    Arguments:
    - graph: NetworkX DiGraph
    - max_length: maximum length of cycle to generate (note that Johnson's algorithm is faster for unbounded length)
    """
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
    """Return the best generator for simple cycles in a DiGraph based on the desired maximum length."""
    if max_length:
        # Liu & Wang's algorithm must be used to avoid enumerating excessively long cycles
        return generatecycles(graph, max_length)
    else:
        # Johnson's algorithm is faster for enumerating all cycles
        return nx.simple_cycles(graph)
