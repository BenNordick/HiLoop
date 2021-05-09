import networkx as nx
import sys

def loadadjmatrix(array, names=None):
    """
    Create a DiGraph from a signed adjacency matrix.

    Each value in the array specifies the interaction between the row and column nodes.
    Zero means no edge, a positive number means activation, and a negative number means repression.

    Arguments:
    - array: list of lists, first index indicates source, second index indicates target
    - names: optional list of names for the nodes (same order as array)
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(range(len(array)))
    if names is None:
        for n in graph.nodes:
            graph.nodes[n]['name'] = str(n)
    else:
        if len(names) != len(array):
            raise ValueError("there must be exactly one name per gene (matrix row)")
        for n, name in enumerate(names):
            graph.nodes[n]['name'] = name
    for src, row in enumerate(array):
        for dst, interaction in enumerate(row):
            if interaction == 0:
                continue
            graph.add_edge(src, dst, repress=(interaction < 0))
    return graph

def loadreader(matrix, names=None):
    """
    Create a DiGraph from a reader of a tab-separated adjacency matrix text document.

    For each entry, the row specifies the source and the column specifies the target.
    The main matrix document must contain only the body of the matrix, not node names, which can be provided in the names reader.
    Comment lines are indicated with percent signs. Text on a name line after an equals sign is also ignored.
    """
    array = []
    for line in matrix:
        array.append([int(part) for part in line.strip().split('\t')])
    names_array = None
    if names is not None:
        names_array = []
        for line in names:
            if line.startswith('%') or line.strip() == '':
                continue
            names_array.append(line.split('=')[0].strip())
    return loadadjmatrix(array, names_array)

if __name__ == "__main__":
    with open(sys.argv[2]) as matrix:
        if len(sys.argv) > 3:
            with open(sys.argv[3]) as names:
                graph = loadreader(matrix, names)
        else:
            graph = loadreader(matrix)
    nx.write_graphml(graph, sys.argv[1])
