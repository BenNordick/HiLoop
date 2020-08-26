import networkx as nx
import sys

def loadadjmatrix(array, names=None):
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
    array = []
    for line in matrix:
        array.append([int(part) for part in line.strip().split('\t')])
    names_array = None
    if names is not None:
        names_array = []
        for line in names:
            if line.startswith('%') or line.strip() == '':
                continue
            if '=' not in line:
                raise ValueError('missing equals sign in non-comment line ' + repr(line))
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
