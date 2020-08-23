import minimumtopologies
import networkx as nx

def summarizenetwork(graph):
    simple_cycles = len([cycle for cycle in nx.algorithms.cycles.simple_cycles(graph) if minimumtopologies.ispositive(graph, cycle)])
    atlas = minimumtopologies.reducetopologies(graph)
    type1 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (True, False) and len(atlas.out_edges(n)) == 0])
    type2 = len([n for n in atlas.nodes if atlas.nodes[n]['motifs'] == (False, True) and len(atlas.out_edges(n)) == 0])
    return {'pfls': simple_cycles, 'type1': type1, 'type2': type2, 'sum12': type1 + type2}
