from networkx.drawing.nx_agraph import to_agraph

def rendergraph(graph, filename):
    graph = graph.copy()
    for n in graph.nodes:
        graph.nodes[n]['label'] = graph.nodes[n]['name']
    for e in graph.edges:
        edge = graph.edges[e[0], e[1]]
        if (edge['repress']):
            edge['arrowhead'] = 'tee'
    ag = to_agraph(graph)
    ag.layout('dot')
    ag.draw(filename)

def colorcycle(graph, cycle, color, mix_color):
    for i in range(len(cycle)):
        edge = graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]
        edge['color'] = mix_color if ('color' in edge) else color
