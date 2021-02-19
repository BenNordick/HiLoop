import networkx as nx
from networkx.drawing.nx_agraph import to_agraph

def rendergraph(graph, filename, in_place=False):
    if not in_place:
        graph = graph.copy()
    for n in graph.nodes:
        graph.nodes[n]['label'] = graph.nodes[n]['name']
    for e in graph.edges:
        edge = graph.edges[e]
        if (edge['repress']):
            edge['arrowhead'] = 'tee'
    ag = to_agraph(graph)
    ag.layout('dot')
    ag.draw(filename)

def colorcycle(graph, cycle, color, mix_color):
    for i in range(len(cycle)):
        edge = graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]
        if 'color' in edge:
            if isinstance(mix_color, dict):
                edge['color'] = mix_color[edge['color']]
            else:
                edge['color'] = mix_color
        else:
            edge['color'] = color

def colorcycles(graph, cycles):
    multigraph = nx.MultiDiGraph(graph)
    colored_edges = set()
    for cycle_info in cycles:
        if len(cycle_info) < 3:
            cycle, color = cycle_info
            style = 'solid'
        else:
            cycle, color, style = cycle_info
        for i in range(len(cycle)):
            src = cycle[i]
            dst = cycle[(i + 1) % len(cycle)]
            if (src, dst) in colored_edges:
                route = multigraph.add_edge(src, dst, **graph.edges[src, dst])
                multigraph.edges[src, dst, route]['color'] = color
                multigraph.edges[src, dst, route]['style'] = style
            else:
                multigraph.edges[src, dst, 0]['color'] = color
                multigraph.edges[src, dst, 0]['style'] = style
                colored_edges.add((src, dst))
    return multigraph

def colorneighborhood(graph, start_node, colors, color_last_edges=True):
    current_ring = set([start_node])
    next_ring = set()
    level = 0
    while level < len(colors):
        for n in current_ring:
            graph.nodes[n]['color'] = colors[level]
            if color_last_edges or level < len(colors) - 1:
                for neighbor in graph.successors(n):
                    graph.edges[n, neighbor]['color'] = colors[level]
                    if 'color' not in graph.nodes[neighbor] and neighbor not in current_ring:
                        next_ring.add(neighbor)
        current_ring = next_ring
        next_ring = set()
        level += 1
