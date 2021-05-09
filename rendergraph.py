import networkx as nx
from networkx.drawing.nx_agraph import to_agraph

# Support functions for drawing network diagrams

def rendergraph(graph, filename, in_place=False):
    """
    Draw a graph to an image file.

    Arguments:
    - graph: NetworkX DiGraph
    - filename: where to save the image (format auto-detected by Graphviz)
    - in_place: whether to make a copy of the graph before adding Graphviz attributes
    """
    ag = graphvizify(graph, in_place=in_place)
    ag.draw(filename)

def graphvizify(graph, in_place=False, layout='dot'):
    """
    Turn a NetworkX DiGraph into a PyGraphviz graph with drawing attributes.

    Arguments:
    - graph: NetworkX DiGraph
    - in_place: whether to make a copy before adding attributes
    - layout: Graphviz layout engine to use for node positions (or None to delay layout)

    Returns a PyGraphviz AGraph.
    """
    if not in_place:
        graph = graph.copy()
    for n in graph.nodes:
        graph.nodes[n]['label'] = graph.nodes[n]['name']
        graph.nodes[n]['fontname'] = 'DejaVuSerif'
    for e in graph.edges:
        edge = graph.edges[e]
        if (edge['repress']):
            edge['arrowhead'] = 'tee'
    ag = to_agraph(graph)
    if layout is not None:
        ag.layout(layout)
    return ag

def highlightedge(edge_attrs):
    """Set edge attributes on a dict to enlarge the edge."""
    edge_attrs['penwidth'] = 2.0
    edge_attrs['arrowsize'] = 1.2

def colorcycle(graph, cycle, color, mix_color):
    """
    Color and highlight a cycle of a DiGraph in-place.

    Arguments:
    - graph: NetworkX DiGraph
    - cycle: list of node IDs in the cycle
    - color: name of color to apply to yet-uncolored edges in the cycle
    - mix_color: name of color to apply to already colored edges, or dict of existing color to new color
    """
    for i in range(len(cycle)):
        edge = graph.edges[cycle[i], cycle[(i + 1) % len(cycle)]]
        highlightedge(edge)
        if 'color' in edge:
            if isinstance(mix_color, dict):
                edge['color'] = mix_color[edge['color']]
            else:
                edge['color'] = mix_color
        else:
            edge['color'] = color

def colorcycles(graph, cycles):
    """
    Color and highlight cycles of a DiGraph without color mixing, producing a multigraph.

    Arguments:
    - graph: NetworkX DiGraph
    - cycles: list of lists of node IDs in each cycle

    Returns a NetworkX MultiDiGraph.
    """
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
            else:
                route = 0
                colored_edges.add((src, dst))
            multigraph.edges[src, dst, route]['color'] = color
            multigraph.edges[src, dst, route]['style'] = style
            highlightedge(multigraph.edges[src, dst, route])
    return multigraph

def colorneighborhood(graph, start_node, colors, color_last_edges=True):
    """
    Color nodes and edges according to how far they are from a specific starting node.

    Likely only useful for debugging. Arguments:
    - graph: NetworkX DiGraph to color in-place
    - start_node: node ID of starting node
    - colors: list of colors for nodes and outgoing edges at each level (starting node first)
    - color_last_edges: whether to color outgoing edges of the last level of colored nodes
    """
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
