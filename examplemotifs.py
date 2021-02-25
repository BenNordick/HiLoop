import argparse
import io
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition, ishalfinhibition
import networkx as nx
import permutenetwork
from PIL import Image
import pygraphviz
import random
import rendergraph

def colorsubgraph(graph, r, y, b):
    def cyclestyle(graph, cycle):
        if len(cycle) == 0:
            return 'solid'
        return 'solid' if ispositive(graph, cycle) else 'dashed'
    return rendergraph.colorcycles(graph, [(r, 'red', cyclestyle(graph, r)), (y, 'gold', cyclestyle(graph, y)), (b, 'blue', cyclestyle(graph, b))])

def logobase(**kwargs):
    ag = pygraphviz.AGraph(bgcolor='#D0D0D0', strict=False, directed=True, ranksep=0.3, **kwargs)
    ag.edge_attr['arrowsize'] = 0.8
    return ag

def logoimage(ag):
    ag.layout('dot')
    png_bytes = ag.draw(format='png')
    return Image.open(io.BytesIO(png_bytes))

def logo_fpnp(fusion_nodes):
    ag = logobase()
    ag.add_node('X', label=',\n'.join(fusion_nodes), fontsize=9.0, width=0.1, height=0.1, margin=0.05)
    ag.add_edge('X', 'X', 0, color='red', style='dashed', arrowhead='tee', headport='ne', tailport='se')
    ag.add_edge('X', 'X', 1, color='blue', headport='sw', tailport='nw')
    return logoimage(ag)

def logo_3fused(fusion_nodes, positivities):
    styles = [('solid' if positive else 'dashed') for positive in positivities]
    tips = [('normal' if positive else 'tee') for positive in positivities]
    ag = logobase(nodesep=0.6)
    ag.add_node('L3', shape='point', width=0.001, color='blue')
    ag.add_node('X', label=',\n'.join(fusion_nodes), fontsize=9.0, width=0.1, height=0.1, margin=0.05)
    ag.add_edge('L3', 'X', arrowhead=tips[2], color='blue', style=styles[2])
    ag.add_edge('X', 'L3', arrowhead='none', color='blue', style=styles[2])
    ag.add_node('L1', shape='point', width=0.001, color='red')
    ag.add_edge('X', 'L1', arrowhead='none', color='red', style=styles[0])
    ag.add_edge('L1', 'X', arrowhead=tips[0], color='red', style=styles[0])
    ag.add_node('L2', shape='point', width=0.001, color='gold')
    ag.add_edge('X', 'L2', arrowhead='none', color='gold', style=styles[1])
    ag.add_edge('L2', 'X', arrowhead=tips[1], color='gold', style=styles[1])
    return logoimage(ag)

def logo_2bridged(fusion1, fusion2, positivities, half12_positive, half21_positive):
    styles = [('solid' if positive else 'dashed') for positive in positivities]
    tips = [('normal' if positive else 'tee') for positive in positivities]
    ag = logobase()
    ag.add_node('C1', label=',\n'.join(fusion1), fontsize=9.0, width=0.1, height=0.1, margin=0.05)
    ag.add_node('C2', label=',\n'.join(fusion2), fontsize=9.0, width=0.1, height=0.1, margin=0.05)
    ag.add_edge('C1', 'C1', color='red', style=styles[0], arrowhead=tips[0])
    ag.add_edge('C1', 'C2', color='gold', style=styles[1], arrowhead=('normal' if half12_positive else 'tee'))
    ag.add_edge('C2', 'C1', color='gold', style=styles[1], arrowhead=('normal' if half21_positive else 'tee'))
    ag.add_edge('C2', 'C2', color='blue', style=styles[2], arrowhead=tips[2])
    return logoimage(ag)

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='input GraphML file')
parser.add_argument('-1', '--find1', type=int, default=0, help='how many Type I examples to find')
parser.add_argument('-2', '--find2', type=int, default=0, help='how many Type II examples to find')
parser.add_argument('-x', '--findexcitable', type=int, default=0, help='how many excitable examples to find')
parser.add_argument('-m', '--findmisa', type=int, default=0, help='how many MISA Type II examples to find')
parser.add_argument('--findnegative1', type=int, default=0, help='how many negative Type I examples to find')
parser.add_argument('--findnegative2', type=int, default=0, help='how many negative Type II examples to find')
parser.add_argument('-p', '--findfpnp', type=int, default=0, help='how many fused PFL/NFL pair examples to find')
parser.add_argument('--images', type=str, help='output image file pattern {id, type, edges}')
parser.add_argument('--logo', action='store_true', help='add motif summary/logo to saved images')
parser.add_argument('--networks', type=str, help='output GraphML file pattern {id, type}')
parser.add_argument('--printnodes', action='store_true', help='print distinct node names')
parser.add_argument('--maxedges', type=int, help='maximum number of unique edges in examples')
parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
parser.add_argument('--maxcycle', type=int, help='maximum number of nodes in a cycle')
parser.add_argument('--maxsharing', type=int, help='maximum number of nodes in common with an already selected subnetwork')
parser.add_argument('--reduceedges', action='store_true', help='randomly drop some extra edges')
parser.add_argument('--requirenodes', nargs='+', type=str, help='node(s) that must be present in the subnetwork')
args = parser.parse_args()

graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
type1, type2, excitable, misa, negative1, negative2, fpnp = 0, 0, 0, 0, 0, 0, 0
seen = []
seen_edgesets = []
printed_nodes = set()

def shouldcheck3fused():
    return type1 < args.find1 or excitable < args.findexcitable or negative1 < args.findnegative1

def shouldcheckbridged():
    return type2 < args.find2 or misa < args.findmisa or negative2 < args.findnegative2

def printnewnodes(nodes):
    if not args.printnodes:
        return
    for n in nodes.difference(printed_nodes):
        printed_nodes.add(n)
        print(graph.nodes[n]['name'])

def reduceedges(subgraph, cycles):
    if args.reduceedges:
        trimmed = nx.DiGraph()
        trimmed.add_nodes_from(subgraph.nodes)
        for n in trimmed.nodes:
            trimmed.nodes[n]['name'] = subgraph.nodes[n]['name']
        necessary_edges = set()
        for cycle in cycles:
            for i in range(len(cycle)):
                edge = (cycle[i], cycle[(i + 1) % len(cycle)])
                necessary_edges.add(edge)
                trimmed.add_edge(edge[0], edge[1], repress=subgraph.edges[edge]['repress'])
        for e in subgraph.edges:
            if (e not in necessary_edges) and random.uniform(0, 1) < (2 / 3):
                trimmed.add_edge(e[0], e[1], repress=subgraph.edges[e]['repress'])
        return trimmed
    else:
        return subgraph

def createimage(graph, filename, logo_func):
    if args.logo:
        logo = logo_func()
        logo_w, logo_h = logo.size
        main_ag = rendergraph.graphvizify(graph, in_place=True, layout=None)
        main_ag.add_node('Logo', label='', shape='box', style='filled', color='#D0D0D0', width=(logo_w / 96), height=(logo_h / 96))
        main_ag.layout('dot')
        main = Image.open(io.BytesIO(main_ag.draw(format='png')))
        main_w, main_h = main.size
        found_placeholder = False
        for y in range(main_h):
            for x in range(main_w):
                if main.getpixel((x, y))[:3] == (0xD0, 0xD0, 0xD0) == main.getpixel((x + 5, y + 5))[:3]:
                    main.paste(logo, (x, y))
                    found_placeholder = True
                    break
            if found_placeholder:
                break
        if not found_placeholder:
            print('Could not find logo placeholder')
        main.save(filename)
    else:
        rendergraph.rendergraph(graph, filename, in_place=True)

cycles = None
if args.maxnodes is None:
    cycles = [(cycle, ispositive(graph, cycle)) for cycle in liuwangcycles.cyclesgenerator(graph, args.maxcycle)]
    if len(cycles) < 3:
        print('Insufficient cycles for high feedback')

def pickcycles(count):
    if count > len(cycles):
        return None
    chosen_cycles = random.sample(cycles, count)
    cycle_sets = [frozenset(cycle) for cycle, positive in chosen_cycles]
    used_nodes = cycle_sets[0]
    for cs in cycle_sets[1:]:
        used_nodes = used_nodes.union(cs)
    if args.requirenodes is not None:
        used_names = {graph.nodes[n]['name'] for n in used_nodes}
        if not used_names.issuperset(args.requirenodes):
            return None
    if args.maxsharing is not None:
        for ns in seen:
            if len(ns.intersection(used_nodes)) > args.maxsharing:
                return None
    subgraph = reduceedges(graph.subgraph(used_nodes), [cycle for cycle, sign in chosen_cycles])
    if args.maxedges is not None and len(subgraph.edges) > args.maxedges:
        return None
    edges_set = set(subgraph.edges)
    if edges_set in seen_edgesets:
        return None
    seen_edgesets.append(edges_set)
    return chosen_cycles, cycle_sets, used_nodes, subgraph

while shouldcheck3fused() or shouldcheckbridged() or fpnp < args.findfpnp:
    if args.maxnodes is not None:
        feasible = permutenetwork.randomsubgraph(graph, args.maxnodes)
        cycles = [(cycle, ispositive(feasible, cycle)) for cycle in liuwangcycles.cyclesgenerator(feasible, args.maxcycle)]
    if fpnp < args.findfpnp:
        pick_results = pickcycles(2)
        if pick_results is not None:
            chosen_cycles, cycle_sets, used_nodes, subgraph = pick_results
            if chosen_cycles[0][1] != chosen_cycles[1][1] and not cycle_sets[0].isdisjoint(cycle_sets[1]):
                seen.append(used_nodes)
                printnewnodes(used_nodes)
                fpnp += 1
                if args.networks:
                    nx.write_graphml(subgraph, args.networks.format(fpnp, 'fpnp'))
                if args.images:
                    red, blue = (chosen_cycles[1][0], chosen_cycles[0][0]) if chosen_cycles[0][1] else (chosen_cycles[0][0], chosen_cycles[1][0])
                    colored = colorsubgraph(subgraph, red, [], blue)
                    shared_nodes = cycle_sets[0].intersection(cycle_sets[1])
                    for n in shared_nodes:
                        colored.nodes[n]['penwidth'] = 2.0
                    logo_func = lambda: logo_fpnp([subgraph.nodes[n]['name'] for n in shared_nodes])
                    createimage(colored, args.images.format(fpnp, 'fpnp', len(colored.edges)), logo_func)
    if len(cycles) < 3 or not (shouldcheck3fused() or shouldcheckbridged()):
        continue
    pick_results = pickcycles(3)
    if pick_results is None:
        continue
    chosen_cycles, cycle_sets, used_nodes, subgraph = pick_results
    loop_signs = [positive for cycle, positive in chosen_cycles]
    if shouldcheck3fused():
        intersection = cycle_sets[0].intersection(cycle_sets[1]).intersection(cycle_sets[2])
        if len(intersection) > 0:
            kind = None
            current_id = None
            if all(loop_signs):
                if type1 < args.find1:
                    kind = 'type1'
                    type1 += 1
                    current_id = type1
            elif not any(loop_signs):
                if negative1 < args.findnegative1:
                    kind = 'negative1'
                    negative1 += 1
                    current_id = negative1
            else:
                if excitable < args.findexcitable:
                    kind = 'excite'
                    excitable += 1
                    current_id = excitable
            if kind is None:
                continue
            seen.append(used_nodes)
            printnewnodes(used_nodes)
            if args.networks:
                nx.write_graphml(subgraph, args.networks.format(current_id, kind))
            if args.images:
                colored = colorsubgraph(subgraph, *[cycle for cycle, sign in chosen_cycles])
                for n in intersection:
                    colored.nodes[n]['penwidth'] = 2.0
                logo_func = lambda: logo_3fused({subgraph.nodes[n]['name'] for n in intersection}, loop_signs)
                createimage(colored, args.images.format(current_id, kind, len(colored.edges)), logo_func)
            continue
    if shouldcheckbridged():
        for connector_c, connector in enumerate(cycle_sets):
            other1 = cycle_sets[(connector_c + 1) % 3]
            other2 = cycle_sets[(connector_c + 2) % 3]
            if connector.isdisjoint(other1) or connector.isdisjoint(other2):
                continue
            if other1.isdisjoint(other2):
                kind = None
                current_id = None
                if all(loop_signs):
                    if type2 < args.find2 or misa < args.findmisa:
                        kind = 'type2'
                        type2 += 1
                        current_id = type2
                        if misa < args.findmisa and ismutualinhibition(subgraph, chosen_cycles[connector_c][0], other1, other2):
                            kind = 'type2misa'
                            misa += 1
                            current_id = misa
                elif not any(loop_signs):
                    if negative2 < args.findnegative2:
                        kind = 'negative2'
                        negative2 += 1
                        current_id = negative2
                if kind is None:
                    continue
                seen.append(used_nodes)
                printnewnodes(used_nodes)
                if args.networks:
                    nx.write_graphml(subgraph, args.networks.format(current_id, kind))
                if args.images:
                    colored = colorsubgraph(subgraph, chosen_cycles[(connector_c + 1) % 3][0], chosen_cycles[connector_c][0], chosen_cycles[(connector_c + 2) % 3][0])
                    def logo_func():
                        other1_names = [subgraph.nodes[n]['name'] for n in connector.intersection(other1)]
                        other2_names = [subgraph.nodes[n]['name'] for n in connector.intersection(other2)]
                        positivities = [loop_signs[(connector_c + 1) % 3], loop_signs[connector_c], loop_signs[(connector_c + 2) % 3]]
                        positive12 = not ishalfinhibition(subgraph, chosen_cycles[connector_c][0], other1, other2)
                        positive21 = not ishalfinhibition(subgraph, chosen_cycles[connector_c][0], other2, other1)
                        return logo_2bridged(other1_names, other2_names, positivities, positive12, positive21)
                    createimage(colored, args.images.format(current_id, kind, len(colored.edges)), logo_func)
                break
    