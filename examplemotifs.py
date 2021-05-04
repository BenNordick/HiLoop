import argparse
from countmotifs import hasrepression
import io
import liuwangcycles
from minimumtopologies import ispositive, ismutualinhibition
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
    return rendergraph.colorcycles(graph, [(r, 'red', cyclestyle(graph, r)), (y, 'gold2', cyclestyle(graph, y)), (b, 'blue', cyclestyle(graph, b))])

labelnodeparams = {'fontsize': 9.0, 'width': 0.1, 'height': 0.1, 'margin': 0.05, 'fontname': 'DejaVuSerif'}

def logobase(**kwargs):
    ag = pygraphviz.AGraph(bgcolor='#D0D0D0', strict=False, directed=True, ranksep=0.3, **kwargs)
    ag.edge_attr['arrowsize'] = 0.8
    return ag

def logo_excitable(fusion_nodes):
    ag = logobase()
    ag.add_node('X', label=',\n'.join(fusion_nodes), **labelnodeparams)
    ag.add_edge('X', 'X', 0, color='red', style='dashed', arrowhead='tee', headport='ne', tailport='se')
    ag.add_edge('X', 'X', 1, color='blue', headport='sw', tailport='nw')
    return ag

def logo_missa(fusion_nodes, missa_targets):
    ag = logobase()
    ag.add_node('C1', label=',\n'.join(fusion_nodes), **labelnodeparams)
    ag.add_node('C2', label=',\n'.join(missa_targets), **labelnodeparams)
    ag.add_edge('C1', 'C1', color='blue', headport='nw', tailport='ne')
    ag.add_edge('C1', 'C2', color='red', arrowhead='tee')
    ag.add_edge('C2', 'C1', color='red', arrowhead='tee')
    return ag

def logo_3fused(fusion_nodes, positivities):
    styles = [('solid' if positive else 'dashed') for positive in positivities]
    tips = [('normal' if positive else 'tee') for positive in positivities]
    ag = logobase(nodesep=0.6)
    ag.add_node('L3', shape='point', width=0.001, color='blue')
    ag.add_node('X', label=',\n'.join(fusion_nodes), **labelnodeparams)
    ag.add_edge('L3', 'X', arrowhead=tips[2], color='blue', style=styles[2])
    ag.add_edge('X', 'L3', arrowhead='none', color='blue', style=styles[2])
    ag.add_node('L1', shape='point', width=0.001, color='red')
    ag.add_edge('X', 'L1', arrowhead='none', color='red', style=styles[0])
    ag.add_edge('L1', 'X', arrowhead=tips[0], color='red', style=styles[0])
    ag.add_node('L2', shape='point', width=0.001, color='gold')
    ag.add_edge('X', 'L2', arrowhead='none', color='gold', style=styles[1])
    ag.add_edge('L2', 'X', arrowhead=tips[1], color='gold', style=styles[1])
    return ag

def logo_2bridged(fusion1, fusion2, positivities, half12_positive, half21_positive):
    styles = [('solid' if positive else 'dashed') for positive in positivities]
    tips = [('normal' if positive else 'tee') for positive in positivities]
    ag = logobase()
    ag.add_node('C1', label=',\n'.join(fusion1), **labelnodeparams)
    ag.add_node('C2', label=',\n'.join(fusion2), **labelnodeparams)
    ag.add_edge('C1', 'C1', color='red', style=styles[0], arrowhead=tips[0])
    ag.add_edge('C1', 'C2', color='gold', style=styles[1], arrowhead=('normal' if half12_positive else 'tee'))
    ag.add_edge('C2', 'C1', color='gold', style=styles[1], arrowhead=('normal' if half21_positive else 'tee'))
    ag.add_edge('C2', 'C2', color='blue', style=styles[2], arrowhead=tips[2])
    return ag

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='input GraphML file')
parser.add_argument('-1', '--find1', type=int, default=0, help='how many Type I examples to find')
parser.add_argument('-2', '--find2', type=int, default=0, help='how many Type II examples to find')
parser.add_argument('-x', '--findmixed', type=int, default=0, help='how many mixed-sign high-feedback examples to find')
parser.add_argument('-m', '--findmisa', type=int, default=0, help='how many MISA Type II examples to find')
parser.add_argument('--findnegative1', type=int, default=0, help='how many negative Type I examples to find')
parser.add_argument('--findnegative2', type=int, default=0, help='how many negative Type II examples to find')
parser.add_argument('-e', '--findexcitable', type=int, default=0, help='how many excitable examples to find')
parser.add_argument('-s', '--findmissa', type=int, default=0, help='how many MISSA examples to find')
parser.add_argument('-u', '--findumissa', type=int, default=0, help='how many mini-MISSA examples to find')
parser.add_argument('--images', nargs='+', type=str, help='output image file pattern(s) {id, type, edges}')
parser.add_argument('--dpi', type=float, default=96.0, help='DPI of output images')
parser.add_argument('--logo', action='store_true', help='add motif summary/logo to saved images')
parser.add_argument('--networks', type=str, help='output GraphML file pattern {id, type}')
parser.add_argument('--printnodes', action='store_true', help='print distinct node names')
parser.add_argument('--maxedges', type=int, help='maximum number of unique edges in examples')
parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
parser.add_argument('--maxcycle', type=int, help='maximum number of nodes in a cycle')
parser.add_argument('--maxsharing', type=int, help='maximum number of nodes in common with an already selected subnetwork')
parser.add_argument('--reduceedges', action='store_true', help='randomly drop some extra edges')
parser.add_argument('--requirenodes', nargs='+', type=str, help='node(s) that must be present in the subnetwork')
parser.add_argument('--usesubgraph', nargs='+', type=str, help='nodes whose induced subgraph to search')
args = parser.parse_args()

if args.images and args.logo and any(name.endswith('.svg') for name in args.images):
    raise RuntimeError('Logos are not supported when rendering to SVG')

graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
type1, type2, mixed, misa, negative1, negative2, excitable, missa, umissa = 0, 0, 0, 0, 0, 0, 0, 0, 0
seen = []
seen_edgesets = []
printed_nodes = set()

if args.usesubgraph:
    node_ids = {graph.nodes[n]['name']: n for n in graph.nodes}
    used_ids = [node_ids[name] for name in args.usesubgraph]
    graph = nx.relabel_nodes(graph.subgraph(used_ids), {key: i for i, key in enumerate(used_ids)})

def shouldcheck2fused():
    return excitable < args.findexcitable or missa < args.findmissa or umissa < args.findumissa

def shouldcheck3fused():
    return type1 < args.find1 or mixed < args.findmixed or negative1 < args.findnegative1

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

def createimage(graph, filename_placeholders, logo_func):
    if args.logo:
        logo_ag = logo_func()
        logo_ag.graph_attr['dpi'] = args.dpi            
        logo_ag.layout('dot')
        logo_bytes = logo_ag.draw(format='png')
        logo = Image.open(io.BytesIO(logo_bytes))
        logo_w, logo_h = logo.size
        main_ag = rendergraph.graphvizify(graph, in_place=True, layout=None)
        main_ag.graph_attr['dpi'] = args.dpi
        main_ag.add_node('Logo', label='', shape='box', style='filled', color='#D0D0D0', width=(logo_w / args.dpi), height=(logo_h / args.dpi))
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
        for pattern in args.images:
            main.save(pattern.format(*filename_placeholders))
    else:
        graph.graph['dpi'] = args.dpi
        for pattern in args.images:
            rendergraph.rendergraph(graph, pattern.format(*filename_placeholders), in_place=True)

cycles = None
if args.maxnodes is None:
    cycles = [(cycle, ispositive(graph, cycle), hasrepression(graph, cycle)) for cycle in liuwangcycles.cyclesgenerator(graph, args.maxcycle)]
    if len(cycles) < 3:
        print('Insufficient cycles for high feedback')

def pickcycles(count):
    if count > len(cycles):
        return None
    chosen_cycles = random.sample(cycles, count)
    cycle_sets = [frozenset(cycle) for cycle, _, _ in chosen_cycles]
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
    subgraph = reduceedges(graph.subgraph(used_nodes), [cycle for cycle, _, _ in chosen_cycles])
    if args.maxedges is not None and len(subgraph.edges) > args.maxedges:
        return None
    edges_set = set(subgraph.edges)
    if edges_set in seen_edgesets:
        return None
    def consume():
        seen.append(used_nodes)
        printnewnodes(used_nodes)
        seen_edgesets.append(edges_set)
    return chosen_cycles, cycle_sets, used_nodes, subgraph, consume

while shouldcheck2fused() or shouldcheck3fused() or shouldcheckbridged():
    if args.maxnodes is not None:
        feasible = permutenetwork.randomsubgraph(graph, args.maxnodes)
        cycles = [(cycle, ispositive(feasible, cycle), hasrepression(feasible, cycle)) for cycle in liuwangcycles.cyclesgenerator(feasible, args.maxcycle)]
    if shouldcheck2fused():
        pick_results = pickcycles(2)
        if pick_results is not None:
            chosen_cycles, cycle_sets, used_nodes, subgraph, consume = pick_results
            intersection = cycle_sets[0].intersection(cycle_sets[1])
            if len(intersection) > 0:
                kind = None
                current_id = None
                if chosen_cycles[0][1] != chosen_cycles[1][1]:
                    if excitable < args.findexcitable:
                        kind = 'excitable'
                        excitable += 1
                        current_id = excitable
                        red, blue = (chosen_cycles[1][0], chosen_cycles[0][0]) if chosen_cycles[0][1] else (chosen_cycles[0][0], chosen_cycles[1][0])
                elif chosen_cycles[0][1] and chosen_cycles[1][1] and (chosen_cycles[0][2] != chosen_cycles[1][2]):
                    if missa < args.findmissa or umissa < args.findumissa:
                        kind = 'missa'
                        missa += 1
                        current_id = missa
                        red, blue = (chosen_cycles[1][0], chosen_cycles[0][0]) if chosen_cycles[1][2] else (chosen_cycles[0][0], chosen_cycles[1][0])
                        if umissa < args.findumissa:
                            if len(cycle_sets[0]) == 1 or len(cycle_sets[1]) == 1:
                                kind = 'minimissa'
                                umissa += 1
                                current_id = umissa
                            elif args.findmissa == 0:
                                kind = None
                if kind is not None:
                    consume()
                    if args.networks:
                        nx.write_graphml(subgraph, args.networks.format(current_id, kind))
                    if args.images:
                        colored = colorsubgraph(subgraph, red, [], blue)
                        for n in intersection:
                            colored.nodes[n]['penwidth'] = 2.0
                        if kind == 'excitable':
                            logo_func = lambda: logo_excitable([subgraph.nodes[n]['name'] for n in intersection])
                        else:
                            def logo_func():
                                intersection_i = red.index(next(iter(intersection))) # The mutual inhibition cycle is red
                                i = intersection_i
                                mi_len = len(red)
                                while True:
                                    next_i = (i + 1) % mi_len
                                    if subgraph.edges[red[i], red[next_i]]['repress']:
                                        mutual_inhibition_target_name = subgraph.nodes[red[next_i]]['name']
                                        break
                                    if next_i == intersection_i:
                                        raise AssertionError('Could not find mutual inhibition for MISSA')
                                    i = next_i
                                return logo_missa([subgraph.nodes[n]['name'] for n in intersection], [mutual_inhibition_target_name])
                        createimage(colored, (current_id, kind, len(colored.edges)), logo_func)
    if len(cycles) < 3 or not (shouldcheck3fused() or shouldcheckbridged()):
        continue
    pick_results = pickcycles(3)
    if pick_results is None:
        continue
    chosen_cycles, cycle_sets, used_nodes, subgraph, consume = pick_results
    loop_signs = [positive for cycle, positive, _ in chosen_cycles]
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
                if mixed < args.findmixed:
                    kind = 'mixed'
                    mixed += 1
                    current_id = mixed
            if kind is None:
                continue
            consume()
            if args.networks:
                nx.write_graphml(subgraph, args.networks.format(current_id, kind))
            if args.images:
                colored = colorsubgraph(subgraph, *[cycle for cycle, _, _ in chosen_cycles])
                for n in intersection:
                    colored.nodes[n]['penwidth'] = 2.0
                logo_func = lambda: logo_3fused({subgraph.nodes[n]['name'] for n in intersection}, loop_signs)
                createimage(colored, (current_id, kind, len(colored.edges)), logo_func)
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
                        is_misa = ismutualinhibition(subgraph, chosen_cycles[connector_c][0], other1, other2)
                        kind = 'type2'
                        type2 += 1
                        current_id = type2
                        if misa < args.findmisa:
                            if is_misa:
                                kind = 'type2misa'
                                misa += 1
                                current_id = misa
                            elif args.find2 == 0:
                                kind = None
                elif not any(loop_signs):
                    if negative2 < args.findnegative2:
                        kind = 'negative2'
                        negative2 += 1
                        current_id = negative2
                if kind is None:
                    continue
                consume()
                if args.networks:
                    nx.write_graphml(subgraph, args.networks.format(current_id, kind))
                if args.images:
                    colored = colorsubgraph(subgraph, chosen_cycles[(connector_c + 1) % 3][0], chosen_cycles[connector_c][0], chosen_cycles[(connector_c + 2) % 3][0])
                    def logo_func():
                        connector_cycle = chosen_cycles[connector_c][0]
                        cc_len = len(connector_cycle)
                        i = 0
                        while not (connector_cycle[i] in other1 and connector_cycle[(i + 1) % cc_len] not in other1):
                            i += 1
                        cycle1_rep = subgraph.nodes[connector_cycle[i % cc_len]]['name']
                        positive12 = True
                        hit_c2 = False
                        while True:
                            positive12 ^= subgraph.edges[connector_cycle[i % cc_len], connector_cycle[(i + 1) % cc_len]]['repress']
                            next_in_c2 = connector_cycle[(i + 1) % cc_len] in other2
                            i += 1
                            if hit_c2 and not next_in_c2:
                                if is_misa:
                                    raise AssertionError('Could not find mutual repression for MISA example')
                                i -= 1 # For Type II, go back into the other cycle and accept the quasi-MISA
                                break
                            if next_in_c2:
                                hit_c2 = True
                                if loop_signs[connector_c] and positive12 == is_misa:
                                    continue
                                break
                        cycle2_rep = subgraph.nodes[connector_cycle[i % cc_len]]['name']
                        positive21 = positive12 ^ (not loop_signs[connector_c])
                        positivities = [loop_signs[(connector_c + 1) % 3], loop_signs[connector_c], loop_signs[(connector_c + 2) % 3]]
                        return logo_2bridged([cycle1_rep], [cycle2_rep], positivities, positive12, positive21)
                    createimage(colored, (current_id, kind, len(colored.edges)), logo_func)
                break
    