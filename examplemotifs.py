import argparse
import liuwangcycles
from minimumtopologies import ispositive
import networkx as nx
import permutenetwork
import random
import rendergraph

def colorsubgraph(graph, r, y, b):
    return rendergraph.colorcycles(graph, [(r, 'red'), (y, 'gold'), (b, 'blue')]) 

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help='input GraphML file')
parser.add_argument('count1', type=int, help='how many Type I examples to find')
parser.add_argument('count2', type=int, help='how many Type II examples to find')
parser.add_argument('--images', type=str, help='output image file pattern {id, edges, type}')
parser.add_argument('--networks', type=str, help='output GraphML file pattern {id, type}')
parser.add_argument('--printnodes', action='store_true', help='print distinct node names')
parser.add_argument('--maxedges', type=int, help='maximum number of unique edges in examples')
parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
parser.add_argument('--maxcycle', type=int, help='maximum number of nodes in a cycle')
parser.add_argument('--maxsharing', type=int, help='maximum number of nodes in common with an already drawn subnetwork')
parser.add_argument('--negative', action='store_true', help='find fused negative feedback loops (not PFLs)')
parser.add_argument('--reduceedges', action='store_true', help='randomly drop some extra edges')
args = parser.parse_args()

graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
type1, type2 = 0, 0
seen = []
seen_edgesets = []
printed_nodes = set()

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

while type1 < args.count1 or type2 < args.count2:
    feasible = graph if args.maxnodes is None else permutenetwork.randomsubgraph(graph, args.maxnodes)
    cycles = [cycle for cycle in liuwangcycles.cyclesgenerator(feasible, args.maxcycle) if ispositive(feasible, cycle) ^ args.negative]
    if len(cycles) < 3:
        if not args.maxnodes:
            break
        continue
    chosen_cycles = random.sample(cycles, 3)
    cycle_sets = [frozenset(cycle) for cycle in chosen_cycles]
    used_nodes = cycle_sets[0].union(cycle_sets[1]).union(cycle_sets[2])
    if args.maxsharing is not None:
        duplicate = False
        for ns in seen:
            if len(ns.intersection(used_nodes)) > args.maxsharing:
                duplicate = True
                break
        if duplicate:
            continue
    subgraph = reduceedges(graph.subgraph(used_nodes), chosen_cycles)
    if args.maxedges is not None and len(subgraph.edges) > args.maxedges:
        continue
    edges_set = set(subgraph.edges)
    if edges_set in seen_edgesets:
        continue
    seen_edgesets.append(edges_set)
    if type1 < args.count1:
        intersection = cycle_sets[0].intersection(cycle_sets[1]).intersection(cycle_sets[2])
        if len(intersection) > 0:
            seen.append(used_nodes)
            type1 += 1
            printnewnodes(used_nodes)
            if args.networks:
                nx.write_graphml(subgraph, args.networks.format(type1, 1))
            if args.images:
                colored = colorsubgraph(subgraph, *chosen_cycles)
                for n in intersection:
                    colored.nodes[n]['penwidth'] = 2.0
                rendergraph.rendergraph(colored, args.images.format(type1, len(colored.edges), 1), in_place=True)
            continue
    if type2 < args.count2:
        for connector_c, connector in enumerate(cycle_sets):
            other1 = cycle_sets[(connector_c + 1) % 3]
            other2 = cycle_sets[(connector_c + 2) % 3]
            if connector.isdisjoint(other1) or connector.isdisjoint(other2):
                continue
            if other1.isdisjoint(other2):
                seen.append(used_nodes)
                type2 += 1
                printnewnodes(used_nodes)
                if args.networks:
                    nx.write_graphml(subgraph, args.networks.format(type2, 2))
                if args.images:
                    colored = colorsubgraph(subgraph, chosen_cycles[(connector_c + 1) % 3], chosen_cycles[connector_c], chosen_cycles[(connector_c + 2) % 3])
                    rendergraph.rendergraph(colored, args.images.format(type2, len(colored.edges), 2), in_place=True)
                break
    