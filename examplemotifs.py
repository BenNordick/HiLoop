import argparse
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
parser.add_argument('--images', type=str, help='output file pattern {id, edges, type}')
parser.add_argument('--printnodes', action='store_true', help='print distinct node names')
parser.add_argument('--maxedges', type=int, help='maximum number of unique edges in examples')
parser.add_argument('--maxnodes', type=int, help='maximum size of network to attempt cycle detection for')
parser.add_argument('--maxsharing', type=int, help='maximum number of nodes in common with an already drawn subnetwork')
args = parser.parse_args()

graph = nx.convert_node_labels_to_integers(nx.read_graphml(args.file))
type1, type2 = 0, 0
seen = []
printed_nodes = set()

def printnewnodes(nodes):
    if not args.printnodes:
        return
    for n in nodes.difference(printed_nodes):
        printed_nodes.add(n)
        print(graph.nodes[n]['name'])

while type1 < args.count1 or type2 < args.count2:
    feasible = graph if args.maxnodes is None else permutenetwork.randomsubgraph(graph, args.maxnodes)
    cycles = [(cycle, ispositive(feasible, cycle)) for cycle in nx.algorithms.simple_cycles(feasible)]
    if len(cycles) < 3:
        continue
    chosen = random.sample(cycles, 3)
    if [entry[1] for entry in chosen] != [True] * 3:
        continue
    chosen_cycles = [entry[0] for entry in chosen]
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
    subgraph = graph.subgraph(used_nodes)
    if args.maxedges is not None and len(subgraph.edges) > args.maxedges:
        continue
    if type1 < args.count1:
        intersection = cycle_sets[0].intersection(cycle_sets[1]).intersection(cycle_sets[2])
        if len(intersection) > 0:
            seen.append(used_nodes)
            type1 += 1
            printnewnodes(used_nodes)
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
                if args.images:
                    colored = colorsubgraph(subgraph, chosen_cycles[(connector_c + 1) % 3], chosen_cycles[connector_c], chosen_cycles[(connector_c + 2) % 3])
                    rendergraph.rendergraph(colored, args.images.format(type2, len(colored.edges), 2), in_place=True)
                break
    