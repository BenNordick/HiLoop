import argparse
import networkx as nx

# Simple command-line utility to take the subgraph of a large network induced by specified nodes

parser = argparse.ArgumentParser()
parser.add_argument('network', type=str, help='input GraphML file')
parser.add_argument('output', type=str, help='output subnetwork GraphML file')
parser.add_argument('nodes', nargs='+', type=str, help='nodes for induced subnetwork')
args = parser.parse_args()

network = nx.read_graphml(args.network)
node_ids = {network.nodes[n]['name']: n for n in network.nodes}
selected_node_ids = [node_ids[n] for n in args.nodes]
subnetwork = nx.convert_node_labels_to_integers(network.subgraph(selected_node_ids))
nx.write_graphml(subnetwork, args.output)
