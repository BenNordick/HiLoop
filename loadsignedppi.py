import networkx as nx
import pandas as pd
import sys

def loadsignedppi(df):
    graph = nx.DiGraph()
    for row in df.itertuples():
        graph.add_edge(row.Symbol1, row.Symbol2, repress=(row.Edge_Sign == '-'))
    for node in graph.nodes:
        graph.nodes[node]['name'] = node
    return nx.convert_node_labels_to_integers(graph)

if __name__ == "__main__":
    df = pd.read_excel(sys.argv[1])
    nx.write_graphml(loadsignedppi(df), sys.argv[2])
