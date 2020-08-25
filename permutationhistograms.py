import json
import networkx as nx
import permutenetwork
import summarizenetwork
import sys

def permutationhistograms(base, runs):
    counts = {'pfls': dict(), 'type1': dict(), 'type2': dict(), 'sum12': dict()}
    for _ in range(runs):
        permutation = permutenetwork.permutenetwork(base)
        result = summarizenetwork.summarizenetwork(permutation)
        for k, v in result.items():
            counts[k][v] = counts[k][v] + 1 if v in counts[k] else 1
    histograms = dict()
    for metric, occurrences in counts.items():
        histogram = [0] * (max(occurrences.keys()) + 1)
        for value, count in occurrences.items():
            histogram[value] = count
        histograms[metric] = histogram
    return histograms

if __name__ == "__main__":
    base = nx.read_graphml(sys.argv[1])
    runs = int(sys.argv[2])
    histograms = permutationhistograms(base, runs)
    print(json.dumps(histograms))
