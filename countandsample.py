from countmotifs import countmotifs
from liuwangcycles import cyclesgenerator
import networkx as nx
import os
import sampledpvalue

header = 'nodes,edges,repressions,samples,maxcyclelen,maxmotifsize,cycles,pfls,count1,count2,sample1,sample2\n'

def countandsample(graph, output, samples=10000, max_cycle_length=None, max_motif_size=None):
    cycles = 0
    for _ in cyclesgenerator(graph, max_cycle_length):
        cycles += 1
    pfls, actual_type1, actual_type2 = countmotifs(graph, max_cycle_length, max_motif_size)
    _, type1_samples, type2_samples = sampledpvalue.summarize(graph, samples, max_motif_size, max_cycle_length)
    repressions = len([0 for edge in graph.edges if graph.edges[edge]['repress']])
    result = f'{len(graph.nodes)},{len(graph.edges)},{repressions},{samples},{max_cycle_length},{max_motif_size}'
    result += f',{cycles},{pfls},{actual_type1},{actual_type2},{type1_samples},{type2_samples}'
    if output is None:
        print(result)
    elif isinstance(output, str):
        existing_file = os.path.isfile(output)
        with open(output, 'a') as f:
            if not existing_file:
                f.write(header)
            f.write(result)
            f.write('\n')
    else:
        output.write(result)
        output.write('\n')
