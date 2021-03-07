from countmotifs import countmotifs
from liuwangcycles import cyclesgenerator
import networkx as nx
import os
import sampledpvalue

header = 'nodes,edges,repressions,samples,maxcyclelen,maxmotifsize,cycles,pfls,count1,count2,countmisa,countumisa,sample1,sample2,samplemisa,sampleumisa\n'

def countandsample(graph, output, samples=10000, max_cycle_length=None, max_motif_size=None):
    cycles = 0
    for _ in cyclesgenerator(graph, max_cycle_length):
        cycles += 1
    pfls, actual_type1, actual_type2, actual_misa, actual_umisa, _, _ = countmotifs(graph, max_cycle_length, max_motif_size)
    _, _, type1_sampled, type2_sampled, misa_sampled, _, _, _, _, umisa_sampled, _ = sampledpvalue.summarize(graph, samples, max_motif_size, max_cycle_length)
    repressions = len([0 for edge in graph.edges if graph.edges[edge]['repress']])
    result = f'{len(graph.nodes)},{len(graph.edges)},{repressions},{samples},{max_cycle_length},{max_motif_size}'
    result += f',{cycles},{pfls},{actual_type1},{actual_type2},{actual_misa},{actual_umisa},{type1_sampled},{type2_sampled},{misa_sampled},{umisa_sampled}'
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
