from countmotifs import countmotifs
from liuwangcycles import cyclesgenerator
import networkx as nx
import os
import sampledpvalue

# PyPy virtual environment recommended for performance

header = 'nodes,edges,repressions,samples,maxcyclelen,maxmotifsize,cycles,pfls,count1,count2,countmisa,countmissa,countumissa,sample1,sample2,samplemisa,samplemissa,sampleumissa\n'

def countandsample(graph, output, samples=10000, max_cycle_length=None, max_motif_size=None):
    """
    Quantify high-feedback motifs in a network by both systematic counting and cycle tuple sampling, producing CSV output.
    Useful for calibrating and testing the estimation of counts from samples.

    Arguments:
    - graph: NetworkX DiGraph representing the network
    - output: how to report the CSV row: None to print, string to create or append to file with the specified name, or writable object to append to
    - samples: how many cycle tuples to sample for the sampling approach
    - max_cycle_length: maximum length of cycle to enumerate (needed for moderately large networks), or None to enumerate all cycles
    - max_motif_size: maximum number of nodes in a motif
    """
    cycles = 0
    for _ in cyclesgenerator(graph, max_cycle_length):
        cycles += 1
    pfls, actual_type1, actual_type2, actual_misa, actual_missa, actual_umissa, _, _ = countmotifs(graph, max_cycle_length, max_motif_size)
    _, _, type1_sampled, type2_sampled, misa_sampled, _, _, _, _, missa_sampled, _, umissa_sampled = sampledpvalue.summarize(graph, samples, max_motif_size, max_cycle_length)
    repressions = len([0 for edge in graph.edges if graph.edges[edge]['repress']])
    result = f'{len(graph.nodes)},{len(graph.edges)},{repressions},{samples},{max_cycle_length},{max_motif_size}'
    result += f',{cycles},{pfls},{actual_type1},{actual_type2},{actual_misa},{actual_missa},{actual_umissa}'
    result += f',{type1_sampled},{type2_sampled},{misa_sampled},{missa_sampled},{umissa_sampled}'
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
