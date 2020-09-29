import argparse
from itertools import product
import json
import networkx as nx
import numpy as np
import tellurium as te

# Adapted from a script written by Tian Hong 9/21/2020
def findmultistability(runner, n_pts1d=5, n_psets=1000, min_attractors=2, time=50, dt=5, print_results=False):

    # Define initial conditions for simulations
    # A uniform grid is used. When number of genes is large, consider using Latin Hypercube Sampling
    n_genes = runner.sv().size
    pts1d = np.linspace(0.0, 3.3, n_pts1d)
    ptsnd = np.tile(pts1d, n_genes).reshape(n_genes, -1)
    ini_combs = np.array(list(product(*ptsnd))) # all combinations of initial conditions

    # Ranges of parameter values to sample
    doms = {'K': (0.05, 4.5), 'k': [3.0, 3.3], 'r': [0.9, 0.99], 'n': [1, 6]}
    rands = np.random.uniform(size=(n_psets, len(runner.ps()))) # Random numbers for perturbing parameters

    results = {'species_names': runner.fs(), 'parameter_names': runner.ps(), 'psets': []}
    for i in range(n_psets):
        for ip, p in enumerate(runner.ps()):
            dom = doms[p[0]]
            runner[p] = dom[0] + rands[i, ip] * (dom[1] - dom[0]) # Set parameter value
            if p[0] == 'n':
                runner[p] = np.round(runner[p], decimals=0) # n is usually an integer
        sols = None
        for ii, ini_comb in enumerate(ini_combs):
            for iv, v in enumerate(runner.fs()):
                runner[v] = ini_comb[iv] # Set initial condition
            runner.simulate(start=0, end=time, points=round(time / dt))
            ss = runner.sv() # Get the steady state solution. Note: the variable names (unsorted) are r.fs()
            if sols is None:
                sols = ss
            else:
                if len(sols.shape) > 1:
                    dist = np.linalg.norm(sols - ss, axis=1)
                else:
                    dist = np.linalg.norm(sols - ss)
                if not np.any(dist < 5E-1): # Make sure we only save non-redundant steady states
                    sols = np.vstack((sols, ss))
        if len(sols.shape) > 1 and sols.shape[0] >= min_attractors: # Multiple attractors with a single parameter set
            result = {'parameters': [runner[p] for p in runner.ps()], 'attractors': [list(state) for state in sols]}
            results['psets'].append(result)
            if print_results:
                print(sols, f'(pset {i})')
    return results

def networksb(network):
    '''Converts a NetworkX directed graph to an Antimony string.'''
    def safenodename(node):
        return network.nodes[node]['name'].replace('-', '').replace('.', '')
    parts = []
    for i, node in enumerate(network.nodes):
        parts.append(f'J{i}: -> X_{safenodename(node)}; k_{safenodename(node)} * ((1 - r_{safenodename(node)}) + r_{safenodename(node)}')
        for regulator in network.predecessors(node):
            interaction_id = f'{safenodename(node)}_{safenodename(regulator)}'
            exp_term = f'(X_{safenodename(regulator)} / K_{interaction_id})^n_{interaction_id}'
            if network.edges[regulator, node]['repress']:
                parts.append(f' * 1 / (1 + {exp_term})')
            else:
                parts.append(f' * {exp_term} / (1 + {exp_term})')
        parts.append(f') - X_{safenodename(node)}\nX_{safenodename(node)} = 0.1\nk_{safenodename(node)} = 3.0\nr_{safenodename(node)} = 0.99\n')
    for regulator, target in network.edges:
        interaction_id = f'{safenodename(target)}_{safenodename(regulator)}'
        parts.append(f'\nK_{interaction_id} = 1.0; n_{interaction_id} = 4;')
    return ''.join(parts)
    
def networkmodel(network):
    '''Converts a NetworkX directed graph to a Tellurium model.'''
    return te.loada(networksb(network))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='input GraphML network (.gxml) or Antimony script (.sb)')
    parser.add_argument('output', type=str, help='output JSON file path')
    parser.add_argument('--psets', type=int, default=1000, help='number of parameter sets to try')
    parser.add_argument('--attractors', type=int, default=2, help='minimum number of attractors to report multistability')
    parser.add_argument('--time', type=int, default=50, help='length of simulation')
    parser.add_argument('--dt', type=float, default=5.0, help='time step length for runner.simulate')
    args = parser.parse_args()
    if args.input.endswith('.sb'):
        with open(args.input) as f:
            r = te.loada(f.read())
    else:
        r = networkmodel(nx.read_graphml(args.input))
    result = findmultistability(r, n_psets=args.psets, min_attractors=args.attractors, time=args.time, dt=args.dt, print_results=True)
    with open(args.output, 'w') as f:
        f.write(json.dumps(result))
