import argparse
from itertools import product
import json
import numpy as np
import tellurium as te

# Adapted from a script written by Tian Hong 9/21/2020

def findmultistability(runner, n_pts1d=5, n_psets=1000, min_attractors=2, print_results=False):

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
            runner.simulate(start=0, end=50, points=10) # For a quick simulation. May need to confirm the accuracy when multiple attractors are obtained
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=str, help='output JSON file path')
    parser.add_argument('--psets', type=int, default=1000, help='number of parameter sets to try')
    parser.add_argument('--attractors', type=int, default=2, help='minimum number of attractors to report multistability')
    args = parser.parse_args()
    model_str = '''
        // Reactions:
        // X1: TP53;    X2: EZH2;   X3: KLF2;   X4: PPARG
        J00: -> X1; k1*((1-r1) + r1 * (X1/K11)^n11/(1+(X1/K11)^n11) * 1/(1+(X2/K12)^n12) * (X4/K14)^n14/(1+(X4/K14)^n14) ) - X1
        J01: -> X2; k2*((1-r2) + r2 * 1/(1+(X1/K21)^n21) ) - X2
        J02: -> X3; k3*((1-r3) + r3 * 1/(1+(X1/K31)^n31) ) - X3
        J03: -> X4; k4*((1-r4) + r4 * 1/(1+(X3/K43)^n43) ) - X4

        // Variable Initial Conditions:
        X1=0.1; X2=0.1; X3=0.1; X4=0.4

        // Parameter Values:
        k1=3.0; k2=3.0; k3=3.0; k4=3.0;
        r1=0.99; r2=0.99; r3=0.99; r4=0.99;
        K11=1.0; n11=4; K12=1.0; n12=4; K14=1.0; n14=4;
        K21=1.0; n21=4; K31=1.0; n31=4; K43=1.0; n43=4;
    ''' # Based on Type I example 14
    r = te.loada(model_str)
    result = findmultistability(r, n_psets=args.psets, min_attractors=args.attractors, print_results=True)
    with open(args.output, 'w') as f:
        f.write(json.dumps(result))
