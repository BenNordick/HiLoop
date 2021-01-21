import argparse
from itertools import product
import json
import networkx as nx
import numpy as np
import re
import scipy.signal as ssig
import tellurium as te

def coalesce_adjacent(points):
    if len(points) == 0:
        return points
    result = []
    last = points[0]
    for p in points[1:]:
        if p > last + 1:
            result.append(last)
        last = p
    result.append(last)
    return result

def describe_attractor(result, dt, print_results):
    recent = result[-round(10 / dt):, 1:]
    if np.linalg.norm(np.max(recent, axis=0) - np.min(recent, axis=0)) <= (4e-5) * recent.shape[1]:
        # Steady state
        return result[-1][1:]
    else:
        norms1 = np.linalg.norm(result[:, 1:] - result[-1][1:], axis=1)
        norms2 = np.linalg.norm(result[:, 1:] - result[-4][1:], axis=1)
        norm_cutoff = 0.03 * dt * recent.shape[1]
        matching_points = coalesce_adjacent(np.flatnonzero((norms1 < norm_cutoff) | (norms2 < norm_cutoff)))
        if len(matching_points) > 5:
            # Oscillation
            if matching_points[0] > result.shape[0] / 2:
                if print_results:
                    print('Warning: late-starting oscillation')
                return None
            last_quarter = result[-round(result.shape[0] / 4):, 1:]
            lq_range = np.max(last_quarter, axis=0) - np.min(last_quarter, axis=0)
            midpoint = round(result.shape[0] / 2)
            last_half = result[midpoint:, 1:]
            lh_range = np.max(last_half, axis=0) - np.min(last_half, axis=0)
            if np.any(lh_range - lq_range > 0.02):
                if print_results:
                    print('Warning: interrupted dampening')
                return None
            species_info = []
            actual_oscillation = False
            for species in result.colnames[1:]:
                series = result[species][midpoint:]
                series_min = np.min(series)
                series_max = np.max(series)
                if series_max - series_min < 0.01:
                    species_info.append({'min': series_min, 'max': series_max, 'ftpeaks': {}})
                    continue
                ft = np.abs(np.fft.rfft(series))
                peaks, props = ssig.find_peaks(ft, prominence=0.25, wlen=max(round(5 / dt), 5))
                peak_dict = dict(zip(peaks, props['prominences']))
                species_info.append({'min': series_min, 'max': series_max, 'ftpeaks': peak_dict})
                if len(peaks) > 0:
                    actual_oscillation = True
            if actual_oscillation:
                norms_to_next = np.linalg.norm(result[1:, 1:] - result[:-1, 1:], axis=1)
                orbit_end = result.shape[0] - 2
                left_start = False
                while norms_to_next[orbit_end - 1] < norms1[orbit_end] or not left_start:
                    if norms1[orbit_end] > 0.05:
                        left_start = True
                    orbit_end = orbit_end - 1
                return {'species': species_info, 'orbit': result[orbit_end:, 1:]}
            else:
                return None
        else:
            if print_results:
                print('Warning: unstable endpoint without oscillation')
            return None

def equivalent_attractors(a, b):
    if isinstance(a, dict) != isinstance(b, dict):
        return False
    if isinstance(a, dict):
        # Oscillatory attractor
        for i in range(len(a['species'])):
            if np.linalg.norm([a['species'][i]['min'] - b['species'][i]['min'], a['species'][i]['max'] - b['species'][i]['max']]) > 1e-1:
                return False
            unmatched_peaks = {k for k, v in b['species'][i]['ftpeaks'].items() if v > 1.5}
            for peak_pos, peak_prom in a['species'][i]['ftpeaks'].items():
                match = False
                for offset in range(-1, 2):
                    cand_pos = peak_pos + offset
                    if cand_pos in b['species'][i]['ftpeaks']:
                        if cand_pos in unmatched_peaks:
                            unmatched_peaks.remove(cand_pos)
                        cand_prom = b['species'][i]['ftpeaks'][cand_pos]
                        match = abs(cand_prom - peak_prom) < 3 or abs((cand_prom - peak_prom) / max(cand_prom, peak_prom)) < 0.1
                        break
                if (not match) and peak_prom > 1.5:
                    return False
            if len(unmatched_peaks) > 0:
                return False
        return True
    else:
        # Single-point attractor
        return np.linalg.norm(a - b) < 5e-1

def serialize_attractor(info):
    if isinstance(info, dict):
        for species in info['species']:
            ftpeaks = list(species['ftpeaks'].items())
            species['peaks'] = [int(fp[0]) for fp in ftpeaks]
            species['prominences'] = [float(fp[1]) for fp in ftpeaks]
            del species['ftpeaks']
        info['orbit'] = [[float(x) for x in r] for r in info['orbit']]
        return info
    else:
        return list(info)

# Adapted from a script written by Tian Hong 9/21/2020
def findmultistability(runner, n_pts1d=5, n_psets=1000, min_attractors=2, min_oscillators=None, time=50, dt=5, fix_params=None, ignore_ptypes='g', print_results=False):

    # Define initial conditions for simulations
    # A uniform grid is used. When number of genes is large, consider using Latin Hypercube Sampling
    n_genes = runner.sv().size
    pts1d = np.linspace(0.0, 3.3, n_pts1d)
    ptsnd = np.tile(pts1d, n_genes).reshape(n_genes, -1)
    ini_combs = np.array(list(product(*ptsnd))) # All combinations of initial conditions
    points = round(time / dt) + 1

    # Ranges of parameter values to sample
    doms = {'K': [0.05, 4.5], 'k': [3.0, 3.3], 'r': [0.9, 0.99], 'n': [1, 6], 'g': [-1, 1]}
    rands = np.random.uniform(size=(n_psets, len(runner.ps()))) # Random numbers for perturbing parameters

    results = {'species_names': runner.fs(), 'parameter_names': runner.ps(), 'psets': [], 'ftpoints': points >> 1, 'tested_psets': n_psets}
    for i in range(n_psets):
        for ip, p in enumerate(runner.ps()):
            if p[0] in ignore_ptypes:
                continue
            dom = doms[p[0]]
            runner[p] = dom[0] + rands[i, ip] * (dom[1] - dom[0]) # Set parameter value
            if p[0] == 'n':
                runner[p] = np.round(runner[p], decimals=0) # n is usually an integer
            elif p[0] == 'g':
                runner[p] = 10.0 ** runner[p] # Uniform log distribution
        if fix_params:
            for k, v in fix_params.items():
                runner[k] = v
        sols = []
        for ii, ini_comb in enumerate(ini_combs):
            for iv, v in enumerate(runner.fs()):
                runner[v] = ini_comb[iv] # Set initial condition
            sim_points = runner.simulate(start=0, end=time, points=points)
            attractor = describe_attractor(sim_points, dt, print_results)
            if attractor is not None and all(not equivalent_attractors(attractor, e) for e in sols):
                sols.append(attractor)
        n_oscillatory_sols = sum(1 for a in sols if isinstance(a, dict))
        if len(sols) >= min_attractors or (min_oscillators is not None and n_oscillatory_sols >= min_oscillators):
            if print_results:
                print([(s['species'] if isinstance(s, dict) else s) for s in sols], f'(pset {i})')
            result = {'parameters': [runner[p] for p in runner.ps()], 'attractors': [serialize_attractor(a) for a in sols]}
            results['psets'].append(result)
    return results

def networksb(network):
    '''Converts a NetworkX directed graph to an Antimony string.'''
    def safenodename(node):
        return network.nodes[node]['name'].replace('-', '').replace('.', '')
    parts = []
    for i, node in enumerate(network.nodes):
        parts.append(f'J{i}: -> X_{safenodename(node)}; g_{safenodename(node)} * (k_{safenodename(node)} * ((1 - r_{safenodename(node)}) + r_{safenodename(node)}')
        for regulator in network.predecessors(node):
            interaction_id = f'{safenodename(node)}_{safenodename(regulator)}'
            exp_term = f'(X_{safenodename(regulator)} / K_{interaction_id})^n_{interaction_id}'
            if network.edges[regulator, node]['repress']:
                parts.append(f' * 1 / (1 + {exp_term})')
            else:
                parts.append(f' * {exp_term} / (1 + {exp_term})')
        parts.append(f') - X_{safenodename(node)})\n')
        parts.append(f'X_{safenodename(node)} = 0.1\nk_{safenodename(node)} = 3.0\nr_{safenodename(node)} = 0.99\ng_{safenodename(node)} = 1.0\n')
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
    parser.add_argument('--oscillators', type=int, help='minimum number of oscillatory attractors to report')
    parser.add_argument('--time', type=int, default=50, help='length of simulation')
    parser.add_argument('--dt', type=float, default=5.0, help='time step length for runner.simulate')
    parser.add_argument('--concs', type=int, default=5, help='number of concentrations to test for each gene')
    parser.add_argument('--fix', type=str, help='file to get fixed parameters from')
    parser.add_argument('--fixfilter', type=str, help='regex to filter fixed parameter names by')
    parser.add_argument('--ignoretypes', type=str, default='g', help='parameter type letters to keep at defaults')
    parser.add_argument('--quiet', '-q', action='store_true', help='do not print results')
    args = parser.parse_args()
    if args.input.endswith('.sb'):
        with open(args.input) as f:
            r = te.loada(f.read())
    else:
        r = networkmodel(nx.read_graphml(args.input))
    fixed_params = None
    if args.fix:
        fixed_params = {}
        with open(args.fix) as f:
            for line in f:
                item, value = re.split(r"\s+", line, maxsplit=1)
                if (args.fixfilter is None) or re.match(args.fixfilter, item):
                    fixed_params[item] = float(value.rstrip())
    result = findmultistability(r, n_pts1d=args.concs, n_psets=args.psets, min_attractors=args.attractors, min_oscillators=args.oscillators, time=args.time, dt=args.dt, 
        fix_params=fixed_params, ignore_ptypes=args.ignoretypes, print_results=(not args.quiet))
    with open(args.output, 'w') as f:
        f.write(json.dumps(result))
