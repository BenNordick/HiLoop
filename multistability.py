import argparse
import collections
from itertools import product
import json
import warnings
import networkx as nx
import numpy as np
import re
import scipy.signal as ssig
import tellurium as te

# CPython virtual environment required

supported_models = ['multiplicative_hill', 'additive_hill', 'multiplicative_activation']

def coalesce_adjacent(points):
    """
    Remove all but the last number in each run of consecutive integers in an ascending list of integers.

    This is used by describe_attractor to detect when a time course has the same value at different times (necessary for oscillations).
    """
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
    """
    Categorize and characterize the endpoint of a simulated time course.

    Arguments:
    - result: Tellurium output (2D NumPy array where the first column is the time)
    - dt: reporting timestep
    - print_results: whether to display warnings when the simulation was too short or coarse to characterize the attractor

    For point attractors, returns a 1D NumPy array of concentrations.
    For oscillators, returns a dictionary:
    - "species": list of dicts (one per species), "min" is the minimum concentration, "max" is the maximum, "ftpeaks" is a dict of Fourier peak locations to magnitudes
    - "orbit": 2D NumPy array of concentrations covering one orbit of the oscillation
    If the attractor cannot be characterized reliably, returns None.
    """
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
                    warnings.warn('late-starting oscillation')
                return None
            last_quarter = result[-round(result.shape[0] / 4):, 1:]
            lq_range = np.max(last_quarter, axis=0) - np.min(last_quarter, axis=0)
            midpoint = round(result.shape[0] / 2)
            last_half = result[midpoint:, 1:]
            lh_range = np.max(last_half, axis=0) - np.min(last_half, axis=0)
            if np.any((lh_range - lq_range > 0.02) | (lh_range / lq_range > 1.5)):
                if print_results:
                    warnings.warn('interrupted dampening')
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
                warnings.warn('unstable endpoint without oscillation')
            return None

def equivalent_attractors(a, b):
    """Determine whether two attractors (in describe_attractor format) are equivalent."""
    if isinstance(a, dict) != isinstance(b, dict):
        # One attractor is an oscillator, but the other is a single point
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
    """Turn an attractor (in describe_attractor format) into an object tree that can be serialized to JSON."""
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

def findattractors(runner, ic_sets, time, dt, print_warnings=False):
    """
    Find the distinct attractors that can be produced by the given model under the given initial conditions.

    Arguments:
    - runner: Tellurium model (parameters will not be changed)
    - ic_sets: iterable of initial condition vectors to try
    - time: simulation time length
    - dt: reporting timestep
    - print_warnings: whether to have describe_attractor print warnings when simulation endpoints cannot be characterized

    Returns a list of distinct attractors.
    """
    sols = []
    points = round(time / dt) + 1
    for ii, ini_comb in enumerate(ic_sets):
        for iv, v in enumerate(runner.fs()):
            runner[v] = ini_comb[iv] # Set initial condition
        sim_points = runner.simulate(start=0, end=time, points=points)
        attractor = describe_attractor(sim_points, dt, print_warnings)
        if attractor is not None and all(not equivalent_attractors(attractor, e) for e in sols):
            sols.append(attractor)
    return sols

def findmultistability(runner, n_pts1d=5, n_psets=1000, min_attractors=2, min_oscillators=None, time=50, dt=5, fix_params=None, ignore_ptypes='g', print_results=False):
    """
    Test many random parameterizations of a model and report attractors produced by each.

    Arguments:
    - runner: Tellurium model (parameters will be changed)
    - n_pts1d: how many initial concentrations to try per gene/FS
    - n_psets: how many parameter sets to try
    - min_attractors: minimum number of attractors to report the system
    - min_oscillators: if set, minimum number of oscillatory attractors to report the system
    - time: simulation time length (increase if interested in oscillations)
    - dt: Tellurium reporting timestep (decrease if interested in oscillations, but does not affect accuracy for single-point attractors)
    - fix_params: dict of parameters to keep fixed: {parameter: value}
    - ignore_ptypes: character vector of parameter type IDs (K, k, r, n, g) to not change
    - print_results: whether to print desired multiattractor/oscillatory systems and warnings

    Returns a dictionary that can be serialized to JSON:
    - "species_names": list of species names in the order used by concentration vectors
    - "parameter_names": list of parameter names in the order used by parameter set reports
    - "psets": list of parameterizations that gave rise to the desired dynamics, each a dict:
      - "parameters": parameter values
      - "attractors": list of distinct attractors
    - "ftpoints": how many points were used for Fourier transforms in characterizing oscillations
    - "tested_psets": how many parameterizations were tested
    """

    # Define initial conditions for simulations
    # A uniform grid is used - if adapting the script for a large number of genes, consider replacing this with Latin Hypercube Sampling
    n_genes = runner.sv().size
    pts1d = np.linspace(0.0, 3.3, n_pts1d)
    ptsnd = np.tile(pts1d, n_genes).reshape(n_genes, -1)
    ini_combs = np.array(list(product(*ptsnd))) # All combinations of initial conditions
    points = round(time / dt) + 1

    # Ranges of parameter values to sample
    doms = {'K': [0.05, 4.5], 'k': [3.0, 3.3], 'r': [0.9, 0.99], 'n': [1, 6], 'g': [-1, 1], 'b': [0, 1]}
    rands = np.random.uniform(size=(n_psets, len(runner.ps()))) # Random numbers for picking/scaling parameters

    results = {'species_names': runner.fs(), 'parameter_names': runner.ps(), 'psets': [], 'ftpoints': points >> 1, 'tested_psets': n_psets}
    for i in range(n_psets):
        term_weight_groups = collections.defaultdict(list)
        for ip, p in enumerate(runner.ps()):
            if p[0] in ignore_ptypes:
                continue
            if p[0] == 'b':
                # Total term weights for each target in additive_hill must be normalized to 1
                group_key = p.split('_')[1] if '_' in p else p[1]
                term_weight_groups[group_key].append(p)
            dom = doms[p[0]]
            runner[p] = dom[0] + rands[i, ip] * (dom[1] - dom[0]) # Set parameter value
            if p[0] == 'n':
                runner[p] = np.round(runner[p], decimals=0) # n is usually an integer
            elif p[0] == 'g':
                runner[p] = 10.0 ** runner[p] # Uniform log distribution
        if fix_params:
            for k, v in fix_params.items():
                runner[k] = v
        for group in term_weight_groups.values():
            # Normalize term weights
            total_weight = sum(runner[p] for p in group)
            for p in group:
                if total_weight > 0:
                    runner[p] /= total_weight
                else:
                    runner[p] = 1 / len(group)
        sols = findattractors(runner, ini_combs, time, dt, print_results)
        n_oscillatory_sols = sum(1 for a in sols if isinstance(a, dict))
        if len(sols) >= min_attractors or (min_oscillators is not None and n_oscillatory_sols >= min_oscillators):
            if print_results:
                print([(s['species'] if isinstance(s, dict) else s) for s in sols], f'(pset {i})')
            result = {'parameters': [runner[p] for p in runner.ps()], 'attractors': [serialize_attractor(a) for a in sols]}
            results['psets'].append(result)
    return results

def networksb(network, model_form='multiplicative_hill'):
    """
    Convert a NetworkX directed graph to an Antimony string.

    By default, the multiplicative form of Hill functions is used.
    "additive_hill" adds (weighted) Hill terms instead of multiplying them.
    "multiplicative_activation" can be used if trying to reproduce a study that did not include Hill cross-terms or repressions.
    """
    if model_form not in supported_models:
        raise ValueError('Unsupported model_form')
    def safenodename(node):
        return network.nodes[node]['name'].replace('-', '').replace('.', '')
    def interactionid(node, regulator):
        return f'{safenodename(node)}_{safenodename(regulator)}'
    def expterm(node, regulator):
        interaction_id = interactionid(node, regulator)
        return f'(X_{safenodename(regulator)} / K_{interaction_id})^n_{interaction_id}'
    parts = []
    for i, node in enumerate(network.nodes):
        parts.append(f'J{i}: -> X_{safenodename(node)}; g_{safenodename(node)} * (k_{safenodename(node)} * ((1 - r_{safenodename(node)}) + r_{safenodename(node)}')
        if '_hill' in model_form:
            hill_terms = []
            interaction_ids = []
            for regulator in network.predecessors(node):
                if network.edges[regulator, node]['repress']:
                    hill_terms.append(f'1 / (1 + {expterm(node, regulator)})')
                else:
                    hill_terms.append(f'{expterm(node, regulator)} / (1 + {expterm(node, regulator)})')
                interaction_ids.append(interactionid(node, regulator))
            if model_form == 'multiplicative_hill':
                for term in hill_terms:
                    parts.append(' * ')
                    parts.append(term)
            elif len(hill_terms) > 0:
                parts.append(' * (')
                for i, (term, interaction) in enumerate(zip(hill_terms, interaction_ids)):
                    if i > 0:
                        parts.append(' + ')
                    parts.append(f'b_{interaction} * ({term})')
                parts.append(')')
        else:
            for regulator in network.predecessors(node):
                if network.edges[regulator, node]['repress']:
                    raise ValueError('Repression is not supported in multiplicative_activation form')
            activator_product = ' * '.join([expterm(node, regulator) for regulator in network.predecessors(node)])
            if len(activator_product) == 0:
                raise ValueError(f'{node}: Must have at least one activator in multiplicative_activation form')
            parts.append(f' * ({activator_product}) / (1 + {activator_product})')
        parts.append(f') - X_{safenodename(node)})\n')
        parts.append(f'X_{safenodename(node)} = 0.1\nk_{safenodename(node)} = 3.0\nr_{safenodename(node)} = 0.99\ng_{safenodename(node)} = 1.0\n')
    for regulator, target in network.edges:
        parts.append(f'\nK_{interactionid(target, regulator)} = 1.0; n_{interactionid(target, regulator)} = 4;')
        if model_form == 'additive_hill':
            node_regulators = network.in_degree(target)
            if node_regulators > 0:
                parts.append(f' b_{interactionid(target, regulator)} = {1 / node_regulators};')
    return ''.join(parts)
    
def networkmodel(network, model_form='multiplicative_hill'):
    """Convert a NetworkX directed graph to a Tellurium model."""
    return te.loada(networksb(network, model_form))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='input GraphML network (.gxml) or Antimony script (.sb)')
    parser.add_argument('output', type=str, help='output JSON file path')
    parser.add_argument('--form', type=str, choices=supported_models, default='multiplicative_hill', help='equation form for modeling networks')
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
        r = networkmodel(nx.read_graphml(args.input), model_form=args.form)
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
