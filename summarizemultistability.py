import argparse
import collections
import copy
import json
import matplotlib.colors as mplcolors
import matplotlib.lines as mplline
import matplotlib.patches as mplpatch
import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
from sklearn import decomposition

def isoscillator(attractor):
    '''Determines whether the given attractor value is an oscillatory attractor.'''
    return isinstance(attractor, dict)

def caricatureattractor(attractor):
    '''Turn the given attractor information value (which might be an oscillation) into a single vector, for comparison.'''
    if isoscillator(attractor):
        return [(s['max'] + s['min']) / 2 for s in attractor['species']]
    else:
        return attractor

def caricatureattractors(attractors):
    '''Caricature each species in the given attractor set (nested list).'''
    return [caricatureattractor(a) for a in attractors]

def summarizeattractors(pset_report):
    '''Get a 2-tuple summarizing a set of attractors: attractor count, monotonic species count.'''
    attractors = caricatureattractors(pset_report['attractors'])
    species = len(attractors[0])
    correlated_species = set()
    most_monotonic_species = 0
    for i in range(species):
        if i in correlated_species:
            continue
        sorted_attractors = sorted(attractors, key=lambda a: a[i])
        correlated_species.add(i)
        monotonic_species = 1
        for j in set(range(species)).difference(correlated_species):
            attractor_concs = [a[j] for a in sorted_attractors]
            if attractor_concs == sorted(attractor_concs) or attractor_concs == sorted(attractor_concs, reverse=True):
                monotonic_species += 1
                correlated_species.add(j)
        most_monotonic_species = max(most_monotonic_species, monotonic_species)
    return len(attractors), most_monotonic_species

def categorizeattractors(report):
    '''Get a dictionary of attractor summary tuples to lists of their occurrences.'''
    summary_occurrences = collections.defaultdict(list)
    pset_list = report['psets'] if isinstance(report, dict) else report
    for pset in pset_list:
        summary_occurrences[summarizeattractors(pset)].append(pset)
    return summary_occurrences

def plotmultistability(report, label_counts=False, colorbar=True):
    '''Set up a multistability heatmap in the current pyplot.'''
    summary_occurrences = categorizeattractors(report)
    max_attractors = max(s[0] for s in summary_occurrences.keys())
    min_attractors = min(s[0] for s in summary_occurrences.keys())
    max_monotonic = len(report['species_names'])
    min_monotonic = 1
    width = max_attractors - min_attractors + 1
    x_range = range(min_attractors, max_attractors + 1)
    height = max_monotonic - min_monotonic + 1
    y_range = reversed(range(min_monotonic, max_monotonic + 1))
    heatmap_pixels = np.zeros((height, width), dtype=int)
    oscillators = np.zeros((height, width), dtype=int)
    for summary, occurrences in summary_occurrences.items():
        x = summary[0] - min_attractors
        y = max_monotonic - summary[1]
        heatmap_pixels[y][x] = len(occurrences)
        oscillators[y][x] = sum(1 for oc in occurrences if any(isoscillator(at) for at in oc['attractors']))
    fig, ax = plt.subplots()
    im = ax.imshow(heatmap_pixels, norm=mplcolors.LogNorm(vmax=heatmap_pixels.max()))
    if colorbar:
        fig.colorbar(im)
    ax.set_xticks(range(width))
    ax.set_yticks(range(height))
    ax.set_xticklabels([str(n) for n in x_range])
    ax.set_yticklabels([str(n) for n in y_range])
    ax.set_xlabel('Attractors')
    ax.set_ylabel('Monotonically correlated species')
    if label_counts:
        for y in range(height):
            for x in range(width):
                if heatmap_pixels[y][x] > 0:
                    text = str(heatmap_pixels[y][x])
                    if oscillators[y][x] > 0:
                        text = f'{text}\n({oscillators[y][x]} osc.)'
                    ax.text(x, y, text, ha='center', va='center', color='gray')

def plotattractors(report, reduction, connect_psets='none', filter_attractors=None, filter_correlated_species=None, downsample=None):
    reduction.prepare(report)
    summary_occurrences = categorizeattractors(report)
    filtered_psets = []
    random.seed(1)
    for summary, occurrences in summary_occurrences.items():
        attractors, monotonic = summary
        if filter_attractors is not None and attractors != filter_attractors:
            continue
        if filter_correlated_species is not None and monotonic != filter_correlated_species:
            continue
        if downsample is not None and attractors in downsample:
            filtered_psets.extend(o for o in occurrences if random.uniform(0, 1) < downsample[attractors])
        else:
            filtered_psets.extend(occurrences)
    xlabel, ylabel = reduction.labels()
    ax_main = None
    if connect_psets == 'line':
        fig, ax_main = plt.subplots()
        for pset in filtered_psets:
            pset_matrix = np.array(caricatureattractors(pset['attractors']))
            pset_xy = reduction.reduce(pset_matrix)
            sorted_attractors = pset_xy[pset_xy[:, 0].argsort()]
            line = ax_main.plot(sorted_attractors[:, 0], sorted_attractors[:, 1])
            pset_color = line[0].get_color()
            point_mask = [not isoscillator(a) for a in pset['attractors']]
            ax_main.scatter(pset_xy[point_mask, 0], pset_xy[point_mask, 1], color=pset_color)
            for osc in (a for a in pset['attractors'] if isoscillator(a)):
                vertices = np.array(osc['orbit'])
                projected_vertices = reduction.reduce(vertices)
                if projected_vertices.shape[0] >= 3:
                    projected_vertices = np.vstack((projected_vertices, projected_vertices[0, :]))
                polygon = mplpatch.Polygon(projected_vertices, color=pset_color, linewidth=1.5, linestyle='--', fill=False)
                ax_main.add_patch(polygon)
    else:
        points = reduction.reduce(psets_matrix(filtered_psets))
        cmap = copy.copy(plt.get_cmap('viridis'))
        cmap.set_under('white', 1.0)
        hex_args = {'linewidths': 0.2, 'norm': mplcolors.LogNorm(vmin=2), 'cmap': cmap}
        if connect_psets == 'arc':
            fig = plt.figure()
            grid = fig.add_gridspec(nrows=1, ncols=2, width_ratios=(6, 1), wspace=0.05)
            ax_main = fig.add_subplot(grid[0, 0])
            ax_main.hexbin(points[:, 0], points[:, 1], **hex_args)
            ax_arcs = fig.add_subplot(grid[0, 1], sharey=ax_main)
            ax_arcs.tick_params(labelbottom=False, labelleft=False, bottom=False)
            color_cycle = ax_arcs._get_lines.prop_cycler
            for pset in filtered_psets:
                pset_matrix = np.array(pset['attractors'])
                pset_xy = reduction.reduce(pset_matrix)
                sorted_ys = sorted(pset_xy[:, 1])
                height = random.uniform(0.2, 1.8)
                color = next(color_cycle)['color']
                for i in range(len(sorted_ys) - 1):
                    a, b = sorted_ys[i:(i + 2)]
                    ax_arcs.add_patch(mplpatch.Arc((0, (a + b) / 2), height, b - a, 180.0, 90.0, 270.0, edgecolor=color, linewidth=0.5))
        else:
            fig, ax_main = plt.subplots()
            ax_main.hexbin(points[:, 0], points[:, 1], **hex_args)
    ax_main.set_xlabel(xlabel)
    ax_main.set_ylabel(ylabel)

def psets_matrix(psets):
    full_matrix = None
    for pset in psets:
        numeric_attractors = np.array(caricatureattractors(pset['attractors']))
        if full_matrix is None:
            full_matrix = numeric_attractors
        else:
            full_matrix = np.vstack((full_matrix, numeric_attractors))
    return full_matrix

class PCA2D():
    def __init__(self):
        self.pca = decomposition.PCA(n_components=2)
    def prepare(self, report):
        self.pca.fit(psets_matrix(report['psets']))
    def reduce(self, matrix):
        return self.pca.transform(matrix)
    def labels(self):
        return 'PC1', 'PC2'

class AverageLog():
    def __init__(self, settings=None):
        self.settings = settings
    def prepare(self, report):
        self.names = report['species_names']
        if self.settings is None:
            raise NotImplementedError('You must specify genes for reduction axes')
        else:
            x, y = self.settings.split('/')
            self.x_components = [self._parsecomponent(report, c.strip()) for c in x.split(',')]
            self.y_components = [self._parsecomponent(report, c.strip()) for c in y.split(',')]
    def reduce(self, matrix):
        return np.stack((self._componentwisereduce(matrix, self.x_components), self._componentwisereduce(matrix, self.y_components)), 1)
    def labels(self):
        return ', '.join((self._componentname(c) for c in self.x_components)), ', '.join((self._componentname(c) for c in self.y_components))
    def _parsecomponent(self, report, text):
        if text.startswith('-'):
            text = text[1:]
            factor = -1
        else:
            factor = 1
        index = self.names.index(text) if text in self.names else self.names.index(f'X_{text}')
        return index, factor
    def _componentwisereduce(self, matrix, components):
        results = None
        for index, factor in components:
            component_log = np.log(matrix[:, index]) * factor
            if results is None:
                results = component_log
            else:
                results += component_log
        return np.exp(results / len(components))
    def _componentname(self, component):
        index, factor = component
        prefix = '-' if factor < 0 else ''
        name = self.names[index]
        return prefix + (name[2:] if name.startswith('X_') else name)

def plotheatmap(report, arcs=False, downsample=None, osc_orbits=1, fold_dist=None):
    gene_names = [(n[2:] if n.startswith('X_') else n) for n in report['species_names']]
    summary_occurrences = categorizeattractors(report)
    filtered_psets = []
    random.seed(1)
    for summary, occurrences in summary_occurrences.items():
        attractors, monotonic = summary
        if downsample is not None and attractors in downsample:
            filtered_psets.extend(o for o in occurrences if random.uniform(0, 1) < downsample[attractors])
        else:
            filtered_psets.extend(occurrences)
    if arcs:
        filtered_pset_types = categorizeattractors(filtered_psets)
        dendrogram_ratio = 3 / (13 + 2 * len(filtered_pset_types))
    else:
        dendrogram_ratio = 0.2
    matrix = None
    unique_fingerprints = None
    row_redundancies = None
    for pset in filtered_psets:
        pset['indexes'] = []
        for attractor in pset['attractors']:
            caricature = np.array(caricatureattractor(attractor))
            if fold_dist is not None:
                if isoscillator(attractor):
                    fingerprint = [100]
                    for species in attractor['species']:
                        fingerprint.append(species['min'])
                        fingerprint.append(species['max'])
                else:
                    fingerprint = [0]
                    for conc in attractor:
                        fingerprint.extend([conc, conc])
                fingerprint = np.array(fingerprint)
            if matrix is None:
                matrix = np.vstack((caricature, ))
                if fold_dist is not None:
                    unique_fingerprints = np.vstack((fingerprint, ))
                    row_redundancies = [1]
                pset['indexes'].append(0)
            elif fold_dist is None:
                pset['indexes'].append(matrix.shape[0])
                matrix = np.vstack((matrix, caricature))
            else:
                existing_indexes, = np.where(np.linalg.norm(unique_fingerprints - fingerprint, axis=1) < fold_dist * 2)
                if len(existing_indexes) > 0:
                    pset['indexes'].append(existing_indexes[0])
                    row_redundancies[existing_indexes[0]] += 1
                else:
                    pset['indexes'].append(matrix.shape[0])
                    matrix = np.vstack((matrix, caricature))
                    unique_fingerprints = np.vstack((unique_fingerprints, fingerprint))
                    row_redundancies.append(1)
    cg = sns.clustermap(matrix, col_cluster=False, cbar_pos=None, dendrogram_ratio=(dendrogram_ratio, 0), xticklabels=gene_names, yticklabels=False, cmap='seismic')
    matrix_display_ind = {v: k for k, v in enumerate(cg.dendrogram_row.reordered_ind)}
    heatmap_index = 1 if fold_dist is None else 2
    width_ratios = [2, 8]
    if arcs:
        width_ratios = [3, 10] + [2] * len(filtered_pset_types)
    if fold_dist is not None:
        width_ratios.insert(1, width_ratios[1] * len(gene_names) * 0.01)
        width_ratios[2] -= width_ratios[1]
    new_gs = plt.GridSpec(1, len(width_ratios), figure=cg.fig, width_ratios=width_ratios)
    cg.ax_heatmap.set_position(new_gs[heatmap_index].get_position(cg.fig))
    cg.ax_col_dendrogram.set_position(new_gs[0].get_position(cg.fig))
    if arcs:
        for fpt_id, summary in enumerate(sorted(filtered_pset_types.keys(), key=lambda am: am[0] * 100 + am[1], reverse=True)):
            ax_arcs = cg.fig.add_subplot(new_gs[heatmap_index + 1 + fpt_id], sharey=cg.ax_heatmap)
            ax_arcs.tick_params(labelbottom=False, labelleft=False, bottom=False)
            color_cycle = ax_arcs._get_lines.prop_cycler
            for pset_id, pset in enumerate(filtered_pset_types[summary]):
                color = next(color_cycle)['color']
                height = 1.75 - 0.2 * (pset_id % 8) + random.uniform(0, 0.1)
                rows = sorted(matrix_display_ind[i] for i in pset['indexes'])
                for i in range(len(rows) - 1):
                    a, b = rows[i:(i + 2)]
                    if a != b:
                        ax_arcs.add_patch(mplpatch.Arc((0, (a + b) / 2 + 0.5), height, b - a, 180.0, 90.0, 270.0, edgecolor=color, linewidth=0.7))
            ax_arcs.set_xlabel(f'{summary[0]} att.,\n{summary[1]} m.s.')
            for spine in ['top', 'right', 'bottom']:
                ax_arcs.spines[spine].set_visible(False)
    mesh = cg.ax_heatmap.collections[0]
    for pset in filtered_psets:
        for index, attr in zip(pset['indexes'], pset['attractors']):
            if isoscillator(attr):
                display_y = matrix_display_ind[index]
                orbit = np.array(attr['orbit'])
                for x in range(orbit.shape[1]):
                    x_stops = np.linspace(x, x + 1, orbit.shape[0] * osc_orbits)
                    color_stops = np.tile(np.vstack((orbit[:, x], orbit[:, x])), osc_orbits)
                    cg.ax_heatmap.pcolormesh(x_stops, [display_y, display_y + 1], color_stops, shading='gouraud', cmap=mesh.cmap, norm=mesh.norm)
    if fold_dist is not None:
        ax_redundancy = cg.fig.add_subplot(new_gs[1], sharey=cg.ax_heatmap)
        y_stops = np.arange(0, matrix.shape[0] + 1)
        reordered_redundancies = np.zeros((matrix.shape[0], 1))
        for i, redundancy in enumerate(row_redundancies):
            reordered_redundancies[matrix_display_ind[i], 0] = redundancy
        ax_redundancy.pcolormesh([0, 1], y_stops, reordered_redundancies, cmap='inferno')
        ax_redundancy.tick_params(labelbottom=False, labelleft=False, bottom=False)
        for spine in ['top', 'left', 'bottom']:
            ax_redundancy.spines[spine].set_visible(False)

def parse_downsample(arg):
    return {int(n): float(p) for n, p in [part.split(':') for part in arg.split(',')]} if arg else None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=str, help='input JSON report filename')
    parser.add_argument('graph', type=str, help='output graph image filename')
    subcmds = parser.add_subparsers(dest='command', required=True, help='kind of graph to make')
    table_parser = subcmds.add_parser('table')
    table_parser.add_argument('--counts', action='store_true', help='display counts in populated cells')
    table_parser.add_argument('--colorbar', action='store_true', help='show colorbar even when counts are displayed')
    scatterplot_parser = subcmds.add_parser('scatterplot')
    scatterplot_parser.add_argument('--attractors', type=int, help='filter parameter sets by number of attractors')
    scatterplot_parser.add_argument('--correlated', type=int, help='filter parameter sets by number of monotonically correlated species')
    scatterplot_parser.add_argument('--connect', choices=['none', 'line', 'arc'], default='none', help='how to connect attractors from the same parameter set')
    scatterplot_parser.add_argument('--reduction', type=str, help='species for dimensions: X1,X2/Y1,Y2 or "pca" to run PCA')
    scatterplot_parser.add_argument('--downsample', type=str, help='chance of keeping a parameter set with specified attractor count, e.g. 2:0.1,3:0.5')
    heatmap_parser = subcmds.add_parser('heatmap')
    heatmap_parser.add_argument('--arc', action='store_true', help='join multiattractor types with arcs')
    heatmap_parser.add_argument('--downsample', type=str, help='chance of keeping a parameter set with specified attractor count, e.g. 2:0.1,3:0.5')
    heatmap_parser.add_argument('--orbits', type=int, default=1, help='number of orbits to display for oscillatory attractors')
    heatmap_parser.add_argument('--fold', type=float, help='distance under which attractors will be combined into one heatmap row')
    args = parser.parse_args()
    with open(args.report) as f:
        report = json.loads(f.read())
    if args.command == 'table':
        plotmultistability(report, label_counts=args.counts, colorbar=(args.colorbar or not args.counts))
    elif args.command == 'scatterplot':
        reduction = PCA2D() if args.reduction == 'pca' else AverageLog(args.reduction)
        plotattractors(report, reduction, connect_psets=args.connect, filter_attractors=args.attractors, filter_correlated_species=args.correlated, downsample=parse_downsample(args.downsample))
    elif args.command == 'heatmap':
        plotheatmap(report, arcs=args.arc, downsample=parse_downsample(args.downsample), osc_orbits=args.orbits, fold_dist=args.fold)
    plt.savefig(args.graph, dpi=150)
    plt.close()
