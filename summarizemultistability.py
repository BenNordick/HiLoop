import argparse
import collections
import colorsys
import copy
import json
import matplotlib.collections as mplcollect
import matplotlib.colors as mplcolors
import matplotlib.lines as mplline
import matplotlib.patches as mplpatch
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
import mpl_toolkits.axes_grid1.inset_locator as mptinset
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

def specificrulevalue(ruleset, summary, default=None):
    specific = None
    attractors, monotonic = summary
    if summary in ruleset:
        specific = summary
    elif attractors in ruleset:
        specific = attractors
    return ruleset[specific] if specific in ruleset else default

def applydownsample(summary_occurrences, downsample):
    filtered_psets = []
    for summary, occurrences in summary_occurrences.items():
        n_psets = None
        if downsample is not None:
            limit_rule = specificrulevalue(downsample, summary, default=len(occurrences))
            if isinstance(limit_rule, int):
                n_psets = limit_rule
            else:
                percent = float(limit_rule.split('%')[0])
                n_psets = int(np.ceil(percent * len(occurrences) / 100))
        if n_psets is None or n_psets >= len(occurrences):
            filtered_psets.extend(occurrences)
        else:
            filtered_psets.extend(random.sample(occurrences, n_psets))
    return filtered_psets

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

def plotattractors(report, reduction, connect_psets=False, downsample=None, focus=None, focus_osc=False, square=False):
    reduction.prepare(report)
    random.seed(1)
    summary_occurrences = categorizeattractors(report)
    filtered_psets = applydownsample(summary_occurrences, downsample)
    xlabel, ylabel = reduction.labels()
    fig, ax_main = plt.subplots()
    if connect_psets:
        distinct_summaries = list(categorizeattractors(filtered_psets).keys())
        for i, pset in enumerate(filtered_psets):
            pset_matrix = np.array(caricatureattractors(pset['attractors']))
            pset_xy = reduction.reduce(pset_matrix)
            sorted_attractors = pset_xy[pset_matrix[:, 0].argsort(), :]
            point_mask = [not isoscillator(a) for a in pset['attractors']]
            has_oscillator = not all(point_mask)
            z = i
            pset_color = None
            linewidth = None
            oscwidth = 1.5
            dotsize = None
            if focus or focus_osc:
                summary = summarizeattractors(pset)
                hue, sat, lum, hue_vary_width = summaryhsl(distinct_summaries, summary)
                if (focus_osc and has_oscillator) or (focus and specificrulevalue(focus, summary, default=False)):
                    hue += random.uniform(0, hue_vary_width)
                    sat *= random.uniform(0.8, 1.0)
                    lum *= random.uniform(0.85, 1.1)
                    z += len(filtered_psets)
                else:
                    linewidth = 0.8
                    oscwidth = 1.1
                    dotsize = 10.0
                    lum = 1 - (1 - lum) * random.uniform(0.3, 0.5)
                    sat *= random.uniform(0.35, 0.45)
                    hue += random.uniform(0, hue_vary_width)
                pset_color = colorsys.hls_to_rgb(hue, lum, sat)
            line = ax_main.plot(sorted_attractors[:, 0], sorted_attractors[:, 1], lw=linewidth, color=pset_color, zorder=z)
            pset_color = line[0].get_color()
            ax_main.scatter(pset_xy[point_mask, 0], pset_xy[point_mask, 1], s=dotsize, color=pset_color, zorder=z)
            for osc in (a for a in pset['attractors'] if isoscillator(a)):
                vertices = np.array(osc['orbit'])
                projected_vertices = reduction.reduce(vertices)
                if projected_vertices.shape[0] >= 3:
                    projected_vertices = np.vstack((projected_vertices, projected_vertices[0, :]))
                polygon = mplpatch.Polygon(projected_vertices, color=pset_color, linewidth=oscwidth, linestyle='--', fill=False, zorder=z)
                ax_main.add_patch(polygon)
    else:
        points = reduction.reduce(psets_matrix(filtered_psets))
        cmap = copy.copy(plt.get_cmap('viridis'))
        cmap.set_under('white', 1.0)
        hex_args = {'linewidths': 0.2, 'norm': mplcolors.LogNorm(vmin=2), 'cmap': cmap}
        ax_main.hexbin(points[:, 0], points[:, 1], **hex_args)
    ax_main.axis('square' if square else 'equal')
    if reduction.zerobased('x'):
        ax_main.set_xlim(left=0)
    if reduction.zerobased('y'):
        ax_main.set_ylim(bottom=0)
    ax_main.xaxis.set_major_locator(mpltick.MultipleLocator(base=0.5))
    ax_main.yaxis.set_major_locator(mpltick.MultipleLocator(base=0.5))
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
    def zerobased(self, axis):
        return False

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
    def zerobased(self, axis):
        return len(self.x_components if axis == 'x' else self.y_components) == 1
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

def summaryhsl(all_summaries, summary):
    lowest_att = min(att for att, ms in all_summaries)
    highest_att = max(att for att, ms in all_summaries)
    att_range = highest_att - lowest_att + 1
    attractors, monotonic_species = summary
    lowest_ms = min(ms for att, ms in all_summaries if att == attractors)
    highest_ms = max(ms for att, ms in all_summaries if att == attractors)
    ms_range = highest_ms - lowest_ms + 1
    bin_width = 1 / (ms_range + 1) / att_range
    hue = ((highest_att - attractors) / att_range) + (highest_ms - monotonic_species) * bin_width
    variability_squeeze = (2 if att_range > 1 else 1) * (2 if ms_range > 1 else 1)
    return hue, 1, colorsys.ONE_THIRD, bin_width / variability_squeeze

def plotheatmap(report, arcs=None, downsample=None, arc_downsample=None, osc_orbits=1, fold_dist=None, bicluster=False):
    gene_names = [(n[2:] if n.startswith('X_') else n) for n in report['species_names']]
    random.seed(1)
    summary_occurrences = categorizeattractors(report)
    filtered_psets = applydownsample(summary_occurrences, downsample)
    if arcs:
        filtered_pset_types = categorizeattractors(filtered_psets)
        arc_pset_types = categorizeattractors(applydownsample(filtered_pset_types, arc_downsample)) if arc_downsample else filtered_pset_types
        dendrogram_ratio = 3 / (13 + 2 * len(arc_pset_types))
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
    plt.rc('font', size=18)
    gene_dendrogram_ratio = 0.1 if bicluster else 0
    cg = sns.clustermap(matrix, col_cluster=bicluster, cbar_pos=None, dendrogram_ratio=(dendrogram_ratio, gene_dendrogram_ratio), xticklabels=gene_names, yticklabels=False, cmap='seismic')
    matrix_display_ind = {v: k for k, v in enumerate(cg.dendrogram_row.reordered_ind)}
    gene_display_ind = {v: k for k, v in enumerate(cg.dendrogram_col.reordered_ind)} if bicluster else {n: n for n in range(len(gene_names))}
    heatmap_index = 1 if fold_dist is None else 2
    width_ratios = [2, 8]
    if arcs:
        width_ratios = [3, 10] + [2] * len(arc_pset_types)
    if fold_dist is not None:
        width_ratios.insert(1, width_ratios[1] * len(gene_names) * 0.01)
        width_ratios[2] -= width_ratios[1]
    rows = 2 if bicluster else 1
    main_row = rows - 1
    height_ratios = (1, 9) if bicluster else None
    new_gs = plt.GridSpec(rows, len(width_ratios), figure=cg.fig, width_ratios=width_ratios, height_ratios=height_ratios)
    cg.ax_heatmap.set_position(new_gs[main_row, heatmap_index].get_position(cg.fig))
    if bicluster:
        cg.ax_col_dendrogram.set_position(new_gs[0, heatmap_index].get_position(cg.fig))
    if arcs:
        for fpt_id, summary in enumerate(sorted(arc_pset_types.keys(), key=lambda am: am[0] * 100 + am[1], reverse=True)):
            ax_arcs = cg.fig.add_subplot(new_gs[main_row, heatmap_index + 1 + fpt_id], sharey=cg.ax_heatmap)
            ax_arcs.tick_params(labelbottom=False, labelleft=False, bottom=False)
            color_cycle = ax_arcs._get_lines.prop_cycler
            for pset_id, pset in enumerate(arc_pset_types[summary]):
                if arcs == 'straight':
                    height = 1.85 - 1.6 * pset_id / len(arc_pset_types[summary])
                    steepness = 0.18 * (1 - (height - 0.35) / 1.6)
                else:
                    height = 1.75 - 0.2 * (pset_id % 8) + random.uniform(0, 0.1)
                color = next(color_cycle)['color']
                rows = sorted(matrix_display_ind[i] for i in pset['indexes'])
                for i in range(len(rows) - 1):
                    a, b = rows[i:(i + 2)]
                    if a != b:
                        if arcs == 'straight':
                            segments = [[(0, a + 0.5), (height, a + 0.8 + steepness), (height, b + 0.2 - steepness), (0, b + 0.5)]]
                            lc = mplcollect.LineCollection(segments, colors=color, linewidths=0.8)
                            ax_arcs.add_collection(lc)
                        else:
                            ax_arcs.add_patch(mplpatch.Arc((0, (a + b) / 2 + 0.5), height, b - a, 180.0, 90.0, 270.0, edgecolor=color, linewidth=0.7))
            ax_arcs.set_xlabel(f'{summary[0]} att.,\n{summary[1]} m.s.')
            if arcs == 'straight':
                ax_arcs.set_xlim(0, 2)
            for spine in ['top', 'right', 'bottom']:
                ax_arcs.spines[spine].set_visible(False)
    mesh = cg.ax_heatmap.collections[0]
    max_orbit_len = 0
    for pset in filtered_psets:
        for attr in pset['attractors']:
            if isoscillator(attr):
                max_orbit_len = max(max_orbit_len, len(attr['orbit']))
    orbit_render_len = max_orbit_len * osc_orbits
    for pset in filtered_psets:
        for index, attr in zip(pset['indexes'], pset['attractors']):
            if isoscillator(attr):
                display_y = matrix_display_ind[index]
                orbit = np.array(attr['orbit'])
                for x in range(orbit.shape[1]):
                    display_x = gene_display_ind[x]
                    x_stops = np.linspace(display_x, display_x + 1, orbit_render_len)
                    color_stops = np.tile(np.vstack((orbit[:, x], orbit[:, x])), int(np.ceil(orbit_render_len / orbit.shape[0])))[:, :orbit_render_len]
                    cg.ax_heatmap.pcolormesh(x_stops, [display_y, display_y + 1], color_stops, shading='gouraud', cmap=mesh.cmap, norm=mesh.norm)
    if fold_dist is not None:
        ax_redundancy = cg.fig.add_subplot(new_gs[main_row, 1], sharey=cg.ax_heatmap)
        y_stops = np.arange(0, matrix.shape[0] + 1)
        reordered_redundancies = np.zeros((matrix.shape[0], 1))
        for i, redundancy in enumerate(row_redundancies):
            reordered_redundancies[matrix_display_ind[i], 0] = redundancy
        mesh = ax_redundancy.pcolormesh([0, 1], y_stops, reordered_redundancies, cmap='inferno')
        ax_redundancy.tick_params(labelbottom=False, labelleft=False, bottom=False)
        for spine in ['top', 'left', 'bottom']:
            ax_redundancy.spines[spine].set_visible(False)
        ax_cbar = mptinset.inset_axes(cg.ax_row_dendrogram, width='15%', height='20%', loc='lower left')
        cg.fig.colorbar(mesh, cax=ax_cbar)

def parse_systemtype(system_spec):
    if system_spec == 'else':
        return None
    elif 'att' in system_spec:
        att, ms_rest = system_spec.split('att')
        return (int(att), int(ms_rest.split('ms')[0]))
    else:
        return int(system_spec)

def parse_downsample(arg_list):
    def parse_one(arg):
        column, downsample = arg.split(':')
        if not downsample.endswith('%'):
            downsample = int(downsample)
        return (parse_systemtype(column), downsample)
    return dict(parse_one(arg) for arg in arg_list) if arg_list else None

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=str, help='input JSON report filename')
    parser.add_argument('graph', type=str, help='output graph image filename')
    parser.add_argument('--dpi', type=int, default=150, help='output bitmap image DPI')
    subcmds = parser.add_subparsers(dest='command', required=True, help='kind of graph to make')
    table_parser = subcmds.add_parser('table')
    table_parser.add_argument('--counts', action='store_true', help='display counts in populated cells')
    table_parser.add_argument('--colorbar', action='store_true', help='show colorbar even when counts are displayed')
    scatterplot_parser = subcmds.add_parser('scatterplot')
    scatterplot_parser.add_argument('--line', action='store_true', help='connect attractors from the same parameter set')
    scatterplot_parser.add_argument('--reduction', type=str, help='species for dimensions: X1,X2/Y1,Y2 (negatives allowed) or "pca" to run PCA')
    scatterplot_parser.add_argument('--downsample', nargs='+', type=str, help='chance of keeping a parameter set with specified type, e.g. 2:10% or 4att3ms:0')
    scatterplot_parser.add_argument('--focus', nargs='*', type=str, help='type(s) of parameter sets to focus on, e.g. 3att4ms or 4')
    scatterplot_parser.add_argument('--focus-osc', action='store_true', help='always focus parameter sets containing oscillations')
    scatterplot_parser.add_argument('--square', action='store_true', help='always use square axes')
    heatmap_parser = subcmds.add_parser('heatmap')
    heatmap_parser.add_argument('--connect', type=str, choices=['arc', 'straight'], help='connect attractors from the same parameter set')
    heatmap_parser.add_argument('--connect-downsample', '--cds', nargs='+', help='downsample connectors e.g. 3att4ms:10% or 4att2ms:5')
    heatmap_parser.add_argument('--downsample', nargs='+', type=str, help='chance of keeping a parameter set with specified type, e.g. e.g. 2:10% or 4att3ms:0')
    heatmap_parser.add_argument('--orbits', type=int, default=1, help='number of orbits to display for oscillatory attractors')
    heatmap_parser.add_argument('--fold', type=float, help='distance under which attractors will be combined into one heatmap row')
    heatmap_parser.add_argument('--bicluster', action='store_true', help='also cluster genes')
    args = parser.parse_args()
    with open(args.report) as f:
        report = json.loads(f.read())
    if args.command == 'table':
        plotmultistability(report, label_counts=args.counts, colorbar=(args.colorbar or not args.counts))
    elif args.command == 'scatterplot':
        reduction = PCA2D() if args.reduction == 'pca' else AverageLog(args.reduction)
        focus = {parse_systemtype(spec): True for spec in args.focus} if args.focus else None
        square = args.square or (args.reduction == 'pca')
        plotattractors(report, reduction, connect_psets=args.line, downsample=parse_downsample(args.downsample), focus=focus, focus_osc=args.focus_osc, square=square)
    elif args.command == 'heatmap':
        plotheatmap(report, arcs=args.connect, downsample=parse_downsample(args.downsample), arc_downsample=parse_downsample(args.connect_downsample),
                    osc_orbits=args.orbits, fold_dist=args.fold, bicluster=args.bicluster)
    plt.savefig(args.graph, dpi=args.dpi)
    plt.close()
