import argparse
import collections
import json
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import numpy as np

def summarizeattractors(pset_report):
    '''Get a 2-tuple summarizing a set of attractors: attractor count, monotonic species count.'''
    attractors = pset_report['attractors']
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

def plotmultistability(report, label_counts=False, colorbar=True):
    '''Sets up a multistability heatmap in the current pyplot.'''
    summary_occurrences = collections.defaultdict(int)
    for pset in report['psets']:
        summary_occurrences[summarizeattractors(pset)] += 1
    max_attractors = max(s[0] for s in summary_occurrences.keys())
    min_attractors = min(s[0] for s in summary_occurrences.keys())
    max_monotonic = len(report['species_names'])
    min_monotonic = 1
    width = max_attractors - min_attractors + 1
    x_range = range(min_attractors, max_attractors + 1)
    height = max_monotonic - min_monotonic + 1
    y_range = reversed(range(min_monotonic, max_monotonic + 1))
    heatmap_pixels = np.zeros((height, width), dtype=int)
    for summary, occurrences in summary_occurrences.items():
        heatmap_pixels[max_monotonic - summary[1]][summary[0] - min_attractors] = occurrences
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
                    ax.text(x, y, str(heatmap_pixels[y][x]), ha='center', va='center', color='gray')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('report', type=str, help='input JSON report filename')
    parser.add_argument('graph', type=str, help='output graph image filename')
    parser.add_argument('--counts', action='store_true', help='display counts in populated cells')
    parser.add_argument('--colorbar', action='store_true', help='show colorbar even when counts are displayed')
    args = parser.parse_args()
    with open(args.report) as f:
        report = json.loads(f.read())
    plotmultistability(report, label_counts=args.counts, colorbar=(args.colorbar or not args.counts))
    plt.savefig(args.graph)
    plt.close()
