# HiLoop

Toolbox for exploring high-feedback motifs in gene regulatory networks. 

If you use HiLoop in your work, please cite:

> Benjamin Nordick and Tian Hong (2021). Identification, visualization, statistical analysis and mathematical modeling of high-feedback loops in gene regulatory networks.
> *BMC Bioinformatics* **22**, 481. https://doi.org/10.1186/s12859-021-04405-z

## Setup

HiLoop is provided as a collection of Python scripts.
It is platform-agnostic, but its dependencies are generally easier to set up on Linux.
If you use Windows 10, you can install the convenient [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

Virtual environments are strongly recommended for managing PyPI package dependencies.
For best performance, set up two virtual environments, one using standard Python (CPython) and the other using [PyPy](https://www.pypy.org/).
CPython is required for the modeling/simulation library Tellurium, but PyPy is significantly faster for pure Python code.

For example, to set up the PyPy environment used for most HiLoop functionality:

    $ /PYPY3_EXTRACTION_PATH/bin/pypy3 -m venv ~/venv/hiloop_pypy
    $ source ~/venv/hiloop_pypy/bin/activate
    (hiloop_pypy) $ pip install -r requirements.txt 

If you encounter difficulties installing SciPy for PyPy, you may be missing system packages such as C/C++ or Fortran compilers, or the development resources for OpenBLAS, LAPACK, or Python.
It may be helpful to upgrade or install the pip packages `pip`, `setuptools`, `wheel`, or `pybind11` before attempting to install all packages listed in the requirements file.

To set up the CPython environment used for modeling (`deactivate` any previously activated virtual environment first):

    $ python3 -m venv ~/venv/hiloop_cpy
    $ source ~/venv/hiloop_cpy/bin/activate
    (hiloop_cpy) $ pip install -r requirements_modeling.txt 

Drawing network diagrams requires Graphviz; install `graphviz` with your package manager (or from [its website](https://graphviz.org/download/)).

## Common Tasks

All common tasks can be accomplished by invoking a script from the command line.
Ensure that you have an appropriate virtual environment (indicated in each subsection) activated first.
You can pass `--help` to a script for more details than are covered in this overview.

### Obtaining a network

This repository contains several ready-to-use transcriptional regulatory networks:

* `trrust.gxml`: [TRRUST](https://www.grnpedia.org/trrust/) version 2 full human network
* `huang_emt.gxml`: Epithelial-mesenchymal transition control network examined by [Huang *et al.* 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005456) (includes miRNA)
* `tcell.gxml`: T cell development network curated by [Kueh & Rothenberg 2011](https://onlinelibrary.wiley.com/doi/abs/10.1002/wsbm.162)
* `ye_tcell.gxml`: Alternative interpretation of a smaller part of the T cell development network used by [Ye *et al.* 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006855)

Similar files suffixed with `_scc` are the main strongly connected components of the respective network.
None of the other, minor strongly connected components in any of these networks contain more than one cycle.

If you would like to examine only a *specific* part of a network, you can extract the subnetwork induced by a set of nodes (with either HiLoop virtual environment activated):

    (hiloop_pypy) $ python subnetwork.py trrust.gxml my_small_network.gxml GATA1 CEBPA MYC WT1

`omnipath.gxml` is a very large network containing transcriptional, miRNA post-transcriptional, and protein destabilization post-translational interactions from the
[OmniPath database](https://omnipathdb.org/) ([Türei *et al.* 2021](https://www.embopress.org/doi/full/10.15252/msb.20209923)).
You would likely want to take a subnetwork rather than using it directly.

**Note:** The common networks whose data is provided in HiLoop format are not necessarily under the same license (MIT) as HiLoop.
They may have restrictions on e.g. commercial use.

HiLoop's command-line interface accepts networks in GraphML format.
Each node should have a string `name` attribute for the gene's name; each edge must have a Boolean `repress` attribute indicating whether it represents a repression.
You can [manipulate graphs](https://networkx.org/documentation/stable/reference/classes/digraph.html) and [save them](https://networkx.org/documentation/stable/reference/readwrite/graphml.html) using the NetworkX library or any other tool that supports GraphML.

### Counting high-feedback topologies

`countmotifs.py` (PyPy virtual environment recommended) systematically counts the high-feedback topologies in a network.

For even moderately large networks, computing all cycles is infeasible.
Therefore, you usually must decide on a maximum cycle length to consider.
For example, 5 is the limit of feasibility for the TRRUST2 network.
You may also set an upper limit on the size of topology you would consider interesting.
Because nodes outside the main strongly connected component cannot be involved in any cycle that is fused with a cycle in the SCC, they cannot be involved in any high-feedback topologies.
Therefore, counting the SCC of any of the provided networks produces the same count of motif instances (except PFLs) as the full network, but is modestly more efficient.

To count high-feedback motifs of at most 10 nodes in the TRRUST2 network:

    (hiloop_pypy) $ python countmotifs.py trrust_scc.gxml --maxcycle 5 --maxnodes 10 --checknfl
    Finding cycles
    Creating cycle intersection graphs
    843 cycles, 154968 node sharings
    Searching for Type I and mixed-sign high-feedback motifs
    Checking fused pairs
    Searching for Type II and MISA motifs
    PFL     399         
    Type1   378454
    Type2   337621
    MISA    142931
    MISSA   5185
    uMISSA  224
    MixHF   2058336
    Excite  77426

HiLoop reports counts of these motifs:

|Code    |Name                                    |Definition                                                                                   |
|--------|----------------------------------------|---------------------------------------------------------------------------------------------|
|`PFL`   |Positive feedback loop                  |Cycle with an even number of represssions                                                    |
|`Type1` |Type-I motif                            |Three fused PFLs (i.e. all sharing at least one node)                                        |
|`Type2` |Type-II motif                           |Two non-overlapping PFLs each sharing at least one node with a third PFL                     |
|`MISA`  |Mutual inhibition self-activation       |Subtype of Type-II in which both "bridging" segments of the connecting PFL are net-repressive|
|`MISSA` |Mutual inhibition single-self-activation|Two fused PFLs, exactly one of which contains any repression edges                           |
|`uMISSA`|Mini-MISSA                              |Subtype of MISSA in which the pure-activation cycle is a self-loop                           |
|`MixHF` |Mixed-sign high-feedback                |Three fused cycles, of which at least one is positive and at least one is negative           |
|`Excite`|Excitable                               |A positive feedback loop fused to a negative feedback loop                                   |

If you are not interested in any motifs that involve negative feedback loops, you can omit `--checknfl` to significantly accelerate the process.
Progress logging can be suppressed with `--quiet`.

### Testing for enrichment

`sampledpvalue.py` (PyPy virtual environment recommended) rapidly computes the enrichment of feedback-rich topologies in a network by generating random permutations and estimating the quantity of high-feedback topologies in each with a cycle tuple sampling approach.

You must again consider the maximum cycle length and motif instance size that you would like to test for.
Use the full, original network as input so that the space of possible alternative networks can be more thoroughly explored, providing a more conservative enrichment report.
Specify the number of permutations to consider and how many samples to take from each.
To reduce noise in the estimate of the input network's metrics, use `--basetrials` to specify how many sampling runs of the input to average together.
If the sign of regulations performed by some species should remain constant&mdash;e.g. miRNAs can only repress&mdash;provide a regex to `--fixedsign` to match those species' names.

For example, to test the EMT network for enrichment against 100,000 permutations:

    (hiloop_pypy) $ python sampledpvalue.py huang_emt.gxml 100000 1000 --basetrials 100 --maxcycle 5 --maxmotifsize 10 --fixedsign miR
            FracLE
    PFLs    0.02008
    PFL/FL  0.0
    Type 1  0.09939
    Type 2  0.00062
    MISA    0.08177
    MixHF   1.0
    T1/Fus3 0.0
    T2/Brdg 0.0
    Excite  1.0
    MISSA   0.58516
    MISSA/F 0.92578
    uMISSA  0.16399 

`FracLE` gives the empirical p-value, specifically the proportion of permuted networks that match or exceed the input network in each metric.
HiLoop computes enrichment of the below metrics, including some ratios; you should determine in advance which are meaningful to you.

|Code     |Description                                                                             |
|---------|----------------------------------------------------------------------------------------|
|`PFLs`   |Number of positive feedback loops                                                       |
|`PFL/FL` |Ratio of positive feedback loops to cycles                                              |
|`Type 1` |Number of Type-I topologies                                                             |
|`Type 2` |Number of Type-II topologies                                                            |
|`MISA`   |Number of mutual inhibition self-activation Type-II topologies                          |
|`MixHF`  |Number of mixed-sign high-feedback topologies                                           |
|`T1/Fus3`|Ratio of Type-I topologies to fused cycle triplets                                      |
|`T2/Brdg`|Ratio of Type-II topologies to arrangements of 2 independent cycles "bridged" by a third|
|`Excite` |Number of excitable topologies                                                          |
|`MISSA`  |Number of mutual inhibition single-self-activation topologies                           |
|`MISSA/F`|Ratio of MISSA topologies to fused pairs of PFLs                                        |
|`uMISSA` |Number of mini-MISSA topologies                                                         |

You can specify an output CSV file with `--saveraw` to record each permutation's metrics for further analysis (e.g. histogram plotting).
To only consider strongly connected permutations, use the `--scc` switch and provide a strongly connected network as input.

### Extracting example topologies

`examplemotifs.py` (PyPy virtual environment recommended) randomly selects instances of high-feedback topologies from a larger network.
It can save the selected subnetworks to GraphML files for dynamics investigations and/or render network diagrams highlighting the interconnected cycles.

Once again, a maximum cycle length must be set for large networks; subnetwork size may also be limited.
The `--images` (to save network diagram images) and `--networks` (to save networks as GraphML) options take Python format strings specifying the file name in terms of the motif type and a numeric ID.
For example, `{1}_{0}.png` would save the first Type-II example as `type2_1.png`.
The script runs until the requested number of examples for each requested motif have been found.

For example, to find 5 small (5 nodes at most) examples of Type-I and Type-II motifs in the TRRUST2 network, saving both images and networks:

    (hiloop_pypy) $ mkdir -p examples
    (hiloop_pypy) $ python examplemotifs.py trrust.gxml --find1 5 --find2 5 --maxcycle 5 --maxnodes 5 --images examples/{1}_{0}.png --networks examples/{1}_{0}.gxml

The `--logo` switch adds a summary diagram of the cycle interconnection pattern.
The `--reduceedges` switch allows using the same set of nodes multiple times by dropping some edges not involved in the selected cycles, thereby making the edge sets distinct.
Press Ctrl+C to interrupt the script if it is progressing too slowly with the limitations you applied.

### Searching for multiattractor or oscillatory dynamics

`multistability.py` (CPython virtual environment required) tests random parameterizations of a model determined by the small input network.
It saves multistable and/or oscillatory systems, as requested, to a JSON report for further analysis.
`summarizemultistability.py` creates hexbin plots, scatter-line plots, or heatmaps from such reports.

First, run the parameter sampling and modeling script.
Specify an input GraphML file (of a small network like an extracted example, not a large original network) and an output JSON file.
`--attractors` (default 2) sets the minimum number of distinct attractors for systems you are interested in.
`--psets` specifies how many parameterizations to test.

For example, to test 20,000 parameterizations of the small Ye *et al.* T cell network and report systems producing at least 3 attractors:

    (hiloop_cpy) $ mkdir -p modeling
    (hiloop_cpy) $ python multistability.py ye_tcell.gxml modeling/ye_tcell.json --attractors 3 --psets 20000

This will take some time.

Console output can be suppressed with `--quiet`.
The number of initial concentrations per species can be set with `--concs`.
You may want to decrease this from the default of 5 if your network has more than four nodes/species.
If you are interested in oscillations, set the reporting timestep `--dt` to less than 1.0 and the simulation length `--time` to at least 200.
The default reporting timestep of 5.0 provides better performance and high accuracy for point attractors, and can sometimes detect oscillatory attractors, but the oscillatory orbits will be poorly resolved.
Oscillations are unlikely to be consistently detected with the default simulation length of 50 time units.
`--oscillators` allows reporting parameterizations that produce at least your specified number of oscillatory attractors, regardless of whether the total number of attractors meets the `--attractors` setting.

Once a dynamics report has been produced, visualizations can be created.
For all visualizations, specify the input JSON report and the output image file.
If your network is capable of generating oscillations, but you are not interested in them, you can pass the `--pointonly` switch to ignore all systems containing oscillations.

For hexbin and scatter-line plots, use the `scatterplot` command.
The `--reduction` option to configure the axes accepts `X/Y` format, e.g. `TCF1/PU1`, or `pca` to show the first two principal components.
Continuing the Ye *et al.* example by displaying a hexbin plot of the attractors:

    (hiloop_cpy) $ python summarizemultistability.py modeling/ye_tcell.json modeling/ye_tcell_hex.png scatterplot --reduction TCF1/PU1

To make a scatter-line plot, pass the `--line` switch.
You can specify system types to highlight with the `--focus` option, e.g. `--focus 4` to focus on 4-attractor systems or `--focus 4att4ms` to focus only on 4-attractor systems in which 4 species concentrations' are monotonically correlated to each other.

    (hiloop_cpy) $ python summarizemultistability.py modeling/ye_tcell.json modeling/ye_tcell_scatter.png scatterplot --reduction TCF1/PU1 --line --focus 4

To add contour lines or shading to a hexbin or scatter-line plot respectively, pass the `--contour` switch, optionally specifying the proportion of attractor density to keep outside the lowest contour (default 0.1).

For cluster-heatmaps, use the `heatmap` command.
To add a colorbar for the heatmap concentration values, pass `--colorbar`.
Arc (connector) columns, indicating which attractors are part of the same system, can be added with `--connect arc`.
If you would like to remove a connector column or reduce the number of connectors, use `--connect-downsample` to specify the limit or retention probability for each kind of system.
Each downsampling directive takes the form `systemtype:probability%` or `systemtype:count`.
For example, `3:50%` shows connectors for only 50% of all 3-attractor systems, `3att4ms:10%` shows connectors for only 10% of 3-attractor 4-monotonically-correlated-species systems, and `3att2ms:5` shows connectors for exactly five 3-attractor 2-monotonically-correlated-species systems.
This heatmap displays all connectors for 4-attractor systems, but no connectors for other kinds:

    (hiloop_cpy) $ python summarizemultistability.py modeling/ye_tcell.json modeling/ye_tcell_heatmap.png heatmap --colorbar --connect arc --connect-downsample 4:100% else:0

To cluster species as well as attractors, use the `--bicluster` switch.
To merge similar attractors from different systems, use the `--fold` option to specify the minimum distance between displayed attractors.

Both `scatterplot` and `heatmap` support a `--downsample` option to limit the systems that are shown at all.
It accepts the same syntax as described for connector downsampling.
To color different kinds of systems differently, use the `--color` switch for the `scatterplot` command and the `--color-coordinate` switch for the `heatmap` command.
If both commands used the same downsampling configuration, the heatmap's connector column labels will be color-coordinated to the lines in the scatter-line plot.
