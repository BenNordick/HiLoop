import argparse
import collections
from itertools import product
import json
import numpy as np
import pandas as pd
import tellurium as te
import time

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='input JSON file')
parser.add_argument('output', type=str, help='output CSV file')
parser.add_argument('--time', type=int, default=50, help='timespan to simulate')
parser.add_argument('--dt', type=float, default=0.01, help='time step for runner.simulate')
args = parser.parse_args()

with open(args.input) as f:
    psets = json.load(f)

model = te.loada('''
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
''')
pts1d = np.linspace(0.0, 3.3, 3)
ptsnd = np.tile(pts1d, 4).reshape(4, -1)
ini_combs = list(product(*ptsnd))
xs = [f'X{i + 1}' for i in range(4)]

results = []
start_time = time.process_time()
for pset in psets:
    for k, v in pset.items():
        model[k] = v
    for comb in ini_combs:
        result = collections.OrderedDict(pset)
        for i, x in enumerate(xs):
            model[x] = comb[i]
            result['start' + x] = comb[i]
        model.simulate(start=0, end=args.time, points=round(args.time / args.dt))
        for x in xs:
            result['end' + x] = model[x]
        results.append(result)
end_time = time.process_time()
print(end_time - start_time, 'seconds')

df = pd.DataFrame.from_records(results)
df.to_csv(args.output, index=False)
