
import numpy as np
import math
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from potts.structures import Lattice
from potts.models import HomologicalPercolation
from potts import Chain, Tape, _version
import sys

# Construct lattice object.
field = 2
L = Lattice().fromCorners([20,20], field=field)
print(np.sqrt(field)/(1 + np.sqrt(field)))

# Set up Model and Chain.
homology = 1
rank = math.comb(len(L.corners), homology)
SW = HomologicalPercolation(L, homology=homology)
N = 1000
M = Chain(SW, steps=N)

ratios = []

for (spins, essentials, satisfied) in M.progress():
	ratios.append(essentials[0].sum()/satisfied.sum())

print(np.array(ratios[300:]).mean())
