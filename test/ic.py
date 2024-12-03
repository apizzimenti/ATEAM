
import numpy as np
import math
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from ateam.structures import Lattice
from ateam.models import InvadedCluster
from ateam import Chain, Tape, _version
import sys



# Construct lattice object.
field = 3
L = Lattice().fromCorners([3,3,3,3], field=field)

# Set up Model and Chain.
homology = 2
rank = math.comb(len(L.corners), homology)
SW = InvadedCluster(L, homology=homology)
N = 1000
M = Chain(SW, steps=N)

def chain():
	for (spins, essentials, satisfied) in M.progress():
		pass

