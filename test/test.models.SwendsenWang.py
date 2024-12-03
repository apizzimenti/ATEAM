
from potts.structures import Lattice
from potts.models import SwendsenWang
from potts import Chain

L = Lattice().fromCorners([3,3])
SW = SwendsenWang(L)
M = Chain(SW, steps=100)

for step in M.progress():
	print(step[0])
	print(step[1])
