
from ateam.structures import Lattice
from ateam.models import SwendsenWang
from ateam import Chain

L = Lattice().fromCorners([3,3])
SW = SwendsenWang(L)
M = Chain(SW, steps=100)

for step in M.progress():
	print(step[0])
	print(step[1])
