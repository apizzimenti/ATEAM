
from ateam.structures import Lattice
from ateam.models import Glauber
from ateam import Chain

L = Lattice().fromCorners([3,3,3,3], dimension=2)
SW = Glauber(L)
M = Chain(SW, steps=100)

for step in M.progress():
	print(len(L.boundary[1]), len(step[0]))
	print(len(L.boundary[2]), len(step[1]))
