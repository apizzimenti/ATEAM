
from potts.structures import Lattice
from potts.viz.lattice import points



L = Lattice().fromCorners([3,3])
points(L, 1)

