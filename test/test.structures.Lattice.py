
from potts.structures import Lattice
from potts.models import HomologicalPercolation
from potts import Chain
import numpy as np

homology = 2
L = Lattice().fromCorners([2,2,2,2], field=3)
L.toFile("local.json")
