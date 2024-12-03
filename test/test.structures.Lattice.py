
from ateam.structures import Lattice
from ateam.models import HomologicalPercolation
from ateam import Chain
import numpy as np

homology = 2
L = Lattice().fromCorners([2,2,2,2], field=3)
L.toFile("local.json")
