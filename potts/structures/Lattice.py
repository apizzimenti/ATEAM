
import galois
import scipy as scp
import numpy as np
from numba import jit, njit
from math import comb
from itertools import combinations

from ..utils import coordinates, binaryEncode, subtractMany, increment, binaryUnencode
from .Cell import ReducedCell


class Lattice:
    def __init__(
            self, corners, field=2, dimension=None, periodicBoundaryConditions=True
        ):
        """
        Creates a cell complex with the given corners made of cells of the
        provided dimension.

        Args:
            corners (list): Corners of the lattice; determines the maximal
                cell dimension.
            field (int): Characteristic of finite field from which cells take
                coefficients.
            dimension (int): Maximal cell dimension; if this argument is larger
                than that permitted by the underlying cell structure, this is
                re-set to the maximal legal dimension. Defaults to 1, so the
                lattice is constructed only of vertices (0-cells) and edges
                (1-cells).
            periodicBoundaryConditions (bool): Do we use periodic boundary
                conditions (i.e. are we making a torus)?
        """
        # Defaults.
        self.corners = corners
        self.dimension = dimension if dimension and dimension <= len(corners) else len(corners)
        self.field = galois.GF(field)
        self.periodicBoundaryConditions = periodicBoundaryConditions

        # Construct the lattice.
        self.skeleta = {}
        self._bottomUp()
        
        # Construct the boundary matrix.
        # self._boundaryMatrix()

    def _bottomUp(self):
        """
        Construct the lattice in a bottom-up fashion.
        """
        # Find all the vertices (coordinate points) from the lattice, and build
        # from there. Eventually we'll have to take periodic boundary conditions
        # into account.
        vertices = coordinates(self.corners)
        
        # Construct 0-cells (vertices).
        self.skeleta[0] = {
            v: ReducedCell(v, vertex=True) for v in vertices
        }

        # Construct higher-dimensional cells.
        self._higherDimensionalCells()

    
    def _higherDimensionalCells(self):
        vertices = self.skeleta[0].values()

        for dimension in range(1, self.dimension+1):
            for vertex in vertices:
                # Get the "anchor" vertex; if this vertex is a boundary vertex (i.e.
                # is too close to the edge in any coordinate) we ignore and continue.
                # For example, a 1x1x1 integer lattice should have *precisely* one
                # cube, six squares, 12 edges, and eight vertices.
                _anchor = vertex.vertices[0]
                if not all(_anchor[t] < self.corners[t] for t in range(len(self.corners))): continue

                # Get the integer value of the coordinates. For each dimension from
                # 1 to the dimension of the complex (or the user-specified dimension,
                # whichever is lower) we find the appropriately-sized hypercubes and
                # construct cells. We do so by figuring out how many different ways
                # there are to combine our basis elements!
                anchor = binaryEncode(_anchor)
                basis = [2**k for k in range(len(_anchor))]

                # Figure out the ways we can combine elements to determine faces.
                linearCombinations = combinations(basis, r=dimension)
                print(anchor)
                print(list(linearCombinations))
                print()
