

import galois
import numpy as np
import numba as nb
from ..utils import coordinates
# from rustworkx import PyGraph
# from rustworkx import cartesian_product as product, graph_adjacency_matrix as gam

from .Cell import Cell, Vertex


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
                conditions (i.e. are we operating on a torus)?
        """
        # Defaults.
        self.corners = corners
        self.dimension = dimension if dimension and dimension <= len(corners) else len(corners)
        self.field = galois.GF(field)
        self.periodicBoundaryConditions = periodicBoundaryConditions

        # Construct the lattice.
        self.skeleta = {}
        self._topDown()

    
    def _topDown(self):
        """
        Construct the lattice in a top-down fashion. Recalling that a cubical
        n-cell has 2n (n-1) faces, we can determine the point at the center of
        each n-cell, figure out neighbors, and build out the remainder of the
        lattice. Moreover, we *should* be able to create a graph from this, so
        we can bump up the speeds of things.
        """
        # For each coordinate, create a cell of the appropriate dimension.
        bigSkeleton = [Cell(dimension=self.dimension, corner=c) for c in coordinates(self.corners)]
        self.skeleta[self.dimension] = set(bigSkeleton)

        # Descend into the faces of each Cell, constructing the skeleta. We make
        # sets of these things because we want to minimize the amount of data
        # eaten up by this structure.
        

