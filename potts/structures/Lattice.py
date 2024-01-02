
import galois
import copy
from itertools import combinations
from functools import reduce

from ..utils import coordinates, binaryEncode, subtractMany, increment, binaryUnencode
from .Cell import ReducedCell, BinaryReducedCell


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
        self.skeleta = {}

        # Construct higher-dimensional cells.
        self._higherDimensionalCells()

    
    def _higherDimensionalCells(self):
        # First, create a prototypical Hamming cube of the specified dimension.
        prototype = self._HammingPrototype()

        # Next, based on the dimensions of the lattice, we glue *copies* of
        # the hypercube together. 


    def _HammingPrototype(self):
        """
        Creates a prototypical Hamming cube of the specified dimension.
        """
        # Basis (powers of two) for the Hamming cube.
        basis = [2**k for k in range(self.dimension)]

        # Cosntruct vertices and edges separately.
        faces = { d: {} for d in range(self.dimension+1) }
        faces[0] = { k: BinaryReducedCell(k, True)  for k in range(2**self.dimension) }
        faces[1] = { frozenset([b]): BinaryReducedCell([faces[0][0], faces[0][b]]) for b in basis }
        
        # Get all linear combinations of basis elements, and encode them as
        # frozensets so they're hashable.
        _linearCombinations = [list(combinations(basis, k)) for k in range(2, self.dimension+1)]
        linearCombinations = [frozenset(s) for s in reduce(lambda a, b: a+b, _linearCombinations)]

        for linearCombination in linearCombinations:
            # Determine the dimension of the hypercube we're constructing;
            # determine all the left-hand (i.e. incident to 0) hypercube
            # encodings of one dimension less. 
            dimension = len(linearCombination)
            subsets = [frozenset(c) for c in combinations(linearCombination, dimension-1)]
            components = []

            # For each hypercube encoding, find the corresponding ReducedCell
            # (which, by induction, exists) and "shift" the hypercube over by the
            # basis element *not* included in the encoding. For example, we get
            # the edge (2, 3) by adding 2 to each coordinate of the edge (0, 1).
            # Combine the old components with the new, and we have our new
            # component. Continue until we've exhausted the linear combinations.
            for subset in subsets:
                component = faces[dimension-1][subset]
                shift = set(linearCombination-subset).pop()
                shifted = copy.deepcopy(component)
                shifted.shift(shift)

                components.extend([component, shifted])

            # Add a new face.
            faces[dimension][linearCombination] = BinaryReducedCell(components)

        # The topmost face is the hypercube we aimed to construct: a Hamming
        # cube which contains information about *all* of its constituent faces.
        # All that's left now is to convert the vertices to coordinates, and
        # we're ready to do some shifting.
        HammingCube = copy.deepcopy(faces[self.dimension][linearCombination])
        print(HammingCube)