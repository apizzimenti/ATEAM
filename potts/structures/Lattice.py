
import galois
import json
from itertools import combinations as combs
from functools import reduce

from ..utils import coordinates, binaryEncode, subtractMany, increment, binaryUnencode
from .Cell import ReducedCell, IntegerReducedCell


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
        self.skeleta = { k: {} for k in range(self.dimension+1)}
        self._bottomUp()

    
    def IntegerHammingCube(self):
        """
        Constructs an integer-cornered Hamming cube anchored at the given corner.
        """
        # Skeleton structure.
        skeleta = { k: {} for k in range(self.dimension+1) }

        # Create vertices.
        vertices = [IntegerReducedCell(k, True) for k in range(2**self.dimension)]
        skeleta[0] = { c.encoding: c for c in vertices }

        # Initialize basis elements. We get creative here!
        basis = frozenset(2**k for k in range(self.dimension))
        basisFaces = { k: vertices[k] for k in range(2**self.dimension) }

        # Find the power set of the basis.
        _power = [[frozenset(c) for c in combs(basis, d)] for d in range(1, self.dimension+1)]
        power = reduce(lambda A,B: A+B, _power)

        for subset in power:
            dimension = len(subset)

            # First, find the amounts we'll shift each face by. Leave out the
            # identity shift (i.e. the 0-combination).
            __shifts = [[frozenset(c) for c in combs(basis-subset, d)] for d in range(1, self.dimension+1)]
            _shifts = reduce(lambda A,B: A+B, __shifts)
            shifts = [sum(s) for s in _shifts]

            # Treat edges specially since they're a bit weird.
            if dimension == 1:
                # Since the subsets here contain only base faces — and those
                # faces are vertices — we can just pop one out. We then create
                # an edge from 0 to the base face, and shift the base face
                # around.
                b = set(subset).pop()
                u = skeleta[0][(0, 0)]
                v = skeleta[0][(b, 0)]
                cube = IntegerReducedCell([u, v])
                basisFaces[subset] = cube
                skeleta[1][(u, v, 1)] = cube

                # For each set of shifts, find the encoding and check whether
                # this cube has already been found. This *shouldn't* be the case,
                # but we're doing it just in case.
                for shift in shifts:
                    shiftedEncoding = cube.shiftedEncoding(shift)
                    w, x, _ = shiftedEncoding
                    if not skeleta[1].get(shiftedEncoding, False):
                        shiftedCube = IntegerReducedCell([basisFaces[w], basisFaces[x]])
                        skeleta[1][shiftedEncoding] = shiftedCube
                        
                # Go to the next thing right away.
                continue

            # If the subset *doesn't* denote an edge, we grab the basis faces
            # and do effectively the same thing. First, get all the (dimension-1)-subsets
            # and shift each by the single remaining element.
            subsetBasisElements = [frozenset(c) for c in combs(subset, dimension-1)]
            subsetInnerShifts = [set(subset-c).pop() for c in subsetBasisElements]

            subsetBasisFaces = [basisFaces[c] for c in subsetBasisElements]
            subsetShiftedBasisFaces = []

            for face, innerShift in zip(subsetBasisFaces, subsetInnerShifts):
                shiftedBasisFaceEncoding = face.shiftedEncoding(innerShift)
                subsetShiftedBasisFaces.append(skeleta[dimension-1][shiftedBasisFaceEncoding])

            # Combine the basis faces and the shifted basis faces to get *all* the
            # faces for this hypercube, and create the object.
            shiftedFaces = subsetBasisFaces + subsetShiftedBasisFaces
            newBasisCube = IntegerReducedCell(shiftedFaces)
            basisFaces[subset] = newBasisCube
            skeleta[dimension][newBasisCube.encoding] = newBasisCube

            # Now, shift each of these faces by the remaining basis elements to
            # get our other cubes.
            for shift in shifts:
                shiftedCubeEncoding = newBasisCube.shiftedEncoding(shift)

                if not skeleta[dimension].get(shiftedCubeEncoding, False):
                    shiftedFaceEncodings = [face.shiftedEncoding(shift) for face in newBasisCube.faces]
                    shiftedFaces = [skeleta[dimension-1][e] for e in shiftedFaceEncodings]
                    shiftedCube = IntegerReducedCell(shiftedFaces)

                    skeleta[dimension][shiftedCubeEncoding] = shiftedCube


    def _bottomUp(self):
        """
        Construct the lattice in a bottom-up fashion.
        """
        # Find all the vertices (coordinate points) from the lattice, and build
        # from there. Eventually we'll have to take periodic boundary conditions
        # into account.
        vertices = coordinates(self.corners)
        self.IntegerHammingCube()