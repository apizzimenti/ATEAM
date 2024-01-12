
import galois
import numpy as np
import scipy as scp
from itertools import combinations as combs
from functools import reduce

from ..utils import coordinates, binaryEncode, elementwiseAdd
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

        # Pre-compute stuff that's pre-computable.
        self._initializeComputables()

        # Construct the lattice.
        self.skeleta = { k: {} for k in range(self.dimension+1)}
        self._bottomUp()

        # Construct the boundary matrix.
        self._constructBoundaryMatrix()


    def _initializeComputables(self):
        # Construct the integer basis for a *single* cube.
        self.basis = frozenset(2**k for k in range(self.dimension))

        # Find the power set of the basis elements.
        _power = [[frozenset(c) for c in combs(self.basis, d)] for d in range(1, self.dimension+1)]
        self.power = reduce(lambda A,B: A+B, _power)

        # For each subset in the power set of basis elements, compute the shifts
        # we'll apply to the subfaces admitted by the subset. We store these in
        # a dictionary, since we spend so much time re-computing them.
        self.shifts = {}
        self.subsetBasisElements = {}
        self.innerShifts = {}

        for subset in self.power:
            # First, find the amounts we'll shift each face by. Leave out the
            # identity shift (i.e. the 0-combination).
            __shifts = [[frozenset(c) for c in combs(self.basis-subset, d)] for d in range(1, self.dimension+1)]
            _shifts = reduce(lambda A,B: A+B, __shifts)
            self.shifts[subset] = [sum(s) for s in _shifts]

            dimension = len(subset)
            subsetBasisElements = [frozenset(c) for c in combs(subset, dimension-1)]
            self.subsetBasisElements[subset] = subsetBasisElements
            self.innerShifts[subset] = [set(subset-c).pop() for c in subsetBasisElements]


    def IntegerHammingCube(self):
        """
        Constructs an integer-cornered Hamming cube.

        Returns:
            A dictionary containing the faces of an integer-valued Hamming cube,
            indexed by dimension.
        """
        # Skeleton structure.
        skeleta = { k: {} for k in range(self.dimension+1) }

        # Create vertices.
        vertices = [IntegerReducedCell(k, True) for k in range(2**self.dimension)]
        skeleta[0] = { c.encoding: c for c in vertices }

        # Initialize basis elements. We get creative here!
        basisFaces = { k: vertices[k] for k in range(2**self.dimension) }

        for subset in self.power:
            dimension = len(subset)

            # Check whether we even need to construct an element of this kind.
            if dimension > self.dimension: break

            # First, find the amounts we'll shift each face by. Leave out the
            # identity shift (i.e. the 0-combination).
            shifts = self.shifts[subset]

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
            subsetBasisElements = self.subsetBasisElements[subset]
            subsetInnerShifts = self.innerShifts[subset]

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

        # Return the skeleton.
        return skeleta

    
    def CoordinateHammingCube(self):
        """
        Constructs a coordinate Hamming cube, given an integer-valued one.

        Returns:
            A coordinate-valued Hamming cube.
        """
        # Get the skeleta and convert all the faces to coordinate-valued ones.
        integerSkeleta = self.IntegerHammingCube()
        coordinateSkeleta = { d: {} for d in range(self.dimension+1) }

        integerCube = list(integerSkeleta[self.dimension].values())[0]
        N = self.dimension

        def descendIntoFace(cube):
            if cube.dimension == 0:
                # We translate everything and *reverse* the encoding, otherwise
                # the coordinate orders don't match up.
                translate = ReducedCell(
                    tuple(reversed(binaryEncode(cube.vertices[0], len(self.corners)))),
                    vertex=True
                )
            else:
                translate = ReducedCell([descendIntoFace(f) for f in cube.faces])

            # Check if we've already created this face. If not, chuck it into the
            # dictionary, and return the translation (whether it's been created
            # already or not).
            if not coordinateSkeleta.get(translate.encoding, False):
                coordinateSkeleta[translate.dimension][translate.encoding] = translate

            return coordinateSkeleta[translate.dimension][translate.encoding]
        
        # Convert the integer-cornered Hamming cube, adding its components to
        # the skeleton's lattice.
        coordinateCube = descendIntoFace(integerCube)
        return coordinateCube
    

    @staticmethod
    def translateCell(cell, corner):
        """
        Translates this ReducedCell so it is anchored at the given coordinate
        `corner`.

        Args:
            corner (tuple): Destination coordinate for the bottom-left corner of
                this Cell. 
        """
        def descendToVertices(cube):
            # If the cube is a vertex, we modify the coordinates of the vertex.
            # Otherwise, we just apply this map to all the faces of the cube.
            if cube.dimension == 0:
                cube.vertices = [tuple(elementwiseAdd(cube.vertices[0], corner))]
                # cube.vertices = [tuple(sum(z) for z in zip(cube.vertices[0], corner))]
            else:
                for face in cube.faces: descendToVertices(face)

            # Re-encode and return.
            cube.reEncode()
            return cube
        
        return descendToVertices(cell)

    
    def replaceFaces(self, cube):
        """
        Checks whether the faces of a given cube already exist in this Lattice;
        if they do, replace them; if not, add them.
        """
        def depthFirst(face):
            if not self.skeleta[face.dimension].get(face.encoding, False):
                face.faces = set(depthFirst(f) for f in face.faces)
                self.skeleta[face.dimension][face.encoding] = face
            else:
                face = self.skeleta[face.dimension][face.encoding]

            return face
        
        return depthFirst(cube)


    def _bottomUp(self):
        """
        Construct the lattice in a bottom-up fashion.
        """
        # Find all the vertices (coordinate points) from the lattice, and build
        # from there.
        vertices = coordinates(self.corners)
        
        # Given this prototype, translate it to create the remaining objects in
        # the lattice. If we encounter something we haven't seen before, we add
        # it to the skeleta; if we encounter something we *have* seen before, we
        # replace it with what already exists.
        for vertex in vertices:
            cube = self.CoordinateHammingCube()
            translated = self.translateCell(cube, vertex)
            existing = self.replaceFaces(translated)


    def _constructBoundaryMatrix(self):
        """
        Constructs the (co)boundary matrix for this cubical lattice.
        """
        faces = list(self.skeleta[self.dimension-1].values())
        cubes = list(self.skeleta[self.dimension].values())
        
        # Create indexes, mapping encodings to row or column indices.
        faceIndex = dict(zip(faces, range(len(faces))))
        cubeIndex = dict(zip(cubes, range(len(cubes))))

        # Construct a compressed sparse column (CSC) matrix of the required size.
        faceCount = len(faces)
        cubeCount = len(cubes)
        B = scp.sparse.csc_array((faceCount, cubeCount), dtype=int)

        # Iterate over the cubes, accounting for faces.
        for cube in cubes:
            j = cubeIndex[cube]
            for face in cube.faces:
                i = faceIndex[face]
                B[i,j] = 1

        self.boundary = B
        self.coboundary = B.transpose()
