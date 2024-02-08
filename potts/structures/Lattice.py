
import galois
import numpy as np
import json
import ast
from rustworkx import PyGraph, cycle_basis as basis
from itertools import combinations as combs, product
from functools import reduce

from ..arithmetic import coordinates, binaryEncode, elementwiseAdd
from .Cell import Cell, IntegerCell, ReducedCell


class Index:
    cubes = {}
    faces = {}


class Lattice:
    """
    A class which, post-construction of the giganto-lattice, keeps only records
    of the encodings of cubes and faces. This format is much smaller and far
    easier to serialize and store.
    """
    def __init__(self): pass
    
    def fromCorners(self, corners, field=2, dimension=None, periodicBoundaryConditions=True):
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
        self.dimension = dimension if dimension else len(corners)
        self.periodicBoundaryConditions = periodicBoundaryConditions

        # Construct an initial LargeLattice.
        _L = LargeLattice(
            corners=corners, dimension=self.dimension,
            periodicBoundaryConditions=self.periodicBoundaryConditions
        )

        # Construct the finite field object.
        self.field = galois.GF(field)

        # Reduced cells; they carry only the encoding and, if they're cubes, the
        # list of faces making them up.
        self.faces = list(sorted(ReducedCell(f.encoding[:-1]) for f in _L.skeleta[self.dimension-1].values()))
        _faceLookup = { f.encoding: f for f in self.faces}

        self.cubes = list(sorted(
            ReducedCell(c.encoding[:-1], faces=[_faceLookup[f.encoding[:-1]] for f in c.faces])
            for c in _L.skeleta[self.dimension].values()
        ))

        # Construct indices, boundary matrix, and graph.
        self._index()
        self._constructBoundaryMatrix()
        self._constructGraph()

        # Force garbage collection on the Lattice and the lookup.
        del _faceLookup
        del _L


    def toFile(self, fp):
        """
        JSON-serializes this object and writes it to file so we can reconstruct
        it later.

        Args:
            fp: File object; must be in write mode.
        """
        json.dump({
                "faces": { k: str(f.encoding) for f, k in self.index.faces.items() },
                "cubes": {
                    str(cube.encoding): [self.index.faces[f] for f in cube.faces]
                    for cube in self.cubes
                },
                "field": self.field.order,
                "dimension": self.dimension,
                "periodicBoundaryConditions": int(self.periodicBoundaryConditions)
            }, fp
        )
        

    def fromFile(self, fp):
        """
        Reconstructs a serialized Lattice.

        Args:
            fp: Python file object; must be in read mode.
        """
        # Read file into memory.
        serialized = json.load(fp)

        # Set field and dimension.
        self.field = galois.GF(int(serialized["field"]))
        self.dimension = int(serialized["dimension"])
        self.periodicBoundaryConditions = bool(serialized["periodicBoundaryConditions"])

        # Create faces and cubes.
        faces = {
            int(k): ReducedCell(ast.literal_eval(f)) for k, f in serialized["faces"].items()
        }

        self.faces = list(sorted(faces.values()))
        self.cubes = list(sorted(
            ReducedCell(ast.literal_eval(c), faces=list(sorted(faces[k] for k in indices)))
            for c, indices in serialized["cubes"].items()
        ))

        # Construct indices, boundary matrix, and graph.
        self._index()
        self._constructBoundaryMatrix()
        self._constructGraph()

        del serialized
        del faces

    
    def _index(self):
        """
        Creates an internal index for the faces and cubes of the Lattice.
        """
        # Create an indexing subobject.
        self.index = Index()
        self.index.cubes = { c: k for k, c in enumerate(self.cubes) }
        self.index.faces = { f: k for k, f in enumerate(self.faces) }

    
    def _constructBoundaryMatrix(self):
        """
        Constructs the (co)boundary matrix for this cubical lattice.
        """
        faceCount = len(self.faces)
        cubeCount = len(self.cubes)
        
        # Construct the boundary matrix.
        B = np.zeros((faceCount, cubeCount), dtype=int)
        
        for cube in self.cubes:
            j = self.index.cubes[cube]
            for face in cube.faces:
                i = self.index.faces[face]
                B[i,j] = 1
        
        self.boundary = self.field(B)
        self.coboundary = self.field(B.transpose())

        self.subboundary = None
        self.subcoboundary = None


    def _constructGraph(self):
        """
        Constructs and maintains a graph, where faces are vertices and cubes are
        edges; two cubes are adjacent if and only if they share a face.

        TODO: we're currently doing this naïvely, and it'd be better to do it in
        a more... efficient fashion. (Really we're just hitting the array twice.)
        """
        # Catalogue to which cubes each face belongs.
        adjacencies = { f: set() for f in self.faces }

        for cube in self.cubes:
            for face in cube.faces:
                adjacencies[face].add(cube)

        # Construct a graph (on the cubes) from the adjacencies!
        A = np.zeros((len(self.cubes), len(self.cubes)))

        for face, cubes in adjacencies.items():
            # Get the indices for the cubes; take the product over the list and
            # add ones.
            indices = product([self.index.cubes[c] for c in cubes], repeat=2)
            indices = [(i,j) for i,j in indices if i!=j]
            for coord in indices: A[coord] = 1

        # Construct a graph from the matrix A. The indices of the cubes *should*
        # map onto the indices of the graph. We also attach a pointer to each
        # Cube to each vertex in this graph.
        self.graph = PyGraph().from_adjacency_matrix(A)
        for cube, index in self.index.cubes.items(): self.graph[index] = cube

        # Create a subgraph object for later.
        self.subgraph = None

        

class LargeLattice:
    def __init__(
            self, corners, dimension=None, periodicBoundaryConditions=True
        ):
        """
        Creates a cell complex with the given corners made of cells of the
        provided dimension.

        Args:
            corners (list): Corners of the lattice; determines the maximal
                cell dimension.
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
        self.dimension = len(self.corners)
        self.periodicBoundaryConditions = periodicBoundaryConditions

        # Pre-compute stuff that's pre-computable.
        self._initializeComputables()

        # Construct the lattice.
        self.skeleta = { k: {} for k in range(self.dimension+1)}
        self._bottomUp()


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
        vertices = [IntegerCell(k, True) for k in range(2**(self.dimension))]
        skeleta[0] = { c.encoding: c for c in vertices }

        # Initialize basis elements. We get creative here!
        basisFaces = { k: v for k, v in enumerate(vertices) }

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
                cube = IntegerCell([u, v])
                basisFaces[subset] = cube
                skeleta[1][cube.encoding] = cube

                # For each set of shifts, find the encoding and check whether
                # this cube has already been found. This *shouldn't* be the case,
                # but we're doing it just in case.
                for shift in shifts:
                    shiftedEncoding = cube.shiftedEncoding(shift)
                    w, x, _ = shiftedEncoding
                    if not skeleta[1].get(shiftedEncoding, False):
                        shiftedCube = IntegerCell([basisFaces[w], basisFaces[x]])
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
            newBasisCube = IntegerCell(shiftedFaces)
            basisFaces[subset] = newBasisCube
            skeleta[dimension][newBasisCube.encoding] = newBasisCube

            # Now, shift each of these faces by the remaining basis elements to
            # get our other cubes.
            for shift in shifts:
                shiftedCubeEncoding = newBasisCube.shiftedEncoding(shift)

                if not skeleta[dimension].get(shiftedCubeEncoding, False):
                    shiftedFaceEncodings = [face.shiftedEncoding(shift) for face in newBasisCube.faces]
                    shiftedFaces = [skeleta[dimension-1][e] for e in shiftedFaceEncodings]
                    shiftedCube = IntegerCell(shiftedFaces)

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

        cube = list(integerSkeleta[self.dimension].values())[0]
        N = len(self.corners)

        def descendIntoFace(cube):
            if cube.dimension == 0:
                # We translate everything and *reverse* the encoding, otherwise
                # the coordinate orders don't match up.
                translate = Cell(
                    tuple(reversed(binaryEncode(cube.vertices[0], N))),
                    vertex=True
                )
            else:
                translate = Cell([descendIntoFace(f) for f in cube.faces])

            # Check if we've already created this face. If not, chuck it into the
            # dictionary, and return the translation (whether it's been created
            # already or not).
            if not coordinateSkeleta.get(translate.encoding, False):
                coordinateSkeleta[translate.dimension][translate.encoding] = translate

            return coordinateSkeleta[translate.dimension][translate.encoding]
        
        # Convert the integer-cornered Hamming cube, adding its components to
        # the skeleton's lattice.
        return descendIntoFace(cube)
    

    @staticmethod
    def translateCell(cell, anchor):
        """
        Translates this ReducedCell so it is anchored at the given coordinate
        `corner`.

        Args:
            cell (Cell): Cell into whose faces we'll be delving.
            anchor (tuple): Destination coordinate for the bottom-left corner of
                this Cell.
        """
        def descendToVertices(cube):
            # If the cube is a vertex, we modify the coordinates of the vertex.
            # Otherwise, we just apply this map to all the faces of the cube.
            if cube.dimension == 0:
                cube.vertices = [
                    tuple(elementwiseAdd(cube.vertices[0], anchor))
                ]
            else:
                for face in cube.faces: descendToVertices(face)

            # Re-encode and return.
            cube.reEncode()
            return cube
        
        return descendToVertices(cell)

    
    def replaceFaces(self, cube):
        """
        Checks whether the faces of a given cube already exist in this Lattice;
        if they do, replace them; if not, add them. If we've asserted periodic
        boundary conditions, we replace any boundary face (i.e. any face which
        is part of exactly one cube) with sthe face at the opposite end of the
        Lattice.
        """
        def depthFirst(face):
            # Check for boundary conditions, replacing vertices when necessary.
            # This will propagate up to the higher-dimensional cubes.
            if face.dimension == 0 and self.periodicBoundaryConditions:
                # Find *all* indices where the location of the vertex is incident
                # to the outer face.
                location = face.encoding[:-1][0]
                indices = [
                    i for i, (j, k) in enumerate(zip(location, self.corners))
                    if j == k
                ]

                # If there are any indices in the list, we need to replace the
                # vertex!
                if len(indices) > 0:
                    modifiedEncoding = (
                        tuple(location[i] if i not in indices else 0 for i in range(len(self.corners))),
                        0
                    )
                    face = self.skeleta[0][modifiedEncoding]
        
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
            self.replaceFaces(translated)

        if self.periodicBoundaryConditions:
            # We now scan through the skeleta, deleting faces with improper vertices;
            # we only do this if the user has specified periodic boundary conditions
            # on the Lattice.
            proper = set(f[0] for f in self.skeleta[0].keys())
            
            for dimension in range(1, self.dimension+1):
                good = {}

                for cube in self.skeleta[dimension].values():
                    cube.reEncode()
                    if set(cube.encoding[:-1]).issubset(proper): good[cube.encoding] = cube

                self.skeleta[dimension] = good

        # Set the "cubes" and "faces" of the lattice.
        self.cubes = list(self.skeleta[self.dimension].values())
        self.faces = list(self.skeleta[self.dimension-1].values())
