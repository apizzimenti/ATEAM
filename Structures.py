
from functools import reduce
import numpy as np

class Simplex:
    def __init__(self, coordinates, orientation=1, index=0):
        """
        Initializes a simplex of any dimension.

        Args:
            coordinates (np.array): Will be integers (0-simplex) or Simplex entries
                (k-simplex, k > 0).
            orientation (int): Which way do we traverse the entries in `coordinates`?
            index (int): When constructing matrices, the index to which this simplex
                corresponds.
        """
        # This determines which kind of simplex we're dealing with: if it's a
        # 0-dimensional simplex, then `coordinates` will have integer entries; if
        # it's any higher-dimensional one, its entries will be instances of Simplex.
        if type(coordinates[0]) in {int, np.int32, np.int64}:
            self.dimension = 0
        else:
            self.dimension = coordinates[0].dimension + 1

        self.coordinates = coordinates
        self.orientation = orientation
        self.index = index

    def __repr__(self): return self.__str__()

    def __str__(self):
        """
        Stringification.
        """
        if self.dimension < 1:
            return str(tuple(self.coordinates))
        else:
            return f"{self.dimension}-dimensional simplex at index {self.index}"
        
    def __eq__(self, other):
        """
        Testing for equality; this asks whether the coordinates of each simplex
        are the same.
        """
        return (self.coordinates == other.coordinates).all()
    
    def __lt__(self, other):
        pass

    def __hash__(self):
        """
        Returns a unique property from each simplex.
        """
        return hash(tuple(self.coordinates))


class Lattice:
    """
    Encodes a d-dimensional integer lattice.

    Attributes:
        corners: Bounds of the integer lattice; equivalent to boundary constraints.
        dimension: Dimension of the integer lattice as determined by the number
            of corners (boundary constraints).
        coordinates: In-order integer-valued vectors representing the vertices of
            the lattice.
    """

    def __init__(self, corners, boundaryDimension=1, maxDimension=None):
        """
        Instantiates an integer lattice.

        Args:
            corners (list): An iterable collection of "corners" for the lattice:
                for example, the argument `[1,1,1]` gives the unit cube in a
                3-dimensional integer lattice. More generally, an argument of
                `[c1, ..., cn]` admits an n-dimensional integer lattice with the
                ith copy of Z bounded below by 0 and above by ci (inclusive).
            boundaryDimension (int): Specifies the dimension for which we construct
                the boundary/coboundary operator matrices; this defaults to 1.
            maxDimension (int): The maximum dimension of simplex constructed;
                if nothing is passed, simplices of all dimensions (from 0 to
                the dimension of the lattice) are constructed.
        """
        # Assign corners and dimensionality.
        self.corners = corners
        self.dimension = len(corners)

        # Construct the lattice.
        self._coordinates()

        # Construct higher-dimensional simplices.
        self._constructHigherDimensionalSimplices(dimension=maxDimension)

        # Construct the boundary operator matrix with the provided dimension.
        self.boundaryOperator(dimension=boundaryDimension)


    def _coordinates(self):
        """
        Private method for generating the coordinates (vertices) in the lattice.
        """
        # Procedurally generate all of the coordinates in the lattice bounded
        # by the given corners. First, `basis` will consist of the integer-valued
        # vectors defining each axis; then, these vectors will be combined to
        # generate the remaining coordinates.
        allVectorsByAxis = []

        for axis in range(self.dimension):
            vectorsForAxis = []

            for k in range(self.corners[axis] + 1):
                v = [0]*self.dimension
                v[axis] = k
                vectorsForAxis.append(v)
            
            allVectorsByAxis.append(vectorsForAxis)

        # Initialize the coordinates to be only those in the first axis. Then,
        # we attack this problem in a dynamic-programming way: because we don't
        # know the number of loops required ahead-of-time, we'll progressively
        # pair each vector on the first axis with each vector on the second axis;
        # then, we'll add each vector in the third axis to each of those pairs;
        # we'll continue in this way until we end up with all possible combinations
        # of integer-scaled vectors. Afterwarsd, all we have to do is elementwise-add
        # the vectors in each combination, and we'll have the coordinates outright.
        combinations = [[b] for b in allVectorsByAxis[0]]
        
        for axisIndex in range(1, len(allVectorsByAxis)):
            intermediate = []
            
            for rightAdjoin in allVectorsByAxis[axisIndex]:
                for combination in combinations:
                    ccombination = list(combination)
                    ccombination.append(rightAdjoin)
                    intermediate.append(ccombination)
            
            combinations = intermediate

        # Now that we've computed the combinations, we need only add.
        coordinates = [
            reduce(lambda l, r: np.array([l[i] + r[i] for i in range(len(l))]), combination)
            for combination in combinations
        ]

        # Turn the coordinates into simplices, and then turn the list of
        # coordinates into an n-dimensional array.
        self.coordinates = np.array([
            Simplex(coordinate, index=index) for index, coordinate in enumerate(coordinates)
        ])

        # Create the internal ndarray structure so indexing is easy; this fixes
        # the amount of memory we'll use on the Lattice.
        zeroSimplices = np.ndarray(tuple([c+1 for c in self.corners]), dtype=Simplex)
        for simplex in self.coordinates:
            zeroSimplices[tuple(simplex.coordinates)] = simplex

        self.vertices = zeroSimplices

        # We'll store these structures in a dictionary.
        self.structure = {
            0: self.coordinates
        }


    def _constructHigherDimensionalSimplices(self, dimension=None):
        """
        Given a dimension, construct the simplices of that dimension by building
        them out of *references* to lower-dimensional simplices; for example, the
        1-simplices of the lattice are pairs of vertices with coordinates at (l1-)
        distance 1 from each other.

        Args:
            dimension (int): The dimension of the simplices we'll be identifying;
                if no dimension is passed, we construct simplices of all dimensions
                up to the dimension of the lattice.
        """
        if dimension == 1 or self.structure.get(1, True):
            self._constructEdges()
            return
    
    def _constructEdges(self):
        """
        Private method for constructing the 1-simplices (edges) of the lattice.
        """
        # Create an empty edge set.
        edges = set()

        for vertex in self.vertices.flatten():
            neighbors = set()

            # For each coordinate, we want to get all possible neighbors; since
            # this is an integer lattice, this involves grabbing all vertices at
            # (l1) distance 1 from the current coordinate.
            for axis, coordinate in enumerate(vertex.coordinates):

                # First, if our coordinate is in between the smallest and largest
                # integer values on the axis (exclusive), we just grab neighbors
                # to the left/right.
                if 0 < coordinate < self.corners[axis]: offsets = [-1, 1]

                # Next, if our coordinate is at one of the extremes, we look to
                # the left or the right to determine its neighbors.
                elif coordinate == 0: offsets = [1]
                else: offsets = [-1]

                for offset in offsets:
                    ccoordinate = list(vertex.coordinates)
                    ccoordinate[axis] = ccoordinate[axis] + offset
                    neighbors.add(self.vertices[tuple(ccoordinate)])

            # Add each neighbor in `neighbors` to the set of edges; sort the
            # neighbors in the edges (dictionary-style) first, so we don't add
            # duplicates. These don't have an order.
            for neighbor in neighbors:
                same = neighbor == vertex
                first = (neighbor, vertex) in edges
                second = (vertex, neighbor) in edges
                if not (same or first or second): edges.add((vertex, neighbor))

        # Now that we have a set of edges, we can construct 1-simplices.
        oneSimplices = []
        for index, edge in enumerate(edges):
            oneSimplices.append(Simplex(np.array(list(edge)), index=index))

        self.structure[1] = np.array(oneSimplices, dtype=Simplex)

    def boundaryOperator(self, dimension=1):
        """
        Constructs the `dimension`-dimensional boundary operator matrix (and, by
        extension, the coboundary matrix).

        Args:
            dimension (int): For which dimension are we constructing the matrix?
        """
        # Create the columns: we do so by inducing an orientation on the 
        B = np.zeros((len(self.structure[dimension-1]), len(self.structure[dimension])))

        # Impose an arbitrary orientation on the edges; we just need the vertex
        # labels to cancel.
        for edge in self.structure[1]:
            for vertex, coefficient in zip(edge.coordinates, [-1, 1]):
                B[vertex.index, edge.index] = coefficient

        self.boundary = B
        self.coboundary = B.T
                
