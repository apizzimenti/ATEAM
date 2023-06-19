
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
        if self.dimension < 1:
            return str(self.coordinates)
        else:
            return ""


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

    def __init__(self, corners):
        """
        Instantiates an integer lattice.

        Args:
            corners (list): An iterable collection of "corners" for the lattice:
                for example, the argument `[1,1,1]` gives the unit cube in a
                3-dimensional integer lattice. More generally, an argument of
                `[c1, ..., cn]` admits an n-dimensional integer lattice with the
                ith copy of Z bounded below by 0 and above by ci (inclusive).
        """
        self.corners = corners
        self.dimension = len(corners)   

        # Procedurally generate all of the coordinates in the lattice bounded
        # by the given corners. First, `basis` will consist of the integer-valued
        # vectors defining each axis; then, these vectors will be combined to
        # generate the remaining coordinates.
        allVectorsByAxis = []

        for axis in range(self.dimension):
            vectorsForAxis = []

            for k in range(corners[axis] + 1):
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
        self.structure = np.ndarray(tuple([c+1 for c in corners]), dtype=Simplex)
        for simplex in self.coordinates:
            self.structure[tuple(simplex.coordinates)] = simplex
