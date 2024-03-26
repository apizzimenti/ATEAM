
import galois
from functools import reduce
from rustworkx import PyGraph
from rustworkx import cartesian_product as product


class Vertex:
    def __init__(self, at, spin, index):
        self.at = at
        self.spin = spin
        self.index = index

class Edge:
    def __init__(self, at, spin, index):
        self.at = at
        self.spin = spin
        self.index = index


class Index:
    cubes = {}
    faces = {}


class GraphLattice:
    """
    Encodes a d-dimensional integer lattice.

    Attributes:
        corners: Bounds of the integer lattice; equivalent to boundary constraints.
        dimension: Dimension of the integer lattice as determined by the number
            of corners (boundary constraints).
        field: Galois field from which we sample things.
        graph: Lattice as a graph.
     """

    def __init__(
            self, corners, field=2, dimension=1, periodicBoundaryConditions=True
        ):
        """
        Instantiates an integer lattice.

        Args:
            corners (list): An iterable collection of "corners" for the lattice:
                for example, the argument `[1,1,1]` gives the unit cube in a
                3-dimensional integer lattice. More generally, an argument of
                `[c1, ..., cn]` admits an n-dimensional integer lattice with the
                ith copy of Z bounded below by 0 and above by ci (inclusive).
            field (int): The finite field over which we're working; that is,
                coefficients are taken from the finite field of order `field`.
            dimension (int): Specifies the dimension for which we construct
                the boundary/coboundary operator matrices; this defaults to 1.
            periodicBoundaryConditions (bool): Do we identify the endpoints of
                our cubical lattice to turn it into a torus?
         """
        # Assign corners, dimensionality, boundary conditions.
        self.corners = corners
        self.complexDimension = len(corners)
        self.dimension = dimension if dimension else len(corners)
        self.field = galois.GF(field)
        self.periodicBoundaryConditions = periodicBoundaryConditions

        # Create an initial graph, then re-index and add stuff based on boundary
        # conditions.
        self.graph = reduce(self._reduceProduct, [self._gridFactory(c) for c in self.corners])
        for i, _ in enumerate(self.graph.nodes()): self.graph[i] = Vertex(self.graph[i], 1, i)

        for i, (j, k) in enumerate(self.graph.edge_list()):
            u, v = self.graph[j], self.graph[k]
            self.graph.update_edge_by_index(i, Edge((u, v), 1, i))

        # Construct our structure so it fits in well with the proposal.
        self.faces = self.graph.nodes()
        self.cubes = self.graph.edges()

        self.index = Index()
        self.index.faces = { }


    @staticmethod
    def _reduceProduct(G, H):
        """
        Static method for taking the cartesian product of two *graphs* and making
        the vertex labels nice.

        Args:
            G (rustworkx.PyGraph): PyGraph object to be multiplied on the right.
            H (rustworkx.PyGraph): PyGraph object to be multiplied on the left.

        Returns:
            PyGraph object which represents the product of G and H and labeled
            appropriately. 
        """
        L, _ = product(G, H)

        # Update vertices.
        for i, vertex in enumerate(L.nodes()):
            u, v = vertex
            if not(type(u) is int and type(v) is int): L[i] = (*u, v)

        # Update edges.
        for i, edge in enumerate(L.edge_list()):
            u, v = edge
            L.update_edge_by_index(i, (L[u], L[v]))

        return L
    

    def _gridFactory(self, l):
        """
        Factory for creating new grid graphs.

        Args:
            l (int): Length of the path.

        Returns:
            PyGraph object representing a path of length `l`.
        """
        path = PyGraph()

        # Add the appropriate coordinates. 
        for coordinate in range(l): path.add_node(coordinate)
        for coordinate in range(1, l): path.add_edge(coordinate-1, coordinate, (path[coordinate-1], path[coordinate]))

        # If we're using periodic boundary conditions, we're pac-manning our
        # graph: that is, we draw an edge between the first and last vertices, so
        # our graph is actually a torus.
        if self.periodicBoundaryConditions: path.add_edge(l-1, 0, (path[l-1], path[0]))
        return path
    

    def assign(self, state):
        """
        Assigns spins to vertices.

        Args:
            state (list): States to apply to vertices.
        """
        for index, spin in enumerate(state): self.graph[index].spin = spin

        for edge in self.graph.edges():
            u, v = edge.at
            edge.spin = (0 if u.spin != v.spin else 1)
