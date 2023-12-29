
import galois
import matplotlib.pyplot as plt
from functools import reduce
from rustworkx import PyGraph
from rustworkx import cartesian_product as product, graph_adjacency_matrix as gam

from .Components import Vertex, Edge


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
            self, corners, field=2, boundaryDimension=1, maxDimension=None,
            periodicBoundaryConditions=True
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
            boundaryDimension (int): Specifies the dimension for which we construct
                the boundary/coboundary operator matrices; this defaults to 1.
            maxDimension (int): The maximum dimension of cell constructed;
                if nothing is passed, cells of all dimensions (from 0 to
                the dimension of the lattice) are constructed.
        """
        # Assign corners and dimensionality.
        self.corners = corners
        self.dimension = len(corners)
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
        self.structure = {
            0: self.graph.nodes(),
            1: self.graph.edges()
        }


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
            edge.spin = 1 if u.spin == v.spin else 0
    

    @staticmethod
    def plot(
        graph, vertexStyle=dict(marker="o", markeredgewidth=0),
        edgeStyle=dict(linewidth=1/2, alpha=1/2), edgeAssignment=None,
        vertexAssignment=None, ax=None
    ):
        
        if not ax: _, ax = plt.subplots()
        
        for edge in graph.edges():
            u, v = edge.at
            ax.plot(
                *([u.at[axis], v.at[axis]] for axis in range(len(u.at))),
                color=(edgeAssignment[edge.index] if edgeAssignment else "k"),
                **edgeStyle
            )

        for vertex in graph.nodes():
            ax.plot(
                *vertex.at,
                color=(vertexAssignment[vertex.index] if vertexAssignment else "k"),
                **vertexStyle
            )

        # Set axes to be equal, turn off panes.
        ax.set_aspect("equal")
        ax.set_axis_off()
        return plt.gcf(), ax
