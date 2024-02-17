
import numpy as np
from rustworkx import connected_components as connectedComponents

from ..structures import GraphLattice
from ..stats import constant
from .Model import Model


class GraphSwendsonWang(Model):
    def __init__(
            self, L: GraphLattice, temperatureFunction=constant(-0.6), initial=None
        ):
        """
        Initializes a Swendson-Wang evolution on the Potts model.

        Args:
            temperature (Callable): A function taking an integer argument which
                returns a real-valued temperature. Defaults to the constant
                schedule.
            testing (Bool): Are we testing?
        """
        self.lattice = L
        self.temperatureFunction = temperatureFunction

        self.state = initial if initial else self.initial()
        self.spins = { face: self.state[face.index] for face in self.lattice.faces }
        self.occupied = { cube: 1 for cube in self.lattice.cubes }
    

    def proposal(self, time):
        """
        Proposal scheme for the Swendson-Wang evolution on the Potts model. 

        Args:
            chain (Chain): Chain object which contains all the information we
                need.

        Returns:
            A proposed state.
        """
        # Compute the probability of choosing any individual edge in the graph.
        self.temperature = self.temperatureFunction(time)
        p = 1-np.exp(self.temperature)

        # Get the graph and choose which edges to include.
        G = self.lattice.graph
        include = []

        for edge in G.edges():
            u, v = edge.at

            if u.spin == v.spin:
                q = np.random.uniform()
                
                if q < p:
                    include.append((u.index, v.index))
                    edge.spin = 1
                else:
                    edge.spin = 0

            else:
                edge.spin = 0

        # Do stuff to the subgraph.
        subgraph = G.edge_subgraph(include)
        components = connectedComponents(subgraph)
        
        # For each vertex in each component, assign a spin.
        for vertexset in components:
            q = np.random.randint(low=0, high=self.lattice.field.order)
            for index in vertexset: G[index].spin = q

        return [v.spin for v in G.nodes()]


    def initial(self):
        """
        Generates a random initial state.

        Args:
            lattice (Lattice): The integer lattice we're experimenting on.
            distribution (Callable): A distribution from which we'll select;
                this distribution must be _discrete_ and take _two integer inputs_
                representing the lower and higher endpoints of the interval over
                the integers.

        Returns:
            A cocycle with initial states.
        """
        G = self.lattice.graph
        vertices = G.nodes()

        for vertex in vertices: vertex.spin = int(np.random.randint(0, self.lattice.field.order))
        return [v.spin for v in vertices]


    def assign(self, cocycle):
        for index, spin in enumerate(cocycle): self.lattice.graph[index].spin = spin
        for edge in self.lattice.graph.edges():
            u, v = edge.at
            edge.spin = (0 if u.spin != v.spin else 1)
