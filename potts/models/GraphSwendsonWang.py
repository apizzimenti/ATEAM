
import numpy as np
from rustworkx import connected_components as connectedComponents
from pathlib import Path

from ..stats import constant, uniform
from .Model import Model


class GraphSwendsonWang(Model):
    def __init__(
            self, temperatureFunction=constant(-0.6), testing=False
        ):
        """
        Initializes a Swendson-Wang evolution on the Potts model.

        Args:
            temperature (Callable): A function taking an integer argument which
                returns a real-valued temperature. Defaults to the constant
                schedule.
            testing (Bool): Are we testing?
        """
        self.temperatureFunction = temperatureFunction
        self.testing = testing
        self.log = ""

        self._directorySetup()
    

    def _directorySetup(self):
        if self.testing:
            # Creates an output directory if none exists.
            outroot = Path("./output/")
            if not outroot.exists():
                outroot.mkdir()
                (outroot/"figures").mkdir()
                (outroot/"matrices").mkdir()
    

    def proposal(self, chain):
        """
        Proposal scheme for the Swendson-Wang evolution on the Potts model. 

        Args:
            chain (Chain): Chain object which contains all the information we
                need.

        Returns:
            A proposed state.
        """
        # Compute the probability of choosing any individual edge in the graph.
        self.temperature = self.temperatureFunction(chain.step)
        p = 1-np.exp(self.temperature)

        # Get the graph and choose which edges to include.
        G = chain.lattice.graph
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
            q = np.random.randint(low=0, high=chain.lattice.field.order)
            for index in vertexset: G[index].spin = q

        return [v.spin for v in G.nodes()]


    def initial(self, lattice, distribution=uniform):
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
        G = lattice.graph
        vertices = G.nodes()

        for vertex in vertices: vertex.spin = int(distribution(0, lattice.field.order))
        return [v.spin for v in vertices]


    def energy(self, lattice, state):
        """
        Computes the Hamiltonian (energy) of the lattice in its current state.

        Args:
            lattice (Lattice): The lattice we're working over.
            state (list): Optional argument for computing the energy of an
                arbitrary state instead of the current one in the chain.

        Returns:
            Integer representing the Hamiltonian.
        """
        G = lattice.graph

        s = 0
        for edge in G.edges():
            u, v = edge.at
            s += (1 if state[u.index] == state[v.index] else 0)

        return -s

