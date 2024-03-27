
import numpy as np
import math

from ..stats import constant, uniform
from .Model import Model


class GraphIsing(Model):
    name = "GraphIsing"
    
    def __init__(self, temperatureFunction=constant(-math.log(math.sqrt(2)/(1+math.sqrt(2))))):
        """
        Initializes a Swendson-Wang evolution on the Potts model.

        Args:
            temperature (Callable): A function taking an integer argument which
                returns a real-valued temperature. Defaults to the constant
                schedule.
            testing (Bool): Are we testing?
        """
        self.temperatureFunction = temperatureFunction
    

    def proposal(self, chain):
        """
        Proposal scheme for the Swendson-Wang evolution on the Potts model. 

        Args:
            chain (Chain): Chain object which contains all the information we
                need.

        Returns:
            A proposed state.
        """
        # Uniformly randomly sample a vertex from the graph and change its spin.
        vertices = chain.lattice.structure[0]
        v = vertices[np.random.randint(low=0, high=len(vertices))]
        v.spin = (v.spin+1)%2

        # Return the sequence of spins.
        return [v.spin for v in vertices]


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

        for vertex in vertices: vertex.spin = np.random.randint(low=0, high=2)
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
            s += state[u.index]*state[v.index]

        return -s

