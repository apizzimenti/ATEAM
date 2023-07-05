
import numpy as np
import pandas as pd
import numba as nb
from gerrytools.plotting import districtr
from rustworkx import minimum_spanning_tree as msp
from rustworkx import minimum_spanning_edges as mse
from rustworkx import connected_components as connectedComponents
from pathlib import Path

import matplotlib.pyplot as plt
from rustworkx.visualization import mpl_draw as draw

from ..stats import constant, always, uniform
from .Model import Model


class GraphSwendsonWang(Model):
    def __init__(
            self, temperature=constant(-0.6), accept=always(), testing=False
        ):
        """
        Initializes a Swendson-Wang evolution on the Potts model.

        Args:
            temperature (Callable): A function taking an integer argument which
                returns a real-valued temperature. Defaults to the constant
                schedule.
            accept (Callable): A function which takes a Chain as input and returns
                a boolean.
        """
        self.temperature = temperature
        self.accept = accept
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
        p = 1-np.exp(self.temperature(chain.step))

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


    def energy(self, lattice):
        """
        Computes the Hamiltonian (energy) of the lattice in its current state.

        Args:
            lattice (Lattice): Lattice over which we're working.

        Returns:
            Integer representing the Hamiltonian.
        """
        s = 0
        for edge in lattice.graph.edges():
            u, v = edge.at
            s += (1 if u.spin == v.spin else 0)

        return -s

