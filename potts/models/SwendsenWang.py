 
import numpy as np
from typing import Callable

from ..arithmetic import sampleFromKernel
from ..structures import Lattice
from ..stats import constant
from .Model import Model


class SwendsenWang(Model):
    name = "SwendsenWang"
    
    def __init__(
            self, L: Lattice, temperatureFunction: Callable=constant(-0.6),
            initial=None
        ):
        """
        Initializes Swendsen-Wang evolution on the Potts model.

        Args:
            L: The `Lattice` object on which we'll be running experiments.
            temperatureFunction (Callable): A temperature schedule function which
                takes a single positive integer argument `t`, and returns the
                scheduled temperature at time `t`.
            initial (np.ndarray): A vector of spin assignments to components.
        """
        self.lattice = L
        self.temperatureFunction = temperatureFunction
        self.faceCells = len(self.lattice.skeleta[self.lattice.dimension-1])
        self.cubeCells = len(self.lattice.skeleta[self.lattice.dimension])

        # SW defaults.
        self.state = initial if initial else self.initial()
        # self.spins = { face: self.state[self.lattice.index.faces[face]] for face in self.lattice.faces }
        self.occupied = set()


    def initial(self) -> np.array:
        """
        Computes an initial state for the model's Lattice.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        return self.lattice.field(
            np.random.randint(0, self.lattice.field.order, size=self.faceCells)
        )
    

    def proposal(self, time):
        """
        Proposal scheme for generalized Swendsen-Wang evolution on the Potts model.

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        # Compute the probability of choosing any individual cube in the complex.
        self.temperature = self.temperatureFunction(time)
        p = 1-np.exp(self.temperature)
        assert 0 <= p <= 1

        # Choose cubes to include; in effect, this just does a boatload of indexing.
        uniforms = np.random.uniform(size=self.cubeCells)
        include = (uniforms < p).nonzero()[0]
        boundary = self.lattice.boundary[self.lattice.dimension][include]
        boundaryValues = self.state[boundary-self.faceCells].sum(axis=1)
        zeros = (boundaryValues == 0).nonzero()[0]

        # Uniformly randomly sample a cocycle on the sublattice admitted by the
        # chosen edges; reconstruct the labeling on the entire lattice by
        # subbing in the values of c which differ from existing ones.
        return sampleFromKernel(self.lattice.matrices.coboundary, self.lattice.field, zeros)
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        pass
        # self.spins = { face: cocycle[self.lattice.index.faces[face]] for face in self.lattice.faces }
        
        # Dual graph of sublattice of occupied cubes.
        # self.lattice.subgraph = self.lattice.graph.subgraph(
        #     [self.lattice.index.cubes[cube] for cube in self.lattice.cubes if not cube in self.occupied]
        # )