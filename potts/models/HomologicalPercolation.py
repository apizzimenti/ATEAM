 
import numpy as np
import phat

from ..arithmetic import essentialCyclesBorn
from ..structures import Lattice
from .Model import Model


class HomologicalPercolation(Model):
    name = "HomologicalPercolation"
    
    def __init__(
            self, L: Lattice, homology=1
        ):
        """
        Initializes homological percolation.

        Args:
            L (Lattice): The `Lattice` object on which we'll be running experiments.
            homology (int=1): Computing the `homology`th homology group of the
                complex.
            mesh (int=1000): Number of temperature sample points.
            initial (None): Dummy.
            temperatureFunction (None): Dummy.
        """
        self.lattice = L
        self.homology = homology

        # Change the Lattice's dimension and construct an initial spin configuration.
        # We have to change the lattice's dimension and reconstruct the boundary
        # matrix, though.
        self.spins = self._initialSpins()
        self.lattice.dimension = homology
        self.lattice._constructBoundaryMatrix()

        # Pre-construct the boundary matrix.
        self.phatBoundary = phat.boundary_matrix()
        self.times = set(range(self.lattice.tranches[-1]))
        self.indices = self.lattice.skeleta[homology].copy()
    
    
    def _initialSpins(self) -> np.array:
        return self.lattice.field([0]*len(self.lattice.skeleta[self.homology-1]))


    def initial(self) -> np.array:
        """
        Computes an initial state for the model's Lattice.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        return np.array([np.random.randint(0, self.lattice.field.order) for _ in self.lattice.skeleta[self.homology]])
        

    def proposal(self, time):
        """
        Proposal scheme for homological percolation. Each "step" in the chain
        is the set of spins from which a giant cycle first emerges (i.e. the
        "homological percolation" event).

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """

        spins, occupied, satisfied = essentialCyclesBorn(
            self.phatBoundary,
            self.lattice.matrices.coboundary,
            self.lattice.boundary,
            self.lattice.tranches,
            self.lattice.skeleta,
            self.homology,
            self.lattice.field,
            self.spins,
            self.times,
            self.indices
        )

        return spins, occupied, satisfied
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = cocycle
