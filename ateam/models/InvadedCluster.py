 
import numpy as np
import phat
from itertools import product

from ..arithmetic import essentialCyclesBorn, boundaryMatrix
from ..structures import Lattice
from .Model import Model


class InvadedCluster(Model):
    name = "InvadedCluster"
    
    def __init__(
            self, L: Lattice, homology=1
        ):
        """
        Initializes the plaquette invaded-cluster algorithm.

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
        self.spins = self.initial()
        self.lattice.dimension = homology

        self.coboundary = boundaryMatrix(L.boundary, self.homology, L.field).T

        # Pre-construct the boundary matrix.
        self.phatBoundary = phat.boundary_matrix()
        self.times = set(range(len(self.lattice.flattened)))
        self.indices = np.arange(*self.lattice.tranches[homology])
        self.dimension = max(self.lattice.reindexed.keys())

        # Precompute some other stuff in the interest of speed.
        self.lower = sum([
            list(product([d], np.sort(self.lattice.reindexed[d]).tolist()))
            if d > 0 else [(0, [])]*len(self.lattice.reindexed[d])
            for d in range(homology)
        ], [])

        self.highest = sum([
            list(product([d], np.sort(self.lattice.reindexed[d]).tolist()))
            for d in range(homology+2, self.dimension+1)
        ], [])
    
    
    def initial(self) -> np.array:
        return self.lattice.field.Random(len(self.lattice.boundary[self.homology-1]))
    

    def proposal(self, time):
        """
        Proposal scheme for plaquette invaded-cluster. Each "step" in the chain
        is the set of spins from which a giant cycle first emerges (i.e. the
        "homological percolation" event).

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """

        spins, occupied, satisfied = essentialCyclesBorn(
            self.phatBoundary,
            self.coboundary,
            self.lattice.boundary,
            self.lattice.reindexed,
            self.lattice.tranches,
            self.homology,
            self.lattice.field,
            self.spins,
            self.times,
            self.indices,
            self.lower, self.highest
        )

        return spins, occupied, satisfied
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = cocycle
