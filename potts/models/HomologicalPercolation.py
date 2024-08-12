 
import numpy as np
import phat

from ..arithmetic import essentialCyclesBorn
from ..structures import Lattice
from .Model import Model


class HomologicalPercolation(Model):
    name = "HomologicalPercolation"
    
    def __init__(
            self, L: Lattice, homology=1, mesh=1000, temperatureFunction=None,
            initial=None
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

        # Partition the cells into "lower," "target," "higher," and "other.""
        # We'll only ever be modifying the middle two components of the partition.
        self.tranches = self.lattice.tranches
        self.boundary = sum([list(L.boundary[d]) for d in range(0, len(L.corners)+1)], [])

        self.lowerCells = np.array([
            np.sort(self.boundary[j])
            for j in range(0, self.tranches[homology-1])
        ])

        self.targetCells = np.array([
            np.sort(self.boundary[j])
            for j in range(self.tranches[homology-1], self.tranches[homology])
        ])

        print(len(self.boundary))

        self.higherCells = np.array([
            self.boundary[j]
            for j in range(self.tranches[homology], self.tranches[homology+1])
        ])

        self.higherCellsFlat = self.higherCells.flatten()

        self.otherCells = np.array([
            np.sort(self.boundary[j])
            for j in range(self.tranches[homology+1], self.tranches[-1])
        ]) if homology+1 < len(L.corners)-1 else np.array([])

        self.dimensions = list(range(len(L.corners)+1))

        
        # Static sets of indices to pass to the shuffler.
        self.times = set(range(self.tranches[-1]))
        self.targetRelativeIndices = np.array(range(0, len(self.lattice.skeleta[homology])))
        self.targetAbsoluteIndices = np.array(range(self.tranches[homology-1], self.tranches[homology]))
        self.lowestRelativeTargetIndex = min(self.targetAbsoluteIndices)
        self.highestRelativeTargetIndex = max(self.targetAbsoluteIndices)

        # Pre-construct the boundary matrix.
        self.boundary = phat.boundary_matrix()


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
        return essentialCyclesBorn(
            self.boundary,
            self.targetAbsoluteIndices,
            self.targetRelativeIndices,
            self.lowestRelativeTargetIndex,
            self.highestRelativeTargetIndex,
            self.lowerCells,
            self.targetCells,
            self.higherCellsFlat,
            self.otherCells,
            (len(self.higherCells), -1),
            self.dimensions,
            self.times
        )
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = cocycle
