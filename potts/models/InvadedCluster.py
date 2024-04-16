 
import numpy as np
from typing import Callable

from ..arithmetic import linalg
from ..structures import Lattice
from ..stats import constant
from .Model import Model


class InvadedCluster(Model):
    name = "InvadedCluster"
    
    def __init__(
            self, L: Lattice, stoppingCondition: Callable=None, initial=None
        ):
        """
        Initializes the invaded-cluster (invasion percolation?) MCMC model on
        the given lattice.

        Args:
            L: The `Lattice` object on which we'll be running experiments.
            initial (np.ndarray): A vector of spin assignments to components.
            stoppingCondition (Callable): A rule which tells the proposal
                algorithm to stopconstructing the current spin and edge configurations.
        """
        self.lattice = L
        self.stoppingCondition = stoppingCondition

        # SW defaults.
        self.state = initial if initial else self.initial()
        self.spins = { face: self.state[self.lattice.index.faces[face]] for face in self.lattice.faces }
        self.occupied = set()


    def initial(self) -> np.array:
        """
        Computes an initial state for the model's Lattice.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        return np.array([np.random.randint(0, self.lattice.field.order) for _ in self.lattice.faces])
    

    def proposal(self, time):
        """
        Proposal scheme for generalized invaded-cluster evolution on the Potts
        model.

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        # First, check whether we're reached the stopping condition; if so, we
        # raise a StopIteration error to force the chain to stop.
        if self.stoppingCondition(self): raise StopIteration

        # Next, randomize the order of the cubes. Then, check whether the cubes
        # in the configuration are sent to 0 by the current cochain; those that
        # are are included unconditionally, and those that aren't are ignored.
        randomized = np.shuffle(self.lattice.cubes)
        includeCubes = [
            cube for cube in randomized
            if self.lattice.field([self.spins[face] for face in cube.faces]).sum() == 0
        ]
        includeCubeIndices = [self.lattice.index.cubes[cube] for cube in includeCubes]
        self.occupied = set(includeCubes)

        # Uniformly randomly sample a cocycle on the sublattice admitted by the
        # chosen edges; reconstruct the labeling on the entire lattice by
        # subbing in the values of c which differ from existing ones.
        return linalg.sampleFromKernel(self.lattice.coboundary, self.lattice.field, includeCubeIndices)
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = { face: cocycle[self.lattice.index.faces[face]] for face in self.lattice.faces }