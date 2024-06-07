 
import numpy as np
from typing import Callable

from ..arithmetic import linalg
from ..structures import Lattice
from ..stats import constant
from .Model import Model


class HomologicalPercolation(Model):
    name = "HomologicalPercolation"
    
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
        Proposal scheme for generalized Swendsen-Wang evolution on the Potts model.

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        self.temperature = self.temperatureFunction(time)
        p = 1-np.exp(self.temperature)
        assert 0 <= p <= 1

        # Assign each site "open" or "closed".
        U = np.random.uniform(size=len(self.lattice.faces))
        return (U < p).astype(int)
    

    def assign(self, cocycle: np.array):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = cocycle