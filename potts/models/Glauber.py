 
import numpy as np
from typing import Callable

from ..structures import Lattice
from ..stats import constant
from .Model import Model


class Glauber(Model):
    name = "Glauber"
    
    def __init__(
            self, L: Lattice, temperatureFunction: Callable=constant(-0.6),
            initial=None
        ):
        """
        Initializes Glauber dynamics on the Potts model.

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
        self.state = self.lattice.field(self.state)
        self._shift = self.lattice.field(np.zeros(len(self.state), dtype=int))
        self._shiftIndex = 0
        self._indices = list(range(len(self.state)))

        self.spins = { face: self.state[self.lattice.index.faces[face]] for face in self.lattice.faces }
        self.occupied = set()


    def initial(self) -> np.array:
        """
        Computes an initial state for the model's Lattice.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        return self.lattice.field.Random(len(self.lattice.faces))
    

    def proposal(self, time):
        """
        Proposal scheme for generalized Glauber dynamics on the Potts model:
        uniformly randomly chooses a face in the complex, flips the face's spin,
        and returns the corresponding cocycle.

        Args:
            time (int): Step in the chain.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        # Choose a location to "flip," then revisit the last shifted vertex.
        loc = np.random.randint(len(self.lattice.faces))
        self._shift[self._shiftIndex] = 0
        self._shift[loc] = self.lattice.field.Random()
        self._shiftIndex = loc

        return self.state+self._shift

    def assign(self, cocycle):
        """
        Updates mappings from faces to spins and cubes to occupations.

        Args: 
            cocycle (np.array): Cocycle on the sublattice.
        """
        self.spins = { face: cocycle[self.lattice.index.faces[face]] for face in self.lattice.faces }
        self.state = cocycle
        
        # Dual graph of sublattice of occupied cubes.
        # self.lattice.subgraph = self.lattice.graph.subgraph(
        #     [self.lattice.index.cubes[cube] for cube in self.lattice.cubes if not cube in self.occupied]
        # )