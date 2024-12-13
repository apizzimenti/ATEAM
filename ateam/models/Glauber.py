 
import numpy as np
from typing import Callable

from ..arithmetic import evaluateCochain
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
        self.faces = self.lattice.boundary[self.lattice.dimension-1]
        self.plaquettes = self.lattice.boundary[self.lattice.dimension]
        self._shift = self.lattice.field(np.zeros(len(self.faces), dtype=int))
        self._shiftIndex = 0
        self._indices = list(range(len(self.faces)))

        self.spins = initial if initial is not None else self.initial()
        coboundary = evaluateCochain(self.plaquettes, self.spins)
        self.closed = -(coboundary > 0).nonzero()[0].sum()
        self.satisfied = (coboundary == 0).nonzero()[0]


    def initial(self):
        """
        Computes an initial state for the model's Lattice.

        Returns:
            A NumPy array representing a vector of spin assignments.
        """
        return self.lattice.field.Random(len(self.faces))
    

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
        loc = np.random.randint(len(self.faces))
        self._shift[self._shiftIndex] = 0
        self._shift[loc] = self.lattice.field.Random()
        self._shiftIndex = loc

        spins = self.spins + self._shift
        coboundary = evaluateCochain(self.plaquettes, spins)
        closed = -(coboundary > 0).nonzero()[0].sum()
        satisfied = (coboundary == 0).nonzero()[0]
        energy = np.exp(self.temperatureFunction(time)*(self.closed - closed))

        if np.random.uniform() < min(1, energy):
            self.closed = closed
            self.satisfied = satisfied
            self.spins = spins

        satisfiedConfiguration = self.lattice.field.Zeros(len(self.plaquettes))
        satisfiedConfiguration[self.satisfied] = 1
        
        return self.spins, satisfiedConfiguration


    def assign(self, cochain):
        """
        Updates mappings from faces to spins and cubes to occupations.
        
        Args:
            cochain (galois.FieldArray): Cochain on the lattice.
        
        Returns:
            None.
        """
        self.spins = cochain
    