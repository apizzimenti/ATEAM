
import jsonlines as jsl
import gzip
import numpy as np

from . import Model
from .Chain import Chain
from typing import Self


class Recorder:
    """
    Safely "records" states from a Chain by tracking changes between successive
    states, and writing these changes --- alongside the initial state --- to
    file.
    """
    def __init__(self): pass

    def record(self, M: Chain, fp: str, compressed:bool=True) -> Self:
        """
        Called to configure the recording apparatus for the Chain, and *should*
        be used as a context manager (i.e. with the `with` statement).

        Args:
            M (Chain): `potts.Chain` object to be recorded; `Recorder` needs
                access to the `Model` and `Lattice` subobjects.
            fp (str): Filename.
            compressed (bool): Do we want to compress our outputs? The answer
                should only be "no" if we're debugging something.

        Returns:
            This `Recorder` object.
        """
        # Set up filepath, the actual file I/O object (gzip-compressed if we want
        # compression, otherwise just a standard file), and initialize a JsonLines
        # Writer.
        self._fp = fp
        self._file = gzip.open(f"{self._fp}.gz", "wb") if compressed else open(self._fp, "w")
        self._writer = jsl.Writer(self._file)
        
        # Save the chain for access during iteration; set the "previous" state to
        # be all -1s, since the first state yielded by iterating over the Chain
        # is the initial state, and we need to record the whole thing; fix a list
        # of possible integer states for later.
        self.chain = M
        self.previous = [-1]*len(self.chain.model.lattice.faces)
        self._states = list(range(0, M.model.lattice.field.order))
        
        # Enter the context manager.
        return self.__enter__()


    def store(self, state: np.array) -> None:
        """
        Stores a state yielded by iteration over a `Chain`.

        Args:
            state (np.array): List of spin assignments for the faces of the
                underlying lattice.
        """
        # Find the indices where spin assignments *differ*, then categorize the
        # indices based on their spin value. *Should* save some space when
        # writing to file.
        indices = (~np.equal(self.previous, state)).astype(int).nonzero()[0]
        faces = { s: [] for s in self._states }
        for i in indices: faces[int(state[i])].append(int(i))

        # Record the "delta," i.e. the changes encountered from the previous
        # state.
        delta = {
            "faces": faces,
            "cubes": [self.chain.model.lattice.index.cubes[i] for i in self.chain.model.occupied]
        }

        # Write the delta to file and update the previous state to be the current
        # one.
        self._writer.write(delta)
        self.previous = state

    def __enter__(self):
        """
        Required context management magic method.
        """
        return self


    def __exit__(self, exc_type, exc_value, exc_tb):
        """
        Required context management magic method; kills the writer and the file.
        """
        self._writer.close()
        self._file.close()


class Player():
    """
    Safely "replays" recorded `Chain` output.
    """
    def __init__(self): pass

    def playback(self, S:Model, fp:str, compressed=True) -> Self:
        """
        Initialize playback; context management.

        Args:
            S (Model): `Model` which matches the configuration on which the
                `Chain` was simulated. **If the `Model` does not have the same
                parameters as the `Model` on which the recorded chain was
                simulated, statistical results from this playback may not be
                accurate.
            fp (str): File in which records are stored.
            compressed (bool): Are the data compressed?

        Returns:
            This Player object.
        """
        # Configure filepath, file object, and JsonLines Reader.
        self._fp = fp
        self._file = gzip.open(self._fp, "rb") if compressed else open(self._fp, "r")
        self._reader = jsl.Reader(self._file)
        
        # Save `Model` for later.
        self.model = S

        # Enter context management.
        return self.__enter__()
    

    def __iter__(self):
        """
        This `Player` is a generator which yields states of the recorded chain.
        """
        return self


    def __next__(self):
        """
        Iteration magic method.
        """
        # Get the next configuration by calling __next__() on the reader.
        try: configuration = self._reader.read()
        except: raise StopIteration

        # Get cubes and faces.
        faces = configuration["faces"]
        cubes = configuration["cubes"]

        # Update the faces' spins and the occupied edges.
        update = {}
        for spin, fs in faces.items():
            for f in fs:
                update[self.model.lattice.faces[f]] = spin

        self.model.spins.update(update)
        self.model.occupied = set(self.model.lattice.cubes[c] for c in cubes)

        # Return the current state.
        return np.array([self.model.spins[f] for f in self.model.lattice.faces])


    def __enter__(self):
        """
        Required context management magic method.
        """
        return self


    def __exit__(self, exc_type, exc_value, exc_tb):
        """
        Required context management magic method; kills the reader and the file.
        """
        self._reader.close()
        self._file.close()