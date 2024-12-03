
import jsonlines as jsl
import gzip
import numpy as np

from .models import Model
from .Chain import Chain


class Recorder:
    """
    Safely "records" states from a Chain by tracking changes between successive
    states, and writing these changes --- alongside the initial state --- to
    file.
    """
    def __init__(self): pass

    def record(self, M: Chain, fp: str, outputs: dict, compressed:bool=True):
        """
        Called to configure the recording apparatus for the Chain, and *should*
        be used as a context manager (i.e. with the `with` statement).

        Args:
            M (Chain): `potts.Chain` object to be recorded; `Recorder` needs
                access to the `Model` and `Lattice` subobjects.
            fp (str): Filename.
            outputs (dict): Dictionary mapping integer keys to initial states.
            compressed (bool): Do we want to compress our outputs? The answer
                should only be "no" if we're debugging something.

        Returns:
            This `Recorder` object.
        """
        # Set up filepath, the actual file I/O object (gzip-compressed if we want
        # compression, otherwise just a standard file), and initialize a JsonLines
        # Writer.
        self._fp = fp
        self._file = gzip.open(f"{self._fp}", "wb") if compressed else open(self._fp, "w")
        self._writer = jsl.Writer(self._file)
        
        # Save the chain for access during iteration; set the "previous" state to
        # be all -1s, since the first state yielded by iterating over the Chain
        # is the initial state, and we need to record the whole thing; fix a list
        # of possible integer states for later.
        self.chain = M
        self.outputs = outputs.copy()
        self.previous = tuple(self.outputs.values())
        
        # Enter the context manager.
        return self.__enter__()


    def store(self, state: tuple) -> None:
        """
        Stores a state yielded by iteration over a `Chain`.

        Args:
            state (tuple): Tuple of assignments or other data points.
        """
        delta = { }

        for k, output in zip(self.outputs.keys(), state):
            # Find the indices where spin assignments *differ*, then categorize the
            # indices based on their spin value. *Should* save some space when
            # writing to file.
            indices = (~np.equal(self.previous[k], output)).astype(int).nonzero()[0]
            faces = { int(s): [] for s in set(np.array(output)) }
            for i in indices: faces[int(output[i])].append(int(i))

            # Record the "delta," i.e. the changes encountered from the previous
            # state.
            delta[k] = faces

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

    def playback(self, S:Model, fp:str, outputs: dict, compressed=True, steps=None):
        """
        Initialize playback; context management.

        Args:
            S (Model): `Model` which matches the configuration on which the
                `Chain` was simulated. **If the `Model` does not have the same
                parameters as the `Model` on which the recorded chain was
                simulated, statistical results from this playback may not be
                accurate.
            fp (str): File in which records are stored.
            outputs (dict): Dictionary mapping integer keys to initial states.
            compressed (bool): Are the data compressed?
            steps (int): How many steps does this chain take?

        Returns:
            This Player object.
        """
        # Configure filepath, file object, and JsonLines Reader.
        self._fp = fp
        self._file = gzip.open(self._fp, "rb") if compressed else open(self._fp, "r")
        self._reader = jsl.Reader(self._file)
        self._steps = steps
        self.outputs = outputs.copy()

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

        for k, output in configuration.items():
            for spin, fs in output.items():
                for f in fs:
                    self.outputs[int(k)][int(f)] = int(spin)

        return tuple(self.outputs.values())
    

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

    
    def progress(self):
        from tqdm.auto import tqdm
        return tqdm(self, total=self._steps)
