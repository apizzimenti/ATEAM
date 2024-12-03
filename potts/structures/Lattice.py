
import galois
import numpy as np
import json
from pathlib import Path
from ast import literal_eval as le

from ..arithmetic import cubicalComplex, boundaryMatrix, flatten


class Matrices:
    boundary = None
    coboundary = None


class Lattice:
    """
    A class specifying a cubical complex.
    """
    def __init__(self): pass
    
    def fromCorners(self, corners, field=2, dimension=np.inf, periodicBoundaryConditions=True):
        """
        Creates a cell complex with the given corners made of cells of the
        provided dimension.

        Args:
            corners (list): Corners of the lattice; determines the maximal
                cell dimension.
            field (int): Characteristic of finite field from which cells take
                coefficients.
            dimension (int): Maximal cell dimension; if this argument is larger
                than that permitted by the underlying cell structure, this is
                re-set to the maximal legal dimension. Defaults to 1, so the
                lattice is constructed only of vertices (0-cells) and edges
                (1-cells).
            periodicBoundaryConditions (bool): Do we use periodic boundary
                conditions (i.e. are we making a torus)?
        """
        self._construct(dimension, corners, periodicBoundaryConditions, field)
        return self
    

    def _construct(
            self, dimension, corners, periodicBoundaryConditions, field, data=None
        ):
        self.dimension = min(dimension, len(corners))
        self.periodicBoundaryConditions = periodicBoundaryConditions
        self.corners = corners
        self.field = galois.GF(field)
        self.matrices = Matrices()

        # Specify sparse "boundary matrices," which we'll convert to *actual*
        # matrices for the purposes of our stuff.
        if not data:
            v, b, B = cubicalComplex(self.corners, self.dimension, self.field, self.periodicBoundaryConditions)
            self.vertexMap, self.boundary = v, b
        else:
            self.vertexMap = data.get("vertexMap", None)
            self.boundary = {
                int(t): data["boundary"][t].astype(int) for t in data["boundary"].files
            }
            data["boundary"].close()
            B = boundaryMatrix(self.boundary, self.dimension, self.field)

        # Construct the finite field and boundary matrices.
        self.matrices.boundary = B
        self.matrices.coboundary = B.T
        self.reindexer, self.reindexed, self.flattened = flatten(self.boundary, self.dimension)

        # Get index ranges.
        self.tranches = { 0: (0, len(self.boundary[0])) }

        for d in range(1, self.dimension+1):
            self.tranches[d] = (self.tranches[d-1][1], self.tranches[d-1][1] + len(self.boundary[d]))



    def toFile(self, fp:str):
        """
        JSON-serializes this object and writes it to file so we can reconstruct
        it later.

        Args:
            fp (str): Filepath.
        """
        # Write compressed boundary matrix and vertex maps to file.
        absolute = Path(fp).resolve()
        root = absolute.parent
        stem = absolute.stem
        boundaryFile = root/f".{stem}.lattice.npz"
        np.savez_compressed(boundaryFile, **{str(t): v for t, v in self.boundary.items()})

        with open(fp, "w") as write:
            json.dump(
                {
                    "field": self.field.order,
                    "dimension": self.dimension,
                    "periodicBoundaryConditions": int(self.periodicBoundaryConditions),
                    "corners": self.corners,
                    "boundary": str(boundaryFile),
                    "vertexMap": { str(k): int(v) for k, v in self.vertexMap.items() }
                }, write
            )
        

    def fromFile(self, fp:str, vertexMap=False):
        """
        Reconstructs a serialized Lattice.

        Args:
            fp (str): Filepath.
        """
        with open(fp, "r") as read:
            # Read file into memory.
            serialized = json.load(read)

            # Set field and dimension.
            field = serialized["field"]
            dimension = int(serialized["dimension"])
            periodicBoundaryConditions = bool(serialized["periodicBoundaryConditions"])
            corners = serialized["corners"]
            data = { "boundary": np.load(serialized["boundary"], allow_pickle=True) }
            if vertexMap: data["vertexMap"] = { le(k): int(v) for k, v in serialized["vertexMap"].items() }

            # Construct indices, boundary matrix, and graph.
            self._construct(
                dimension, corners, periodicBoundaryConditions, field, data
            )
