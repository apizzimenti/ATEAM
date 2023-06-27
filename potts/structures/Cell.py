
import numpy as np


class Cell:
    # Default spin.
    spin = 0;

    def __init__(self, coordinates, orientation=1, index=0):
        """
        Initializes a cell of any dimension.

        Args:
            coordinates (list): Will be integers (0-cell) or cell entries
                (k-cell, k > 0).
            orientation (int): Which way do we traverse the entries in `coordinates`?
            index (int): When constructing matrices, the index to which this cell
                corresponds.
        """
        # This determines which kind of cell we're dealing with: if it's a
        # 0-dimensional cell, then `coordinates` will have integer entries; if
        # it's any higher-dimensional one, its entries will be instances of cell.
        if type(coordinates[0]) in {int, np.int32, np.int64}:
            self.dimension = 0
        else:
            self.dimension = coordinates[0].dimension + 1

        self.coordinates = coordinates
        self.orientation = orientation
        self.index = index

    def __repr__(self): return self.__str__()

    def __str__(self):
        """
        Stringification.
        """
        if self.dimension < 1:
            return str(tuple(self.coordinates))
        else:
            return f"{self.dimension}-dimensional cell at index {self.index}"
        
    def __eq__(self, other):
        """
        Testing for equality; this asks whether the coordinates of each cell
        are the same.
        """
        return self.coordinates == other.coordinates
    
    def __lt__(self, other):
        pass

    def __hash__(self):
        """
        Returns a unique property from each cell.
        """
        return hash(tuple(self.coordinates))
