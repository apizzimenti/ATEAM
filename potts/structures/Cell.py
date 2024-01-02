
from functools import reduce
from ..utils import binaryEncode, binaryUnencode


class ReducedCell():
    def __init__(self, faces, vertex=False):
        self.dimension = 0 if vertex else faces[0].dimension+1
        self.vertices = [faces] if vertex else list(sorted(reduce(lambda A, B: A|B, [set(face.vertices) for face in faces])))
        self.encoding = tuple(self.vertices) + (self.dimension,)
        self.faces = set() if self.dimension == 0 else set(faces)
    
    def __eq__(self, other): return self.encoding == other.encoding
    def __str__(self): return str(self.encoding)
    def __hash__(self): return hash(self.encoding)


class BinaryReducedCell():
    def __init__(self, faces, vertex=False):
        self.dimension = 0 if vertex else faces[0].dimension+1
        self.vertices = [faces] if vertex else list(sorted(reduce(lambda A, B: A|B, [set(face.vertices) for face in faces])))
        self.encoding = tuple(self.vertices) + (self.dimension,)
        self.faces = set() if self.dimension == 0 else set(faces)

    def shift(self, k):
        self.vertices = [v+k for v in self.vertices]
        self.encoding = tuple(self.vertices) + (self.dimension,)
    
    def __eq__(self, other): return self.encoding == other.encoding
    def __str__(self): return str(self.encoding)
    def __hash__(self): return hash(self.encoding)


class Vertex:
    # TODO consolidate the GraphLattice and Lattice classes!
    pass
