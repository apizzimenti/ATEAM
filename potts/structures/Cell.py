
from functools import reduce

class ReducedCell():
    def __init__(self, faces, vertex=False):
        self.dimension = 0 if vertex else faces[0].dimension+1
        self.vertices = [faces] if vertex else list(sorted(reduce(lambda A, B: A+B, [face.vertices for face in faces])))
        self.encoding = tuple(self.vertices) + (self.dimension,)
        # self.faces = set(faces)
    
    def __eq__(self, other): return self.encoding == other.encoding
    def __str__(self): return str(self.encoding)
    def __hash__(self): return hash(self.encoding)


class Vertex:
    # TODO consolidate the GraphLattice and Lattice classes!
    pass
