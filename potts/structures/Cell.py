
from itertools import combinations


class Cell:
    def __init__(self, dimension, spin=1, corner=None, vertices=None):
        """
        Initializes a Cell object.

        Args:
            faces (list): Ordered collection of faces.
            spin (Number): Spin (coefficient) of this Cell; typically an integer.
            corner (tuple): Coordinate of the bottom-left corner of this Cell.
        """
        # If there are no faces, we're just a point; otherwise, we're of one
        # dimension *higher* than the first face.
        self.dimension = dimension
        self.vertices = vertices
        self.corner = corner
        self.spin = spin

        # Construct the remaining faces if we're of the right dimension.
        if corner: self._constructMaxDimensionalFaces()
        else: self._constructLowerDimensionalFaces()
        # self.allComponents = self._bfs(self, set())
        # print(self.allComponents)


    @staticmethod
    def _bfs(cell, components):
        if cell not in components: components.add(cell)
        for face in cell.faces: Cell._bfs(face, components)

        return components
    

    @staticmethod
    def _increment(t, at, by=1):
        t = list(t)
        for a in at: t[a] += by
        return tuple(t)
    

    def _constructLowerDimensionalFaces(self):
        """
        Construct the faces of this Cell; if this method is called, we know that
        this cell is *not* a cell of maximal dimension.
        """
        # Handle vertices and edges separately.
        if self.dimension == 0:
            self.corner = self.vertices[0]
            self.faces = None
            return
        
        if self.dimension == 1:
            self.corner = self.vertices[0]
            self.faces = [Cell(dimension=0, vertices=[v]) for v in self.vertices]
            return

        # Set the anchoring corner and create an empty bucket for faces.
        self.corner = self.vertices[0]
        faces = []

        # As before, a "face" of a set of coordinates is found by fixing the
        # value of a single coordinate entry, then finding all other coordinates
        # which have the same fixed value; this is equivalent to finding the
        # (k-1)-dimensional subspaces of a k-dimensional space.
        for fixed in range(len(self.corner)):
            for level in [0, 1]:
                # Create a bucket for the vertices of the face, and set a basepoint
                # (i.e. the "corner," but incrementing the fixed coordinate entry
                # by 1).
                face = []
                basepoint = self._increment(self.corner, [fixed], level)

                # Add all vertices matching on the fixed entry.
                for vertex in self.vertices:
                    if vertex[fixed] == basepoint[fixed]: face.append(vertex)
                
                # We capture extra information when determining faces (i.e. we get
                # the entire boundary or nothing for two cases), so we remove
                # extraneous data.
                if len(face) != 2*(self.dimension-1): continue
                faces.append(face)

        # Construct a *set* of all the vertices, so we can distinguish this Cell
        # from other Cells of the same dimension.
        vertexSets = [set(face) for face in faces]
        vertexSet = set().union(*vertexSets)
        self.vertices = tuple(sorted(list(vertexSet)))

        # For each of the coordinate sets of each face, create a new Cell.
        self.faces = [
            Cell(dimension=self.dimension-1, vertices=c)
            for c in faces
        ]
    

    def _constructMaxDimensionalFaces(self):
        """
        Construct the faces of this Cell; if this method is called, we know that
        this cell is a cell of maximal dimension.
        """
        # Create an empty bucket for faces.
        faces = []

        # Hold one of the coordinates constant and increment each of the others.
        for fixed in range(len(self.corner)):
            # These are the coordinates we're allowed to increment.
            allowed = list(set(range(len(self.corner)))-{fixed})

            # Find all the ways we can change each corner; this cycles us around
            # each face.
            incrementable = [
                [list(c) for c in combinations(allowed, r=k)]
                for k in range(self.dimension)
            ]
            
            # For each "level" (i.e. each value change for the coordinate held
            # fixed) and each of the index sets, increment the appropriate
            # coordinates, generating the vertices composing each of the faces.
            for level in [0,1]:
                # Bucket for a new face.
                face = []

                # For each index set, adjust the appropriate coordinates,
                # generating the vertices bounding a face; add this new face to
                # the list of faces.
                for indexSet in incrementable:
                    for indices in indexSet:
                        basepoint = self._increment(self.corner, [fixed], by=level)
                        face.append(self._increment(basepoint, indices))

                faces.append(face)

        # Construct a *set* of all the vertices, so we can distinguish this Cell
        # from other Cells of the same dimension.
        vertexSets = [set(face) for face in faces]
        vertexSet = set().union(*vertexSets)
        self.vertices = tuple(sorted(list(vertexSet)))

        # For each of the coordinate sets of each face, create a new Cell.
        self.faces = [
            Cell(dimension=self.dimension-1, vertices=c)
            for c in faces
        ]


    def __hash__(self):
        """
        Cells are uniquely determined by their (ordered!) vertices and dimension.

        Returns:
            Hash key for this Cell.
        """
        return hash(self.vertices + (self.dimension,))


    def __str__(self):
        """
        Returns:
            String representation of this Cell.
        """
        return f"dim {self.dimension}\t{self.vertices}"

    
    def __eq__(self, other):
        """
        Determines whether this Cell and `other` are the same.
        """
        return hash(other) == hash(self)
