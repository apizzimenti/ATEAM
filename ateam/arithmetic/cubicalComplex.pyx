
import numpy as np
from itertools import combinations as combs, product
from math import comb


def flatten(Complex, D):
	"""
	Flattens the given Complex for a PHAT-style specification.
	"""
	# Line up indices.
	reindexer = {
		d: np.arange(len(Complex[d])) + sum([len(Complex[t]) for t in range(d)])
		if d > 0 else np.arange(len(Complex[0]))
		for d in range(D+1)
	}

	# Reindex.
	reindexed = {
		d: reindexer[d-1][Complex[d]]
		for d in range(1, D+1)
	}

	reindexed[0] = Complex[0]

	# Flatten.
	flattened = sum([reindexed[t].tolist() for t in range(D+1)], [])

	return reindexer, reindexed, flattened


def boundaryMatrix(complex, d, field):
	"""
	Construct the boundary (operator) matrix given the boundary specification
	"""
	faces = complex[d-1]
	cubes = complex[d]

	# Construct the empty matrix; this will *only* store 0, 1, and -1 (i.e. whatever
	# the additive inverse of 1 is in the field; to be safe, we use the below method)
	# rather than straight subtraction).
	m, n = len(faces), len(cubes)
	B = field(np.zeros((m,n), dtype=int))

	for j in range(len(cubes)):
		indices = cubes[j]
		coefficients = np.resize([-field(1), 1], len(indices))
		B[:,j][indices] = coefficients
	
	return B


def hammingCube(D):
	"""
	Constructs a Hamming cube of dimension D.

	Args:
		D (int): Top dimension.

	Returns:
		A dictionary of NumPy arrays corresponding to (compressed) boundary
		matrices at each dimension.
	"""
	# Bucket for the boundary.
	Complex = {}

	# Specify the vertices and edges, building everything else from there.
	basis = np.array([2**k for k in range(0, D)])
	vertices = np.array(range(0, 2**D))
	Complex[0] = vertices
	limit = max(vertices)

	edges = np.array([
		np.array([v, v+b])
		for v, b in product(vertices, basis)
		if v != b and v+b <= limit and v&(v+b) == v
	])
	Complex[1] = edges

	for d in range(2, D+1):
		# Construct an appropriately-sized bucket for the boundary of the cube.
		cubes = 2**(D-d)*comb(D, d)
		facesPerCube = 2*d
		_boundary = np.zeros((cubes, facesPerCube), dtype=int)

		# Get all the combinations of things in the previous skeleton and
		# filter them; we effectively check that the composition-of-boundary-maps
		# condition holds by requiring every element of the (d-2) skeleton to
		# appear *exactly twice* in whatever we're constructing.
		combinations = np.array(list(combs(range(len(Complex[d-1])), facesPerCube)))\
			.flatten()\
			.reshape((-1, facesPerCube))
		
		filled = 0
		
		for combination in combinations:
			faces = Complex[d-1][combination]
			indices = faces.flatten()
			
			if len(np.unique(indices)) == len(indices)/2:
				_boundary[filled] = combination
				filled += 1

		Complex[d] = _boundary

	return Complex


def cubicalComplex(corners, D, field, periodic=True):
	"""
	Constructs a cubical complex completely specified by a boundary matrix and
	a vertex map.

	Args:
		corners (iterable): An iterable specifying the extent of the lattice. If
			`periodic` is truthy, antipodal vertices are identified; for example,
			if `corners = [3,3]`, then the vertices at indices (0,0), (0,3), (3,0),
			and (3,3) are identified.
		periodic (bool): If truthy, we construct a \(d\)-fold torus, where \(d\)
			is the topmost dimension of the lattice (as determined by the number
			of corners specified). For example, if `periodic` is falsy, then
			`corners = [3,3]` hands back a 4x4 grid (zero-indexed).
	"""
	# Get the boundary specification for a D-dimensional Hamming cube; we use
	# this cube to enforce an ordering on the vertices (and thus edges, etc.) so
	# the boundary/coboundary operations make sense.
	boundary = hammingCube(D)

	# Specify the "basis" Hamming cube and the remaining separately; this
	# helps when we need to shift the cubes around and properly index.
	_basis = list(product(*[range(2) for c in corners]))
	basis = np.array(_basis)

	i = 1
	_extra = list(product(*[range(c+i) for c in corners]))
	extra = np.array([v for v in _extra if v not in _basis])

	# Construct a complete list of coordinates and indices; create a *dict* to
	# make creation a bit easier.
	vertices = np.concatenate([basis, extra])
	indices = np.arange(len(vertices))
	tups = [tuple(v) for v in vertices]

	verticesToIndices = { tups[t]: t for t in indices }

	if periodic:
		# Identify antipodal vertices.
		for v in vertices: verticesToIndices[tuple(v)] = verticesToIndices[tuple(v%corners)]

		# After identification, count the unique values we see, then re-map.
		unique = np.unique(list(verticesToIndices.values()))
		vertices = np.arange(len(unique))
		indices[unique] = vertices
		verticesToIndices = { t: indices[verticesToIndices[t]] for t in tups }
		
	# Explode and do some reshaping. We need everything in terms of their
	# local coordinates, so we'll have to do some re-indexing on the higher-
	# dimensional faces.        
	Complex = {}
	Complex[0] = vertices
	Skeleta = {}

	anchors = np.array(list(product(*[range(c) for c in corners])))
	anchorMapping = { tuple(a): { d: [] for d in range(D+1) } for a in anchors }

	# Handle edges first, then reindex.
	for anchor in anchors:
		stacked = np.array([anchor]*len(basis))
		localIndices = stacked+basis
		localFaces = np.array([verticesToIndices[tuple(v)] for v in localIndices])

		anchorMapping[tuple(anchor)][0] = localFaces
		anchorMapping[tuple(anchor)][1] = localFaces[boundary[1]]

	allEdges = np.concatenate([a[1] for a in anchorMapping.values()])
	uniqueEdges = np.unique(allEdges, axis=0)
	Complex[1] = uniqueEdges
	Skeleta[1] = { tuple(e): j for j, e in enumerate(uniqueEdges) }

	# Now do everything else! We shouldn't get overlaps here.
	for d in range(2, D+1):
		for anchor in anchors:
			localFaces = anchorMapping[tuple(anchor)][d-1]
			localBoundary = np.zeros(boundary[d].shape).astype(int)

			for j in range(len(boundary[d])):
				localized = np.array([Skeleta[d-1][tuple(localFaces[f])] for f in boundary[d][j]])
				localBoundary[j] = localized

			anchorMapping[tuple(anchor)][d] = localBoundary
		
		allBoundaries = np.concatenate([a[d] for a in anchorMapping.values()])
		uniqueBoundaries = np.unique(allBoundaries, axis=0)
		Complex[d] = uniqueBoundaries
		Skeleta[d] = { tuple(f): j for j, f in enumerate(uniqueBoundaries) }

		# Delete repeated keys in the vertex mapping.
		vti = {}
		seen = set()

		for v, k in verticesToIndices.items():
			if k not in seen: vti[v] = k
			seen.add(k)

		verticesToIndices = vti

	return verticesToIndices, Complex, boundaryMatrix(Complex, D, field)
