
import numpy as np
import matplotlib.pyplot as plt
import rustworkx as rx
from functools import reduce
from matplotlib.patches import Rectangle
from itertools import combinations, product


def shortestPath(L, assignment):
    # Construct a graph.
    G = rx.PyGraph(multigraph=False)

    # Add vertices and edges.
    G.add_nodes_from(range(len(L.boundary[0])))
    G.add_edges_from_no_data([
        tuple(e) for i, e in enumerate(L.boundary[1])
        if assignment[i]
    ])

    # Look through all the cycles to see whether they pass through "antipodal"
    # points.
    w, h = L.corners

    antipodes = {
        0: [w-1, w*(h-1)],
        w-1: [0, (w*h)-1],
        w*(h-1): [0, (w*h)-1],
        (w*h)-1: [w-1, w*(h-1)]
    }

    bottom = set(range(1, w-1))
    left = set(range(w, w*(h-1), w))

    for k in bottom: antipodes.update({ k: [k+w*(h-1)] })
    for k in left: antipodes.update({ k: [k+w-1] })

    antipodes = { k: set(v) for k, v in antipodes.items() }
    boundary = set(antipodes.keys())
    interior = set(range(len(G.nodes())))-boundary

    B = rx.cycle_basis(G)
    B = [c for c in B if len(c) >= min(L.corners)]

    # This *kinda* gets the cycles, but not really.
    B = [
        c for c in B
        if set(c).issubset(boundary) or (any(f in boundary and set(c)&antipodes[f] for f in c) and len(set(c)&interior) >= min(L.corners))
    ]

    # Grab the cycle of shortest length (by default) and get edge tuples; then,
    # get the coordinates of the edges' endpoints, and return.
    if B:
        essential = list(sorted(B, key=lambda L: len(L)))[0]
        edgeVertices = list(set(essential))
        edges = list(zip([essential[-1]] + essential[:-1], essential))
        edges = np.array(edges)
        edges.sort(axis=1)

        edgeIncluded = np.array([
            int(list(r) in edges) for r in L.boundary[1]
        ])

        return edgeIncluded, edgeVertices
    return [], []


def points(L, dimension, scale=1, adj=2/3):
    _vertices = [k for k, v in L.vertexMap.items()]
    vertices = [(x*scale, y*scale) for x, y in _vertices]
    edges = list((vertices[u], vertices[v]) for (u,v) in L.boundary[dimension])
    external = [[] for _ in edges]
    adj = adj*scale

    # Move through the list of edges, re-casting edges that are on the boundary
    # of the torus.
    for j, edge in enumerate(edges):
        ((ux, uy), (vx, vy)) = edge
        if ux == vx:
            if abs(uy-vy) > 1:
                external[j] += [
                    ((vx, vy), (vx, vy+adj)),
                    ((ux, uy), (ux, uy-adj))
                ]
        if uy == vy:
            if abs(ux-vx) > 1:
                external[j] += [
                    ((vx, vy), (vx+adj, vy)),
                    ((ux, uy), (ux-adj, uy))
                ]
    
    # Reform the sets of vertices and edges so they're numpy arrays.
    vertices = np.array(vertices).T
    return vertices, (edges, external)


def lattice2D(
        L, assignment,
        padding=0.1,
        shortest=False,
        vertexOccupiedArgs=dict(),
        vertexShortestArgs=dict(),
        edgeOccupiedArgs=dict(),
        edgeShortestArgs=dict(),
        edgeVacantArgs=dict(),
        squareArgs=dict(),
    ):
    """
    Plots the flat torus specified by `L`.

    Args:
        L (Lattice): `potts.structures.Lattice` object.
        assignment (iterable): An iterable of 0s and 1s specifying which edges
            to draw.
        padding (float): Space between the maximum coordinate pair on a given
            axis and the edge of the figure.
        vertexArgs (dict): Arguments passed to `plt.scatter()`.
        edgeOccupiedColor (str): Color of occupied edges.
        edgeOccupiedWidth (float): Width of occupied edges.
        edgeOccupiedColor (str): Color of vacant edges.
        edgeOccupiedWidth (float): Width of vacant edges.
        edgeArgs (dict): Arguments passed to `plt.plot()`.
        squareArgs (dict): Arguments passed to `patches.Rectangle()`.

    Returns:
        `(matplotlib.Figure, matplotlib.Axes)` subplots pair.
    """
    # Defaults.
    _vertexOccupiedArgs=dict(s=6, color="#C19A6B", zorder=10, facecolor="none", lw=1)
    _vertexShortestArgs=dict(s=6, color="#734F96", zorder=10, facecolor="none", lw=1)
    _edgeOccupiedArgs=dict(lw=1.5, color="#C19A6B", zorder=1)
    _edgeShortestArgs=dict(lw=1.6, color="#734F96", zorder=1)
    _edgeVacantArgs=dict(lw=0, color="none")
    _squareArgs=dict(facecolor="#87A96B75", edgecolor="none", zorder=0)

    _vertexOccupiedArgs.update(vertexOccupiedArgs)
    _vertexShortestArgs.update(vertexShortestArgs)
    _edgeOccupiedArgs.update(edgeOccupiedArgs)
    _edgeShortestArgs.update(edgeShortestArgs)
    _edgeVacantArgs.update(edgeVacantArgs)
    _squareArgs.update(squareArgs)

    # Get the coordinates for the shortest path.
    if shortest:
        shortestEdges, shortestEdgesVertices = shortestPath(L, assignment)

    # Create subplots, turn axes off, set axis limits.
    fig, ax = plt.subplots()

    xlim, ylim = L.corners
    ax.set_xlim(-padding, xlim+padding)
    ax.set_ylim(-padding, ylim+padding)
    ax.set_axis_off()
    ax.set_aspect("equal")


    # Create a vertex map which specifies the possible embedded points each coordinate
    # can represent.
    vertexmap = {
        (x, y): list(
            product([x, xlim] if x == 0 else [x], [y, xlim] if y == 0 else [y])
        )
        for (x, y) in [c.encoding[0] for c in L.cells if len(c.encoding) < 2]
    }

    # Plot squares *first*. We need to check whether this is a torus (a periodic
    # cubical complex) as well, otherwise we end up plotting weird shit.
    squares = [c.encoding for c in L.cells[L.tranches[1]:]]

    for square in squares:
        possibleVertices = reduce(lambda A,B: A+B, [vertexmap[v] for v in square])
        possibleSquares = list(combinations(possibleVertices, r=4))

        for possibleSquare in possibleSquares:
            pairs = combinations(possibleSquare, r=2)
            dist = sum(
                1 for ((px, py), (qx, qy)) in pairs
                if (px == qx and abs(py-qy) == 1) or (py == qy and abs(px-qx) == 1)
            )
            
            if dist == 4:
                coordinates = list(sorted(possibleSquare))
                anchor = coordinates[0]
                rect = Rectangle(anchor, width=1, height=1, **_squareArgs)
                ax.add_patch(rect)

    # Plot edges next.
    edges = [c.encoding for c in L.cells[L.tranches[0]:L.tranches[1]]]
    nonzero = (assignment == 1).nonzero()[0]

    # Do it again for the shortest path.
    for j, ((ux, uy), (vx, vy)) in enumerate(edges):
        # No markers for edge ends.
        _edgeVacantArgs.update(dict(marker="none"))
        _edgeOccupiedArgs.update(dict(marker="none"))

        possibleVertices = list(product(vertexmap[(ux, uy)], vertexmap[(vx, vy)]))
        compatibleEdges = [
            ((ux, vx), (uy, vy)) for ((ux, uy), (vx, vy)) in possibleVertices
            if ((ux == vx and abs(uy-vy) == 1 and max(ux, vx) < L.corners[0])
                or (uy == vy and abs(ux-vx) == 1) and max(uy, vy) < L.corners[1])
        ]

        for x, y in compatibleEdges:
            if j in nonzero: ax.plot(x, y, **_edgeOccupiedArgs)
            else: ax.plot(x, y, **_edgeVacantArgs)

    # Do it again for the shortest path.
    if shortest:
        for j, ((ux, uy), (vx, vy)) in enumerate(shortestEdges):
            # No markers for edge ends.
            _edgeShortestArgs.update(dict(marker="none"))

            possibleVertices = list(product(vertexmap[(ux, uy)], vertexmap[(vx, vy)]))
            compatibleEdges = [
                ((ux, vx), (uy, vy)) for ((ux, uy), (vx, vy)) in possibleVertices
                if ((ux == vx and abs(uy-vy) == 1 and max(ux, vx) < L.corners[0])
                    or (uy == vy and abs(ux-vx) == 1) and max(uy, vy) < L.corners[1])
            ]

            for x, y in compatibleEdges: ax.plot(x, y, **_edgeShortestArgs)


    # Plot vertices *last*.
    others = list(set(range(len(L.skeleta[0]))) - (set(shortestEdgesVertices) if shortest else set()))

    px, py = zip(*[c.encoding[0] for c in L.cells[others]])
    ax.scatter(px, py, **_vertexOccupiedArgs)

    if shortest and shortestEdgesVertices:
        ix, iy = zip(*[c.encoding[0] for c in L.cells[shortestEdgesVertices]])
        ax.scatter(ix, iy, **_vertexShortestArgs)

    return fig, ax
