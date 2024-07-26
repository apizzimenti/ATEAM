
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
from matplotlib.patches import Rectangle
from itertools import combinations, product


def lattice2D(
        L, assignment, padding=0.1, vertexArgs=dict(s=40, color="k", zorder=10),
        edgeOccupiedColor="#3B444B", edgeOccupiedWidth=1.5, edgeVacantColor="#3B444B10",
        edgeVacantWidth=1, edgeArgs=dict(zorder=0),
        squareArgs=dict(alpha=1/2, facecolor="#87A96B", edgecolor="none", zorder=0)
    ):
    """
    """
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
    last = np.cumsum(list(L.skeleta.values()))
    squares = [c.encoding for c in L.cells[last[1]:]]

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
                rect = Rectangle(anchor, width=1, height=1, **squareArgs)
                ax.add_patch(rect)

    # Plot edges next.
    edges = [c.encoding for c in L.cells[last[0]:last[1]]]
    nonzero = (assignment == 1).nonzero()[0]

    for j, ((ux, uy), (vx, vy)) in enumerate(edges):
        # No markers for edge ends.
        edgeArgs.update(dict(marker="none"))

        if j in nonzero:
            edgeArgs.update(dict(color=edgeOccupiedColor))
            edgeArgs.update(dict(linewidth=edgeOccupiedWidth))
        else:
            edgeArgs.update(dict(color=edgeVacantColor))
            edgeArgs.update(dict(linewidth=edgeVacantWidth))

        possibleVertices = list(product(vertexmap[(ux, uy)], vertexmap[(vx, vy)]))
        compatibleEdges = [
            ((ux, vx), (uy, vy)) for ((ux, uy), (vx, vy)) in possibleVertices
            if ((ux == vx and abs(uy-vy) == 1) or (uy == vy and abs(ux-vx) == 1))
        ]

        for x, y in compatibleEdges:
            ax.plot(x, y, **edgeArgs)

    # Plot vertices *last*.
    vx, vy = zip(*[c.encoding[0] for c in L.cells[:L.skeleta[0]]])
    ax.scatter(vx, vy, **vertexArgs)

    return fig, ax
