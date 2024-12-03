
from random import sample
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from functools import reduce
from itertools import combinations, product
import phat

vertexmap = {
    0: [(0,0), (3,0), (0,3), (3,3)],
    1: [(1,0), (1,3)],
    2: [(2,0), (2,3)],
    3: [(0,1), (3,1)],
    4: [(1,1)],
    5: [(2,1)],
    6: [(0,2), (3,2)],
    7: [(1,2)],
    8: [(2,2)]
}

vertices = [(0, [])]*9

edges = [
    (1, [0, 1]),
    (1, [0, 3]),
    (1, [0, 2]),
    (1, [0, 6]),
    (1, [1, 2]),
    (1, [1, 4]),
    (1, [1, 7]),
    (1, [2, 5]),
    (1, [2, 8]),
    (1, [3, 4]),
    (1, [3, 6]),
    (1, [3, 5]),
    (1, [4, 5]),
    (1, [4, 7]),
    (1, [5, 8]),
    (1, [6, 7]),
    (1, [6, 8]),
    (1, [7, 8]),
]

squares = [
    (2, [9, 10, 14, 18]),
    (2, [13, 14, 16, 21]),
    (2, [18, 19, 22, 24]),
    (2, [21, 22, 23, 26]),
    (2, [10, 11, 16, 20]),
    (2, [9, 12, 15, 24]),
    (2, [19, 20, 23, 25]),
    (2, [13, 15, 17, 26]),
    (2, [11, 12, 17, 25]),
]

edgeIndices = list(range(len(vertices), len(vertices)+len(edges)))
edgeReordering = list(zip(edgeIndices, edges))
edgeReordering = sample(edgeReordering, len(edgeReordering))
indices = [i for i, _ in edgeReordering]
edges = [e for _, e in edgeReordering]
edgeMap = { i: k for i, k in zip(indices, edgeIndices) }

squares = [ (d, list(sorted([edgeMap[f] for f in faces]))) for d, faces in squares ]


# To compare the resultant persistence pairs, we shuffle the order in which the
# indices of each square's faces are listed. These indices are listed in increasing
# order by default.

# questions = ["original", "shuffled"]
# shufflers = [lambda A: A, lambda B: sample(B, len(B))]
questions = ["original"]
shufflers = [lambda A: A]

for question, shuffler in zip(questions, shufflers):
    # Shuffle the ordering of the squares' faces and construct the cubical complex.
    squares = [(d, shuffler(faces)) for d, faces in squares]
    complex = vertices + edges + squares

    # Construct the boundary matrix and compute the persistence pairs.
    boundary = phat.boundary_matrix()
    boundary.columns = complex
    pairs = boundary.compute_persistence_pairs()
    pairs.sort()
    pairs = list(pairs)

    # Determine the births of essential cycles; these should be invariant under
    # re-orderings of the squares' faces in the specification.
    deaths = set([d for b, d in pairs])
    births = set([b for b, d in pairs])
    times = set(range(len(complex)))
    essential = times-(births|deaths)

    birthmap = { b: d for b, d in pairs }
    deathmap = { d: b for b, d in pairs}

    print(question)
    print(f"births @ {births}")
    print(f"deaths @ {deaths}")
    print(f"essential @ {essential}")
    print()

    for k in range(len(complex)):
        subcomplex = complex[:k+1]

        fig, ax = plt.subplots()
        ax.set_xlim(-0.1, 3.1)
        ax.set_ylim(-0.1, 3.1)
        ax.set_axis_off()

        for j in range(len(subcomplex)):
            d, faces = subcomplex[j]
            ##############
            ## VERTICES ##
            ##############
            if d == 0:
                compatibleVertices = [(x, y) for x, y in vertexmap[j] if x < 3 and y < 3]
                X, Y = zip(*compatibleVertices)
                ax.scatter(X, Y, s=40, color="k", zorder=10000)

            ###########
            ## EDGES ##
            ###########
            if d == 1:
                u, v = faces

                possibleVertices = list(product(vertexmap[u], vertexmap[v]))
                compatibleEdges = [
                    ((ux, vx), (uy, vy)) for ((ux, uy), (vx, vy)) in possibleVertices
                    if ((ux == vx and abs(uy-vy) == 1) or (uy == vy and abs(ux-vx) == 1))
                ]

                for x, y in compatibleEdges:
                    ax.plot(x, y, lw=4, color="blue" if j == k else "red", zorder=0)

            #############
            ## SQUARES ##
            #############
            if d == 2:
                edgeVertices = set(reduce(lambda A,B: A+B, [subcomplex[f][1] for f in faces]))
                possibleVertices = reduce(lambda A,B: A+B, [vertexmap[v] for v in edgeVertices])
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
                        square = Rectangle(anchor, width=1, height=1, edgecolor="none", facecolor="cyan", zorder=-1)
                        ax.add_patch(square)

        if k in births: title = f"born at {k}, dies at {birthmap[k]}"
        if k in deaths: title = f"born at {deathmap[k]}, dies at {k}"
        if k in essential: title = f"essential born at {k}"
        
        ax.set_title(f"step {k} --- {title}")
            
        plt.show()