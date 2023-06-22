
import matplotlib.pyplot as plt
from gerrytools.plotting import districtr

from potts import Lattice, SwendsonWang, Chain

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([3, 3], field=3)
model = SwendsonWang(testing=True)
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=10)

# Colors, for testing.
colors = districtr(L.field.order)

for state in chain:
    # Create color assignments for each.
    vertexAssignment = [colors[spin] for spin in state]
    edgeAssignment = ["red" if e.spin == 1 else "k" for e in L.structure[1]]
    L.plot(
        edgeAssignment=edgeAssignment, vertexAssignment=vertexAssignment,
        vertexStyle=dict(marker="o", ms=20, markeredgewidth=0), edgeStyle=dict(lw=3),
        vertexLabels=True,
    )
    plt.savefig(f"./output/figures/lattice-step-{chain.step-1}.pdf", dpi=300, bbox_inches="tight")
    plt.clf()

# After we're done, write the model's log to file.
with open("./output/log.txt", "w") as w: w.write(model.log)