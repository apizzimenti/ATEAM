
import numpy as np
import matplotlib.pyplot as plt
from gerrytools.plotting import districtr, latex

from potts import GraphLattice, GraphSwendsonWang, GraphPercolation, Chain
from potts.stats import critical

# Create the Lattice, then instantiate the Swendson-Wang model.
L = GraphLattice([10,10], field=2)
model = GraphSwendsonWang(temperatureFunction=critical(L.field.order))
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=4)

# Colors, for testing.
colors = districtr(L.field.order)


for state in chain.progress():
    chain.lattice.assign(state)
    vertexAssignment = { v.index: latex["Amber"] if v.spin == 0 else latex["Cadmium Green"] for v in L.graph.nodes()}
    edgeAssignments = [
        { e.index: "#00000025" for e in L.graph.edges() },
        { e.index: ("k" if (state[e.at[0].index] == state[e.at[1].index] and np.random.uniform() < 0.6) else "#00000025") for e in L.graph.edges() }
    ]

    for index, edgeAssignment in enumerate(edgeAssignments):
        # Enable LaTeX.
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "text.latex.preamble": r"\usepackage{amsfonts}"
        })

        # Plot the lattice.
        ax = L.plot(
            L.graph,
            edgeAssignment=edgeAssignment,
            vertexAssignment=vertexAssignment,
            vertexStyle=dict(marker="o", ms=10, markeredgewidth=0),
            edgeStyle=dict(lw=3),
            #vertexLabels=True,
        )

        # Create a string for the lattice.
        latticeDesignation = ("times".join([str(c) for c in L.corners])).replace("times", "\\times")
        latticeDesignation = "$" + latticeDesignation + "$"
        latticeDesignation += f" sublattice of $\\mathbb{{Z}}^{len(L.corners)}$"

        # Create a string for the temperature designation and the step.
        field = f", coefficients in $\\mathbb{{F}}_{L.field.order}$"
        temperature = f"$p = 1-e^{{-\\ln(3)}} = \\frac 23$"
        step = f"iteration ${chain.step}$"
        indicator = latticeDesignation + field + "\n" + temperature + "\n" + step

        # Add a little text box in the bottom-left.
        #  plt.text(0, -1/2, indicator, ha="left", va="top", fontsize=6)

        # Save the figure.
        plt.savefig(f"./output/figures/lattice-step-{chain.step-1}-{index}.png", dpi=300, bbox_inches="tight", transparent=True)
        plt.clf()

# After we're done, write the model's log to file.
with open("./output/log.txt", "w") as w:
    log = f"{'x'.join([str(t) for t in L.corners])} integer lattice over GF({L.field.order})\n\n"
    w.write(log + model.log)
