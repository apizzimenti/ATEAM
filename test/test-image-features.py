
import numpy as np
import matplotlib.pyplot as plt
from gerrytools.plotting import districtr

from potts import Lattice, SwendsonWang, Chain
from potts.stats import constant

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([20, 20], field=5)
model = SwendsonWang(temperature=constant(-np.log(3)))
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=10)

# Colors, for testing.
colors = districtr(L.field.order)

for state in chain.progress():
    # Create color assignments for each.
    vertexAssignment = [colors[spin] for spin in state]
    edgeAssignment = [e.spin for e in L.structure[1]]

    # Enable LaTeX.
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "text.latex.preamble": r"\usepackage{amsfonts}"
    })

    # Plot the lattice.
    ax = L.plot(
        edgeAssignment=edgeAssignment,
        vertexAssignment=vertexAssignment,
        vertexStyle=dict(marker="o", ms=4, markeredgewidth=0),
        edgeStyle=dict(lw=1, color="k"),
        #vertexLabels=True,
    )

    # Create a string for the lattice.
    latticeDesignation = ("times".join([str(c) for c in L.corners])).replace("times", "\\times")
    latticeDesignation = "$" + latticeDesignation + "$"
    latticeDesignation += f" sublattice of $\\mathbb{{Z}}^{len(L.corners)}$"

    # Create a string for the temperature designation and the step.
    temperature = f"$p = 1-e^{{-\\ln(3)}} = \\frac 23$"
    step = f"iteration ${chain.step}$"
    indicator = latticeDesignation + "\n" + temperature + "\n" + step

    # Add a little text box in the bottom-left.
    plt.text(0, -1/2, indicator, ha="left", va="top", fontsize=6)

    # Save the figure.
    plt.savefig(f"./output/figures/lattice-step-{chain.step-1}.png", dpi=300, bbox_inches="tight")
    plt.clf()

# After we're done, write the model's log to file.
with open("./output/log.txt", "w") as w:
    log = f"{'x'.join([str(t) for t in L.corners])} integer lattice over GF({L.field.order})\n\n"
    w.write(log + model.log)
