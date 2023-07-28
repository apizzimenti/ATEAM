
import matplotlib.pyplot as plt
from gerrytools.plotting import districtr

from potts import GraphLattice, GraphSwendsonWang, Chain
from potts.stats import critical

# Create the Lattice, then instantiate the Swendson-Wang model.
L = GraphLattice([20, 20], field=3)
model = GraphSwendsonWang(temperatureFunction=critical(L.field.order))
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=5)

# Colors, for testing.
colors = districtr(L.field.order)

for state in chain.progress():
    # Create color assignments for each.
    vertexAssignment = { v.index: colors[v.spin] for v in L.graph.nodes() }
    edgeAssignment = { e.index: ("k" if e.spin else "#00000025") for e in L.graph.edges() }

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
        vertexStyle=dict(marker="o", ms=3, markeredgewidth=0),
        edgeStyle=dict(lw=1),
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
    plt.savefig(f"./output/figures/lattice-step-{chain.step-1}.png", dpi=300, bbox_inches="tight")
    plt.clf()

# After we're done, write the model's log to file.
with open("./output/log.txt", "w") as w:
    log = f"{'x'.join([str(t) for t in L.corners])} integer lattice over GF({L.field.order})\n\n"
    w.write(log + model.log)
