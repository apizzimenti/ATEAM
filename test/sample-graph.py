
from potts.models import GraphIsing, GraphSwendsonWang, GraphPercolation
from potts.structures import GraphLattice
from potts.stats import critical, metropolis, linear, randomizedToConstant
from potts import Chain
import math
import numpy as np

steps = 1000
q = 2
GL = GraphLattice([20, 20], field=q)
model = GraphPercolation(
    temperatureFunction=randomizedToConstant(
        -math.log(1/2), steps, distribution=np.random.uniform
    )
)
initial = model.initial(GL)

# Create and run the chain.
chain = Chain(GL, model, initial, accept=metropolis, steps=steps)

for state in chain.progress(): pass
