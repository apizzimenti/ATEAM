# `potts`
A high-performance(-ish) tool for running experiments on the Potts model.

## Installation
Clone this repository by `git clone https://github.com/apizzimenti/potts.git`
and run `sh setup.sh` in your favorite terminal. This installs the `potts` package
(and its dependencies) globally in development mode, so any changes you make to
the source are reflected system-wide.

## Usage
Use the `test/` directory for testing code to go into the package. Sample code
is below. If the `testing=True` flag is passed to the call to `SwendsonWang()`,
then output directories are created (though they should already exist) where the
`SwendsonWang` class expects them to exist.

```python
from potts import Lattice, GraphSwendsonWang, Chain

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([32, 32], field=2)
model = GraphSwendsonWang()
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=10)

# Iterate over the chain.
for state in chain:
    <do whatever>
```

Using the imaging features, we can make pretty outputs, like this simulation of
the Ising model on a 32x32 lattice at _just_ below the critical temperature.

<p align="center">
  <img src="test/output/figures/lattice.gif" alt="animated" />
</p>