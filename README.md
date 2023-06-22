## Potts lattice gauge theory

### Installation
Clone this repository and run `sh setup.sh` in your favorite terminal. This
installs the `potts` package globally.

### Usage
Use the `test/` directory for testing code to go into the package. Sample code
is below. If the `testing=True` flag is passed to the call to `SwendsonWang()`,
then output directories are created (though they should already exist) where the
`SwendsonWang` class expects them to exist.

```python
from potts import Lattice, SwendsonWang, Chain

# Create the Lattice, then instantiate the Swendson-Wang model.
L = Lattice([3, 3], field=3)
model = SwendsonWang()
initial = model.initial(L)

# Create the chain.
chain = Chain(L, model, initial, steps=10)

# Iterate over the chain.
for state in chain:
    <do whatever>
```