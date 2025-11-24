# Circuit

The `Circuit` class defines a genetic circuit and manages simulations.  
It connects TORCphysics workflows (binding, effect, unbinding) and tracks time-series data.

---

## Usage Example

```python
from TORCphysics.src.circuit import Circuit

c = Circuit(
    circuit_filename="circuit.csv",
    site_filename="sites.csv",
    enzyme_filename="enzymes.csv",
    environment_filename="environment.csv",
    frames=500,
    dt=1.0,
)

c.run() 
```

::: TORCphysics.src.circuit
options:
show_root_heading: false
show_source: false