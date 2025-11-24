# Enzyme

The `Enzyme` class defines a DNA-bound molecules. These once were `Environmentals` that became `Effectors/Enzymes` upon binding.  

- These bound molecules have may move along the DNA, act as topological barriers and modify the superhelicity/twist. 
- They inherit the physical propierties of their corresponding environmentals, in terms of size and effect models.

---

## Usage Example 

Defines a custom `Enzyme` and appends it to the current circuit.

Here, we do not provide an Effect Model, hence the bound IHF do not change the superhelical density but acts as a topological barrier.

```python
from TORCphysics import Circuit, Enzyme

my_circuit = Circuit(
    circuit_filename="circuit.csv",
    site_filename="sites.csv",
    enzyme_filename="enzymes.csv",
    environment_filename="environment.csv",
    frames=500,
    dt=1.0,
)
# Note that the twist shpi
IHF_enzyme = Enzyme(e_type='IHF', name='IHF', site='IHF_site', size=100, effective_size=50, position=600,
                         twist=0.0, superhelical=0.0, unbinding_model_name='PoissonUnBinding')
my_circuit.add_new_enzymes([IHF_enzyme])
```

::: TORCphysics.src.enzyme
    options:
      show_root_heading: false
      show_source: false