# Environment

The `Environment` class defines a molecule in the virtual available environment. 

- Environmentals have the capicity to recognise binding sites.
- If they bind a binding site, they become `Effectors/Enzymes`.
- Environmentals have physical properties like sizes and concentrations.
- They also have binding, effect and unbinding models that describe their binding, effect on DNA and unbinding condition respectively.
- Currently, the environment does not deplete. In other words, the concentration remains constant.
---

## Usage Example 

Defines a custom IHF `Environmental` and appends it to the current circuit.

```python
from TORCphysics import Circuit, Environment

my_circuit = Circuit(
    circuit_filename="circuit.csv",
    site_filename="sites.csv",
    enzyme_filename="enzymes.csv",
    environment_filename="environment.csv",
    frames=500,
    dt=1.0,
)

IHF_environment = Environment(e_type='IHF', name='IHF', site_list=my_circuit.site_list, concentration=1.0, size=20,
                                effective_size=10, site_type='IHF_site',
                            unbinding_model='PoissonUnBinding')
my_circuit.add_custom_Environment(IHF_environment)
```

::: TORCphysics.src.environment
    options:
      show_root_heading: false
      show_source: false