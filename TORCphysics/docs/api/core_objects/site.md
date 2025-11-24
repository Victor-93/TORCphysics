# Site

The `Site` class defines a "binding" site. Environmentals have the capacity to bind specific sites.

---

## Usage Example 

Defines a custom site for IHF and adds it to my_circuit object

```python
from TORCphysics import Circuit, Site

my_circuit = Circuit(
    circuit_filename="circuit.csv",
    site_filename="sites.csv",
    enzyme_filename="enzymes.csv",
    environment_filename="environment.csv",
    frames=500,
    dt=1.0,
)
IHF_site = Site(site_type='IHF_site', name='IHF1', start=2500, end=0, k_on=0.0025, binding_model='PoissonBinding')
my_circuit.add_custom_Site(IHF_site)
```

::: TORCphysics.src.site
    options:
      show_root_heading: false
      show_source: false