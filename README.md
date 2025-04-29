# TORCphysics: A Physical Model of DNA-Topology-Controlled Gene Expression

<p align="center">
  <img src="TORCphysics/logo.svg" alt="TORCphysics Logo" width="200"/>
</p>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/torcphysics.svg)](https://pypi.org/project/torcphysics/)
[![GitHub](https://img.shields.io/badge/Source-GitHub-black?logo=github)](https://github.com/Victor-93/TORCphysics)

---
 
**TORCphysics** is a physics-based simulation framework to model gene expression regulated through 
DNA supercoiling. It simulates interactions between DNA-binding proteins such as RNA polymerases and 
topoisomerases, capturing both physical DNA properties and dynamics. 

Transcription is modeled based on the twin-supercoiling domain model, and the simulation supports 
supercoiling-sensitive and non-sensitive promoters, as well as multi-stage binding kinetics. 
Outputs include transcription rates and time-series data, enabling direct comparisons with experimental 
results or kinetics inference via parameter search.

---

## 📦 Installation

### Latest version from PyPI:
```bash
pip install torcphysics
```

### TORCphysics paper version for reproducing results:
```bash
pip install git+https://github.com/Victor-93/TORCphysics.git@TORCphysics_paper
```

---

## 🚀 Quick Start

TORCphysics can be used from the command line or directly in Python scripts.

Each simulation requires four interconnected input files:

| Input                 | Description                                                                                                                                                                        |
|-----------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `circuit.csv`         | General information of gene circuit, such as open (linear, .e.g., chromosomal) or closed (circular e.g., plasmid)  structure,and the initial superhelical density.. **(Required)** |
| `environment.csv`     | List of DNA-binding molecules (e.g., RNAPs, topoisomerases, NAPs), that can bind particular DNA sites such as promoters or protein binding sites.                                  |
| `sites.csv`           | List of binding sites on the DNA (e.g., promoters, protein binding sites).                                                                                                         |
| `enzymes.csv`         | Initial state: DNA-bound molecules present at the beginning of the simulation.                                                                                                     |

> **⚠️ Inputs warnings**  
> The only required input to run a simulation is `circuit.csv`. 
> However, if `sites.csv` is not provided, molecules from the environment will not bind to the DNA. 
> If `environment.csv` is missing, no molecules will bind to the sites. 
> Lastly, if `enzymes.csv` is not provided, the simulation will start without any pre-bound molecules, 
> but it will still run without issues.

---

### 🖥️ Command-Line Usage

Run a simulation with command:

```bash
TORCphysics -c circuit.csv -s sites.csv -e enzymes.csv -n environment.csv -o out -f 3000 -t 1.0 -r
```

This will simulate:
- `3000` frames with `1.0` second time step of the system described by the csv files.
- Outputs will be saved as dataframes and log files with prefix "out".

For more information type:

```bash
TORCphysics --help
```

---

### 🐍 Scripting Usage

TORCphysics can also be used via Python scripting for custom simulations:

```python
from TORCphysics import Circuit

my_circuit = Circuit(
    circuit_filename="circuit.csv",
    sites_filename="sites.csv",
    enzymes_filename="enzymes.csv",
    environment_filename="environment.csv",
    output_prefix="out",
    frames=3000,
    series=True,
    continuation=False,
    dt=1.0
)

my_circuit.print_general_information()
my_circuit.run()
```

---

## 📚 Examples

Example files and Jupyter notebooks are available in the [`Examples/`](TORCphysics/Examples/) directory:

- `Example_1.ipynb` — Single gene simulations and analysis.
- `Example_2.ipynb` — Multiple simulations with statistical analysis.
- `Example_3.ipynb` — Defining custom enzyme/site models using built-in models.

---

## ⚙️ Algorithms

> Detailed documentation and algorithmic explanations coming soon on the [Wiki](https://github.com/Victor-93/TORCphysics/wiki)

---

## 📖 Citation

If you use TORCphysics in your research, please cite the paper below:

TORCphysics paper coming soon!

---

## 📬 Contact

Questions or collaborations?

📧 V.VelascoBerrelleza@sheffield.ac.uk

---

## 🔓 License

This project is licensed under the **GNU General Public License v3.0** (GPLv3).  
See the [LICENSE](./LICENSE) file for more details.
```
