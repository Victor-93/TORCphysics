# 🧪 TORCPhysics – Examples

This folder contains example files and Jupyter notebooks demonstrating how to use the TORCphysics package for simulating and analysing genetic circuits.

---

## 📁 Contents

### 🧬 Inputs

#### Circuit File
- **`circuit.csv`** – Circuit file containing general information about the genetic circuit.

#### Environment Files
- **`environment_1.csv`** – Contains basic RNAP model with constant velocity and unbinding upon reaching the termination site.
- **`environment_2.csv`** – Adds topoisomerases with continuous effects (non-binding) in all topological domains.
- **`environment_3.csv`** – Advanced RNAP model with multiple transcription stages (stalling) and stochastic topoisomerase activity 
(supercoiling-dependant binding). Only compatible with `sites_3.csv`.
- **`environment_4.csv`** – Same RNAP model as environment_3, but with continuous topoisomerase models.

#### Site Files
- **`sites_1.csv`** – Two-gene system in tandem orientation with constant (supercoiling-independent) binding rates.
- **`sites_2.csv`** – Same genes with supercoiling-dependent initiation using a sigmoid modulation function.
- **`sites_3.csv`** – Same system with multi-stage transcription model: `closed → open complex → elongation`. Only compatible with `environment_3.csv` and `environment_4.csv`.

---

### 📓 Jupyter Notebooks

- **`Example_1.ipynb`** – Basic simulation and analysis of a gene circuit.
- **`Example_2.ipynb`** – Multiple simulations with analysis and estimation methods.
- **`Example_3.ipynb`** – Using predefined models and custom parameter sets from TORCphysics.

---

## 📌 Note

These examples are intended for demonstration purposes. They can be used as templates for setting up and analysing your own genetic circuit simulations.
