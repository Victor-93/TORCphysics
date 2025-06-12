# Topo Kinetics Experiments

This folder contains Python scripts for calibrating the stochastic activity of topoisomerases using kinetic parameters from the experimental study by Wang et al., *"Kinetic Study of DNA Topoisomerases by Supercoiling-Dependent Fluorescence Quenching"*, in *E. coli*.

---

## 📋 General Instructions

### 🔧 Calibration

1. *(Optional)* Navigate to `kinetic_curves/` and run `plot_kinetic_curves.py` to visualise the experimental reaction curves.
2. Navigate to `calibration/` and run `calibrate_topoisomerases.py` to start the calibration process.

### 📊 Analysis

1. Navigate to `analysis/`.
2. Run `calculate_averages.py` to estimate the topoisomerase model parameters from calibration results.
3. Run `plot_loss.py` to analyze the distribution of calibration losses.
4. Run `plot_calibration.py` to visualise global superhelical densities, fitting errors, and the binding/effect model outputs.
5. Run `plot_steady-state.py` to visualise the steady-state behavior (number of bound enzymes).

---

## ⚙️ Calibration Information

### 🧪 Pre-processing: Reference Superhelical Density Curves

Before the calibration, we simulate the Michaelis-Menten reaction kinetics to infer global superhelical density. We create four scenarios:

1. Topoisomerase I acting on negatively supercoiled DNA
2. Gyrase acting on relaxed DNA
3. Both enzymes acting on relaxed DNA
4. Both enzymes acting on negatively supercoiled DNA

### 🧮 Calibration Procedure

Using the reference curves, we calibrate the topoisomerase models by:

* Running multiple parallel simulations for each scenario with varying parameter sets.
* Averaging the resulting global superhelical density for each parameter set.
* Comparing simulated curves with experimental references.
* Optimizing parameters using a random search algorithm to minimize combined error across all four scenarios.

### 🔬 Topoisomerase Models

**Topoisomerase I:**
* Binding model* = Recognition Binding
* *Effect model* = Linear effect
* *Unbinding model* = Poisson unbinding

**Gyrase:**
* *Binding model* = Recognition Binding
* *Effect model* = Linear effect
* *Unbinding model* = Poisson unbinding

---

## 📦 Requirements (in addition to TORCphysics)

* `hyperopt`
* `pandas`
* `numpy`

---

## 📁 Repository Structure

| Path                                      | Description                                                                                    |
| ----------------------------------------- |------------------------------------------------------------------------------------------------|
| `kinetic_curves/plot_kinetic_curves.py`   | visualises the experimental kinetic reaction curves for all four scenarios.                    |
| `calibration/calibrate_topoisomerases.py` | Performs the full calibration workflow. Outputs best parameter sets and iteration-wise values. |
| `circuit.csv`                             | Input TORCphysics circuit file.                                                                  |
| `environment.csv`                         | Input TORCphysics environment file.                                                            |
| `analysis/calculate_averages.py`          | Processes losses and averages the top 10% of best-fitting parameter sets.                      |
| `analysis/plot_loss.py`                   | Plots the histogram of losses from the calibration process.                                    |
| `analysis/plot_steady-state.py`           | Re-simulates with average parameter sets and visualises steady-state.                  |
| `analysis/plot_calibration.py`            | Re-simulates and plots global superhelical density curves using estimated parameter sets.      |


