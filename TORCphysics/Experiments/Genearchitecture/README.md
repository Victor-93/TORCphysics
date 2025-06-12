# Gene Architecture Experiments

This folder contains scripts for simulating genetic systems using **TORCphysics**. 
It includes tools for parameterising the melting energy required to form the open complex via the **SIST** (Stress-Induced Duplex Destabilization (SIDD)) algorithm, 
and for simulating gene expression and susceptibility from the single-gene experiment in *Boulas et al.*
(*"Assessing in vivo the impact of gene context on transcription through DNA supercoiling"*)

---

## 🛠️ General Instructions

### 🧹 Preprocessing – SIST Analysis

1. Navigate to `SIST_on_promoters/`.
2. Read the `README` in that folder for specific instructions.
3. Run `sequence_generator.py`, then move and rename the generated sequence files to `promoter/SIST_code/` as `input.fasta`.
4. Run the `run_comp_n.sh` bash scripts *or* copy the commands in the `COMMAND` file into your terminal.
5. Navigate to `filtering/` and execute the bash scripts.
6. Go to `plot/response_curve/` and run `fit_function.py` to fit a sigmoid curve to the melting energy. This function becomes an input to the promoter/gene models.
7. *(Optional)* Return to `SIST_code/` and run `plot_energies.py` and `responses_fig.py` to visualize the energy profiles and fitted responses.

---

### ⚙️ Calibration

1. Navigate to `calibration/`.
2. Run `calibrate-genearch_simple_VN.py` for each model set `V0`, `V1`, and `V2`. Each script performs a random search calibration and outputs:

   * `*-trials.pkl`: Information for each iteration,
   * `*values.csv`: Parameters and corresponding losses.
3. Repeat the calibration for each promoter (`weak`, `medium`, `strong`).

---

### 📊 Analysis

1. Navigate to `analysis/`.
2. Run `plot_loss.py` – visualizes loss distributions from the calibrations.
3. Run `plot_genearch-figs.py` – performs overall analysis, plotting:

   * Susceptibility,
   * Number of bound enzymes,
   * Production rate vs. loss,
   * Superhelical densities.
4. Run `plot_responses.py` – visualizes the promoter response curves for each model set.

---

## 🔄 Calibration Workflow

### Step 1 – Modeling Promoter Melting Energy

Before calibrating, the SIST algorithm is applied to each promoter sequence (weak, medium, strong) to generate destabilization profiles. 

* Synthetic sequences are built for each promoter, with flanked by 250 G's on both sides.
* The SIST algorithm is run at varying superhelical densities.
* The average melting energy is computed within the promoter region.
* A sigmoid curve is fitted to represent melting energy as a function of supercoiling.
* This curve modulates RNAP binding (V0) or the open-complex formation (V1, V2).

### Step 2 – Calibration Against Experimental Data

Using *Boulas et al.* susceptibility data:

* A random search algorithm is used to vary promoter parameterisation for each model set and promoter.
* For each, case, 63 simulations are run per distance.
* Susceptibility is calculated relative to a reference distance.
* The simulation susceptibility is compared with experimental susceptibility and its error.
* A loss function quantifies the deviation.

---

## 🧬 Enzyme Models

### 🧪 V0 Model Set

**Topoisomerases**

* *Binding*: Recognition
* *Effect*: Linear
* *Unbinding*: Poisson

**RNA Polymerase (RNAP)**

* *Binding*: SISTBasedBinding – based on melting energy
* *Effect*: RNAPStall – immediately advances after binding
* *Unbinding*: Simple – unbinds at the termination site

---

### 🧪 V1 Model Set

**Topoisomerases**

* *Binding*: Recognition
* *Effect*: Linear
* *Unbinding*: Poisson

**RNA Polymerase (RNAP)**

* *Binding*: Gaussian – elastic binding
* *Effect*: RNAP Stages + Stalling – multi-step binding with potential stalling during elongation
* *Unbinding*: RNAP Stages + Simple – spontaneous unbinding before elongation or upon reaching termination

---

### 🧪 V2 Model Set

**Topoisomerase I**

* *Binding*: Recognition + RNAP Tracking
* *Effect*: Linear
* *Unbinding*: Poisson

**RNA Polymerase (RNAP)**

* *Binding*: Gaussian – elastic binding 
* *Effect*: RNAP Stages + Stalling – multi-step binding with potential stalling during elongation
* *Unbinding*: RNAP Stages + Simple – spontaneous unbinding before elongation or upon reaching termination

---

## 📦 Requirements (in addition to TORCphysics)

* `numpy`
* `pandas`
* `hyperopt`

---

## 📁 Repository Structure

| Path                                          | Description                                                            |
| --------------------------------------------- |------------------------------------------------------------------------|
| `SIST_on_promoters/`                          | Scripts to perform SIST analysis and extract melting energy (`U_melt`) |
| `promoter_responses/`                         | Preprocessed promoter response curves from SIST                        |
| `junier_data/`                                | Experimental susceptibility data for each promoter                     |
| `calibration/calibrate-genearch_simple_V0.py` | Calibration script for V0 model                                        |
| `calibration/calibrate-genearch_simple_V1.py` | Calibration script for V1 model                                        |
| `calibration/calibrate-genearch_simple_V2.py` | Calibration script for V2 model                                        |
| `calibration_environment_avgx1_dt1.0.csv`     | Environment file for V0/V1 (no RNAP tracking)                          |
| `calibration_environment_avgx2_dt1.0.csv`     | Environment file for V2 (with RNAP tracking)                           |
| `analysis/plot_loss.py`                       | Visualizes calibration loss distributions                              |
| `analysis/plot_genearch-figs.py`              | Performs overall analysis                                              |
| `analysis/plot_responses.py`                  | Visualizes promoter response functions by model and promoter           |


