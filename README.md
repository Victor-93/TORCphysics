# TORCphysics: A Physical Model of DNA-Topology-Controlled Gene Expression
<p align="center">
<img src="TORCphysics/logo.svg" />
</p>

And the other thingyis with the links to PyPI, paper, and something else.

[TORCphysics](https://github.com/Victor-93/TORCphysics) is a physics-based simulation 
platform for simulating the interactions of DNA-binding molecules through DNA supercoiling, with a 
particular emphasis on transcription. Transcription is described according the twin-supercoiling 
domain model, and the model can account for the stochastic binding of DNA 
topoisomerases (topo I  and gyrase). Outputs can be directly compared with experimental results 
by outputing transcription rates, and parameters can be inferred through parameter exploration 
and comparing with experimental data.

## Installation and Versions

Install the latest version of **TORCphysics**:
```ruby
pip install TORCphysics
```
Install the TORCphysics paper version with scripts for reproducing experiments:
```ruby
pip install git+https://github.com/Victor-93/TORCphysics.git@sub_models
```


## Usage

TORCphysics can be quickly used through the command line or scripting, 
providing a higher flexibility and verstility for manioulating simulations.

In either case, to run TORCphysics it needs four inputs that are inter-connected: 
**environment, sites and enzymes**.
* The **Environment** input is composed by the DNA-binding molecules that have the capacity of binding
particular DNA **sites** such as promoters or protein binding sites. 
* **Sites** represent the binding sites of particular molecules, for example, RNA polymerases can bind 
gene start sites (promoters), while proteins bind their indicated binding sites.
* **Enzymes** correspond to the bound molecules (e.g., enzymes and proteins) to the DNA at 
the start of the simulation.
* **Circuit**  input represents overall information about the simulated system, such as
size, open (linear, .e.g., chromosomal) or closed (circular e.g., plasmid)  structure, 
and the initial superhelical density. 

These inputs can be defined through csv files, e.g., **environment.csv, 
sites.csv, enzymes.csv and circuit.csv**, or through scripting using TORCphysics 
libraries (see example REF). 

For more detailed information regarding inputs and outputs can be found in REF.

> **⚠️ Inputs warnings**  
* The only required input to run a simulation is the **circuit**, however if **sites**
are not provided, the molecules in the environment will not bind the DNA, 
and if the **environment** is missing, no molecules will bind sites. 
Lastly, if **enzymes** is not provided, the code will not have any starting
bound molecules but simulation can progress without problems.

#### Command Line Usage
To quickly run TORCphysics, do the following command:

```ruby
TORCphysics -c circuit.csv -s sites.csv -e enzymes.csv -n environment.csv -o out -f 3000 -t 1.0 -r
```

This will produce a simulation using the input system described through the csv files, 
will run a simulation with 3000 frames with 1.0 seconds time-step, and 
outputing dataframes and log files with prefix name "out". You can run this command
in the example directory with the example csv files. 

You can get more information regarding inputs and outputs with:

```ruby
TORCphysics --help
```

#### Scripting Usage


## Examples

Hay que hacer varios ejemplitos con el mismo sistema.

* Ejemplo 1 basico single sim, correr simulacion, cargar DFs, plotear senales, y producir animacion. Hacer para continuum model of topos y stochastic.En este hagamos todo: 
  * 1.- Stochastic initiation, RNAP sin topos ni stalling. RNAP con stalling simple.
  * 2.- Sigmoid initiation, RNAP uniform y topos continuum
  * 3.- Sigmoid inititaion, RNAP stall y stochastic topos
* Ejemplo 2 simulaciones multiples: correr muchas simulaciones, no producir dfs pero mantenerlos en el codigo, y analyzar correlaciones, production rates, local supercoiling, y global supercoiling
* Ejemplo 3 Definir lac model en medio, y hacer lo mismo que el ejemplo 2.

## Algorithms or functions?


## Citations
If you use this code, please site the paper below:

TORCphysics paper coming soon!

## Contact Information

## License