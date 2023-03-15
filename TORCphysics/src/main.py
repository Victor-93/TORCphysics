import argparse
import sys

import params

from TORCphysics import Circuit

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This program simulates a genetic circuit under certain conditions given the inputs:
# circuit.csv, sites.csv, enzymes.csv and environment.csv
# According these inputs, RNAPs will stochastically bind the DNA and will generate
# supercoiling accordingly.
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------

# All parameters are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma
dt = params.dt

# ---------------------------------------------------------------------------------------------------------------------
# INPUTS
# ---------------------------------------------------------------------------------------------------------------------

# Create the parser
# ---------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Version1 of the physical model of transcription-supercoiling")
parser.add_argument("-f", "--frames", type=int, action="store", help="Number of frames (timesteps)", default=5000)
parser.add_argument("-c", "--continuation", action="store_true", help="Continuation of a simulation")
parser.add_argument("-ic", "--input_circuit", action="store", help="Circuit input file", default="../circuit.csv")
parser.add_argument("-ig", "--input_sites", action="store", help="Genome input file", default="../sites.csv")
parser.add_argument("-io", "--input_enzymes", action="store", help="Objects input file", default="../enzymes.csv")
parser.add_argument("-ie", "--input_environment", action="store", help="Environment input file",
                    default="../environment.csv")
parser.add_argument("-s", "--series", action="store_true", help="Print dynamic results per timestep")
parser.add_argument("-t", "--test", action="store_true", help="Run series of tests stored in the test folder")
parser.add_argument("-o", "--output", action="store", help="Output prefix for output files", default="output")
parser.add_argument("-dt", "--timestep", action="store", help="Simulation time step", default=dt)
parser.add_argument("-tm", "--topoisomerase_model", action="store", help="Model for topoisomerase activity",
                    default='continuum')
parser.add_argument("-mm", "--mechanical_model", action="store", help="Model for enzyme mechanics", default='uniform')

# Process terminal commands
# ---------------------------------------------------------------------------------------------------------------------
args = parser.parse_args()

# Let's put this information into variables
frames = args.frames  # Number of frames
circuit_filename = args.input_circuit  # Circuit file
sites_filename = args.input_sites  # Sites file
enzymes_filename = args.input_enzymes  # Enzymes file
environment_filename = args.input_environment  # Environment file
output_prefix = args.output  # Output file prefix
dt = abs(args.timestep)  # The time step
topoisomerase_model = args.topoisomerase_model
mechanical_model = args.mechanical_model

if args.continuation:  # If this is the continuation of a previous run
    continuation = True
else:
    continuation = False
if args.series:  # If series, then it will print dynamic output files (.txt)
    series = True
else:
    series = False
if args.test:  # If true, then we run some tests and stop
    test = True
else:
    test = False

# Pass the command line inputs, read csvs and initialize "my_circuit"
my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                     output_prefix, frames, series, continuation, dt, topoisomerase_model, mechanical_model)

# TODO: Don't forget to do your tests
# TODO: Events - series?
# TODO: Log
# TODO: Output
# TODO: Print output information
# Let's print some info of the system
my_circuit.print_general_information()
print(my_circuit.site_list[0].name)

# Now run
my_circuit.run()
sys.exit()
