from TORCphysics import Circuit
import argparse

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This program runs a single simulation of the TORCphysics platform in the commmand line.
# It is basically just a parser that can run simple simulations taking the input from the command line.

# ---------------------------------------------------------------------------------------------------------------------
def main():
    """
    TORCphysics command-line interface.

    This function parses command-line arguments and runs a single genetic circuit simulation
    using the TORCphysics platform. It expects input files describing the circuit, sites,
    enzymes, and environment, and produces output files with simulation results.

    Command-line arguments:
        -c / --circuit       Path to circuit input CSV file (required)
        -s / --sites         Path to sites input CSV file (required)
        -e / --enzymes       Path to enzymes input CSV file (required)
        -n / --environment   Path to environment input CSV file (required)
        -o / --output        Output prefix for result files (default: "TORCsim")
        -f / --frames        Number of simulation frames (default: 3000)
        -r / --series        Series flag, if included, it prints the output dataframes and log file.
        -t / --dt            Time step in seconds (default: 1.0)

    Example:
        TORCphysics -c circuit.csv -s sites.csv -e enzymes.csv -n environment.csv -o out -f 3000 -t 1.0
    """

    #TODO: We still need to include the outputs.

    parser = argparse.ArgumentParser(description="TORCphysics: a physics-based simulation platform of supercoiling mediated regulation of gene expression in gene circuits.")
    #parser.add_argument("-c", "--continuation", action="store_true", help="Continuation of a simulation")
    parser.add_argument('-c', '--circuit', help='Circuit input csv file.', required=True)
    parser.add_argument('-s', '--sites', help='Sites input csv file.', required=True)
    parser.add_argument('-e', '--enzymes', help='Enzymes input csv file.',required=True)
    parser.add_argument('-n', '--environment', help='Environment input csv file.',required=True)
    parser.add_argument('-o', '--output', help='Output prefix for output files.',default='TORCsim')
    parser.add_argument('-f', '--frames', help='Number of frames (default 3000 frames).',type=int, default=3000)
    parser.add_argument('-r', '--series', action='store_true', help='Print output dataframes and log file.')
    parser.add_argument('-t', '--dt', type=float, help='Timestep (default 1 second).', default=1.0)

    # Parse infor
    args = parser.parse_args()

    continuation=False

    # Define circuit - All the related warnings should be displayed by he code itself
    my_circuit = Circuit(
        args.circuit, args.sites, args.enzymes,
        args.environment, args.output,
        args.frames, args.series, continuation, args.dt
    )

    if args.series:
        print("Series output and logging enabled.")

    # Print information system.
    my_circuit.print_general_information()

    my_circuit.run()

if __name__ == "__main__":
    main()


