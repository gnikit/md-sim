import os
import sys
import getopt
import lxml.etree as ET


USAGE = "usage: generate_xml_file.py " \
        "[-d dir] [-s steps] [-c compress] [-r rdf_bins]" \
        "[-p particles] [-l lattice] [-t tracking] [-w rdf_wait]\n" \
        "                            [-R rho] [-T temperature] [-n=strength]" \
        "[-A softening_param] [-P pair_potential]\n" \
        "                            [-o output_file]\n\n"

EXAMPLE = "generate_xml_file.py -d ./examples/example_data -s 10000 -c false"\
          " -r 500 -p 10 -l SC -t false -w 2000 --rho 0.5 -T 0.5 -n 8 -A 0.5"\
          " --pp BIP -o ./MD-simulation/examples/example\n\n"


HELP = "Options:\n" \
       "  -d, --dir:         Output directory for simulation data.\n" \
       "  -s, --steps:       Number of steps for the simulation.\n" \
       "  -c, --compress:    Whether the fluid will be compressed.\n" \
       "                     Set to true when testing for phase changes.\n" \
       "                     default is false.\n" \
       "  -r, --rdf_bins:    Histogram bins used to resolve the accuracy\n" \
       "                     of the Radial Distribution Function (RDF)\n" \
       "                     default is 500.\n" \
       "  -p, --particles:   Particles per axis of the cubic box.\n" \
       "                     default is 10 per axis.\n" \
       "  -l, --lattice:     Lattice selection for the simulation.\n" \
       "                     Initial lattice selection:\n" \
       "                         \"FCC\": Face Centred Cubic\n" \
       "                         \"BCC\": Body Centred Cubic\n" \
       "                          \"SC\": Simple Cubic\n" \
       "  -t, --track:       true/false Enables/disables particle tracking.\n" \
       "                     set to true and use visualise.py. \n" \
       "                     Default is false.\n" \
       "  -w, --rdf_wait:    Number of steps the RDF will wait until it \n" \
       "                     starts collecting data. The steps input in\n" \
       "                     the `rdf_wait` should be greater than the \n" \
       "                     total number of steps used in the simulation.\n" \
       "  -R, --rho:         The number density of the fluid.\n" \
       "  -T, --temp:        The temperature of the fluid.\n" \
       "  -n, --n:           The pair potential strength of the fluid.\n" \
       "  -A, --A:           The pair potential softening parameter \n" \
       "                     (where applicable).\n" \
       "  -P, --pp:          The type of pair potential used by the fluid.\n" \
       "                     The present types are:\n" \
       "                         \"BIP\": Bounded Inverse Power\n" \
       "                         \"GCM\": Gaussian Core Model\n" \
       "                         \"EXP\": Exponential\n" \
       "  -o, --out:         Output directory for the generated XML file.\n\n"


DIR = "/home/gn/Code/MD-simulation/examples/example_data"
STEPS = 15000
COMPRESSION = "false"
RDF_BINS = 500
PARTICLES_PER_AXIS = 10
LATTICE = "SC"
TRACK_PARTICLES = "false"
RDF_EQUILIBRATE = 2000

RHO = 0.5
TEMPERATURE = 0.5
POWER = 8
A = 0.5
PAIR_POTENTIAL = "BIP"


def create_xml(dir, steps, compression, rdf_bins, particles, lattice, tracking,
               rdf_wait, rho, t, n, a, pp, file_out):

    root = ET.Element("input_variables")
    # Input parameters for the constructor of MD
    input = ET.SubElement(root, "constructor")
    ET.SubElement(input, "output_dir").text = str(dir)
    ET.SubElement(input, "steps").text = str(steps)
    ET.SubElement(input, "compression").text = compression
    ET.SubElement(input, "rdf_bins").text = str(rdf_bins)
    ET.SubElement(input, "particles_per_axis").text = str(particles)
    ET.SubElement(input, "lattice").text = str(lattice)
    ET.SubElement(input, "track_particles").text = tracking
    ET.SubElement(input, "rdf_equilibrate").text = str(rdf_wait)

    # Input parameters for MD::Simulation
    sim = ET.SubElement(root, "simulation_input")

    ET.SubElement(sim, "rho").text = str(rho)
    ET.SubElement(sim, "T").text = str(t)
    ET.SubElement(sim, "n").text = str(n)
    ET.SubElement(sim, "A").text = str(a)
    ET.SubElement(sim, "pp").text = str(pp)

    # Build tree in xml file
    tree = ET.ElementTree(root)
    tree.write(f"{file_out}.xml", pretty_print=True)



def command_line_input(argv):
    try:
        opts, args = getopt.getopt(argv, "hd:s:c:r:p:l:t:w:R:T:n:A:P:o:", [
                                   "dir=", "steps=", "compress=",
                                   "rdf_bins=", "particles=",
                                   "lattice=", "track=", "rdf_wait=",
                                   "rho=", "temp=", "n=", "A=", "pp=",
                                   "out="])
    except getopt.GetoptError:
        print(USAGE)
        print(EXAMPLE)
        sys.exit(2)

    out_dir = None
    steps = None
    compress = None
    rdf_bins = None
    particles = None
    lattice = None
    track_particles = None
    rdf_wait = None
    rho = None
    temperature = None
    n = None
    par_a = None
    pair_potential = None
    output_file = None

    for opt, arg in opts:
        if opt == '-h':
            print(USAGE)
            print(EXAMPLE)
            print(HELP)
            sys.exit()
        elif opt in ("-d", "--dir"):
            out_dir = arg
        elif opt in ("-s", "--steps"):
            steps = arg
        elif opt in ("-c", "--compress"):
            compress = arg
        elif opt in ("-r", "--rdf_bins"):
            rdf_bins = arg
        elif opt in ("-p", "--particles"):
            particles = arg
        elif opt in ("-l", "--lattice"):
            lattice = arg
        elif opt in ("-t", "--track"):
            track_particles = arg
        elif opt in ("-w", "--rdf_wait"):
            rdf_wait = arg
        elif opt in ("-R", "--rho"):
            rho = arg
        elif opt in ("-T", "--temp"):
            temperature = arg
        elif opt in ("-n", "--n"):
            n = arg
        elif opt in ("-A", "--A"):
            par_a = arg
        elif opt in ("-P", "--pp"):
            pair_potential = arg
        elif opt in ("-o", "--out"):
            output_file = arg
        else:
            print("Unrecognised option \"" + opt[0] + "\" entered")

    vars = f"\nInput variables to XML:\n" \
           f" Output file destination: {out_dir}\n" \
           f" Steps: {steps}\n" \
           f" Compression: {compress}\n" \
           f" RDF_bins: {rdf_bins}\n" \
           f" Particles per axis: {particles}\n" \
           f" Lattice selection: {lattice}\n" \
           f" Tracking particles: {track_particles}\n" \
           f" RDF_wait: {rdf_wait}\n" \
           f" -------------------------------------\n" \
           f" RHO: {rho}\n" \
           f" T: {temperature}\n" \
           f" n: {n}\n" \
           f" A: {par_a}\n" \
           f" Pair potential: {pair_potential}\n" \
           f" --------------------------------------\n" \
           f" Output file: {output_file}\n\n"

    print(vars)
    return (out_dir, steps, compress, rdf_bins, particles, lattice, \
            track_particles, rdf_wait, \
            rho, temperature, n, par_a, pair_potential, output_file)


if __name__ == "__main__":
    l = command_line_input(sys.argv[1:])
    print(l)
    if (not all(i is None for i in l)):
        create_xml(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], \
                   l[8], l[9], l[10], l[11], l[12], l[13])
    else:
        print("Error: Uninitialised variables in the arguments list.")