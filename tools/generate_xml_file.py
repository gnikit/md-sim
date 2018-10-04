import sys
import xml.etree.ElementTree as ET

# DIR = sys.argv[1]
# STEPS = sys.argv[2]
# RHO = sys.argv[3]
# TEMPERATURE = sys.argv[4]
# POWER = sys.argv[5]
# A = sys.argv[6]
# TODO: add exception handling 
# TODO: add int __main__ for file execution from the terminal

DIR = "/home/gn/Code/MD-simulation/examples/"
STEPS = 10000
RHO = 0.5
TEMPERATURE = 0.5
POWER = 6
A = 0.5

root = ET.Element("input_variables")

# Input paramters for the constructor of MD
input = ET.SubElement(root, "constructor")

ET.SubElement(input, "dir", name="ouput_dir").text = str(DIR)
ET.SubElement(input, "steps", name="steps_N").text = str(STEPS)

# Input parameters for MD::Simulation
sim = ET.SubElement(root, "simulation_input")

ET.SubElement(sim, "rho", name="density_rho").text = str(RHO)
ET.SubElement(sim, "T", name="temperature_T").text = str(TEMPERATURE)
ET.SubElement(sim, "n", name="potential_strenght_n").text = str(POWER)
ET.SubElement(sim, "A", name="softening_parameter_a").text = str(A)

# Build tree in xml file
tree = ET.ElementTree(root)
tree.write("tools/input_schema.xml")
