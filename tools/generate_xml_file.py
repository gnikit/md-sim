import os
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

DIR = "/home/gn/Code/MD-simulation/examples/example_data"
STEPS = 10000
COMPRESSION = "false"
RHO = 0.5
TEMPERATURE = 0.5
POWER = 6
A = 0.5

root = ET.Element("input_variables")

# Input paramters for the constructor of MD
input = ET.SubElement(root, "constructor")

ET.SubElement(input, "output_dir").text = str(DIR)
ET.SubElement(input, "steps").text = str(STEPS)
ET.SubElement(input, "compression").text = COMPRESSION

# Input parameters for MD::Simulation
sim = ET.SubElement(root, "simulation_input")

ET.SubElement(sim, "rho").text = str(RHO)
ET.SubElement(sim, "T").text = str(TEMPERATURE)
ET.SubElement(sim, "n").text = str(POWER)
ET.SubElement(sim, "A").text = str(A)

# Build tree in xml file
tree = ET.ElementTree(root)
tree.write("schemas/input_schema.xml")
