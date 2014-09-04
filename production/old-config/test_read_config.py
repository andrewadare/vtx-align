'''
Test the reading of alignment configuration file
'''

from argparse import ArgumentParser

parser = ArgumentParser(description = "Test the reading of config files")

parser.add_argument("-f", "--file", dest="filename",
                  help="read config file FILE", metavar="FILE")

parser.set_defaults(filename="config-example.txt")

args = parser.parse_args()


##############################################
# Get the configuration parameters from file
##############################################
print("\n--> Reading in parameters from {}".format(args.filename))
with open(args.filename,'r') as file:
	for line in file:
		if line[0] == '#':
			continue
		if "beamcenter" in line:
			beamcenter = line.rstrip().split(' ')[1:]
			print("beamcenter: ({}, {})".format(beamcenter[0],beamcenter[1]))
		if "east-to-west" in line:
			eastToWest = line.rstrip().split(' ')[1:]
			print("east-to-west offset: ({}, {}, {})".format(eastToWest[0], eastToWest[1], eastToWest[2]))
		if "vtx-to-cnt" in line:
			vtxToCNT = line.rstrip().split(' ')[1:]
			print("VTX-to-CNT offset: ({}, {}, {})".format(vtxToCNT[0], vtxToCNT[1], vtxToCNT[2]))
		if "geomfile" in line:
			parfile = line.rstrip().split(' ')[1]
			print("parfile: {}".format(parfile))