#!/opt/phenix/bin/python
'''
Makes condor submit file for zf alignment production
'''

import os
import sys
from argparse import ArgumentParser

#First get any command line arguments
parser = ArgumentParser(description = """Make condor submit script
										 for zerofield alignment production""")

parser.add_argument("-f", "--file", dest="filename",
                  help="write output to FILE", metavar="FILE")
parser.add_argument("-c", "--config", dest="config",
					help="configuration file", metavar="CONFIG")
parser.add_argument("-r", "--run", dest="run", help="run number")
parser.add_argument("-n", "--nsegments", dest="nsegments", help="number of segments")
parser.add_argument("-o", "--output", dest="outdir",
					help="Write output to DIR", metavar="DIR")
parser.add_argument("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="print status messages to stdout")

parser.set_defaults(filename="condor.job",
					config="/direct/phenix+u/dcm07e/work/vtx-align/production/config/config-ideal.txt",
					run='406541',
					nsegments='24',
					outdir="/direct/phenix+prod01/phnxreco/millepede/test/")

args = parser.parse_args()



#Check the input arguments for compatabiity and fix if necessary
if not args.outdir.endswith('/'):
	args.outdir += '/'

if args.config.startswith('.'):
	print("\nWARNING!! Configuration file requires an absolute path!\n")
	print(" Fixing now based on your current working directory.")
	cwd = os.getcwd()
	args.config = cwd[0:cwd.rfind('/')]  + args.config[args.config.find('/'):]
	print(" Config file now: {}\n".format(args.config))

if not os.path.isfile(args.config):
	print("\nERROR!! Can't find config file - {}\n".format(args.config))
	sys.exit(1)


#Print out arguments
if args.verbose :
	print(" Writting job file to {}".format(args.filename))
	print(" Using config file {}".format(args.config))
	print(" Using run {}".format(args.run))
	print(" Producing {} segments".format(args.nsegments))
	print(" Output of condor jobs will go to {}".format(args.outdir))

#Check whether the output directory specified exists
# if not, create it
os.system("mkdir -p {}condor_logs".format(args.outdir))


#make the file
fo = open(args.filename,'w')

#write the job information
fo.write("Universe        = vanilla\n")
fo.write("Executable      = {}/run_alignment_production.py\n".format(os.getcwd()))
fo.write("Arguments       = {} {} $(Process) {} {}\n".format(args.config,args.run,args.outdir,os.getcwd()[0:-20]))
fo.write("Requirements    = (CPU_Speed >= 1 && CPU_Experiment == \"phenix\")\n")
fo.write("Rank            = CPU_Speed\n")
fo.write("Priority        = +1\n")
fo.write("GetEnv          = True\n")
fo.write("Initialdir      = {}\n".format(os.getcwd()))
fo.write("Output          = {}condor_logs/ran_production_{:0>10}-$(Process).out\n".format(args.outdir,args.run,))
fo.write("Error           = {}condor_logs/ran_production_{:0>10}-$(Process).err\n".format(args.outdir,args.run,))
fo.write("Log             = {}condor_logs/ran_production_{:0>10}-$(Process).log\n".format(args.outdir,args.run,))
fo.write("+Experiment     = \"phenix\"\n")
fo.write("PeriodicHold = (NumJobStarts >= 1 && JobStatus == 1)\n")
fo.write("+Job_Type       = \"cas\"\n")
fo.write("Queue {}".format(args.nsegments))
