'''
Python script for running a small zero field production
for use in the millipede alignment
'''

import os
import sys

##############################################
# Get the command line arguments
##############################################
configFileName = sys.argv[1]
runNumber = int(sys.argv[2])
segNumber = int(sys.argv[3])
outDir = sys.argv[4]

##############################################
# Set some parameters for running
##############################################
nevents = 1000
#condorDir = os.environ["_CONDOR_SCRATCH_DIR"]
condorDir = "/direct/phenix+prod01/phnxreco/millepede/test/"
#outDir = "/phenix/prod01/phnxreco/millepede/ZF_production/"
vtxalignDir = "/direct/phenix+u/dcm07e/work/vtx-align/production/zerofield/"


##############################################
# Get the configuration parameters from file
##############################################
print("\n--> Reading in parameters from {}".format(configFileName))
with open(configFileName,'r') as file:
	for line in file:
		if (line[0] == '#' and line[0] == '\r' and line[0] == '\n'):
			continue
		if "beamcenter" in line:
			beamcenter = line[:-2].split(' ')[1:]
			print("beamcenter: ({}, {})".format(beamcenter[0],beamcenter[1]))
		if "east-to-west" in line:
			eastToWest = line[:-2].split(' ')[1:]
			print("east-to-west offset: ({}, {}, {})".format(eastToWest[0], eastToWest[1], eastToWest[2]))
		if "vtx-to-cnt" in line:
			vtxToCNT = line[:-2].split(' ')[1:]
			print("VTX-to-CNT offset: ({}, {}, {})".format(vtxToCNT[0], vtxToCNT[1], vtxToCNT[2]))
		if "geomfile" in line:
			parfile = line[:-2].split(' ')[1]
			print("parfile: {}".format(parfile))



##############################################
# Set environment variables
##############################################
print("\n--> Setting environment variables")

os.environ["ODBCINI"] = "/opt/phenix/etc/odbc.ini.master"
os.environ["DCACHE_DOOR"] = "phnxdoor1.rcf.bnl.gov:22133"

print("ODBCINI: {}".format(os.environ["ODBCINI"]))
print("DCACHE_DORR: {}".format(os.environ["DCACHE_DOOR"]))

##############################################
# Get the desired PRDFF
##############################################
print("\n--> Getting the desired PRDFF")
os.chdir(condorDir)
print(os.getcwd())

runRangeLow = runNumber - runNumber%1000
runRangeHigh = runRangeLow + 1000
prdfDir = "/pnfs/rcf.bnl.gov/phenix/phnxsink/run14/eventdata/run_{:0>10}_{:0>10}/".format(runRangeLow,runRangeHigh)
prdfFile = "EVENTDATA_P00-{:0>10}-{:0>4}.PRDFF".format(runNumber,segNumber)

print(prdfDir+prdfFile)
print("size: {}".format(os.stat(prdfDir+prdfFile).st_size))
os.system("dccp {} .".format(prdfDir+prdfFile))
os.system("ls -lh {}".format(prdfFile))

##############################################
# Run the production
##############################################
print("\n--> Running the production")

#setup the production directory
os.chdir(condorDir)
print("production directory: " + os.getcwd())

os.system("LuxorLinker.pl -1 {}".format(runNumber))
os.system("ln -sf {}Fun4All_VTX_ZeroField.C .".format(vtxalignDir))
os.system("ln -sf {}TrigSelect.C .".format(vtxalignDir))
os.system("ln -sf {}OutputManager.C .".format(vtxalignDir))
os.system("ln -sf {}rawdatacheck.C .".format(vtxalignDir))

# set up the parameters to be passed to the macro
pixel_refmap = "/direct/phenix+hhj2/dcm07e/run14MiniProd/fieldon/blank_pixel_refmap.txt"
pixel_diffmap = "/direct/phenix+hhj2/dcm07e/vtx/deadmaps_Run14AuAu200/singleruns/pixel_deadmap_run14auau200_run{}.dat".format(runNumber)
pixel_chipmap = "/direct/phenix+hhj2/dcm07e/vtx/deadmaps_Run14AuAu200/singleruns/chip_deadmap_run14auau200_run{}.dat".format(runNumber)
strip_deadchannel = "/direct/phenix+hhj/theok/theok/stability/strip_deadmaps/run14/run{}_strip_hotdeadChannels.txt".format(runNumber)
strip_deadRCC = "/direct/phenix+hhj/theok/theok/stability/strip_deadmaps/run14/run{}_strip_hotdeadReadouts.txt".format(runNumber)

command = "root -b -q \'Fun4All_VTX_ZeroField.C("
command += str(nevents) +","
command += "\"" + prdfFile + "\","
command += "\"" + outDir + "DST_SVX/" + "\","
command += "\"" + vtxalignDir + "geom/" + parfile + "\","
command += str(beamcenter[0]) + "," + str(beamcenter[1]) + ","
command += str(eastToWest[0]) + "," + str(eastToWest[1]) + "," + str(eastToWest[2]) + ","
command += str(vtxToCNT[0]) + "," + str(vtxToCNT[1]) + "," + str(vtxToCNT[2]) + ","
command += "\"" + pixel_refmap + "\","
command += "\"" + pixel_diffmap + "\","
command += "\"" + pixel_chipmap + "\","
command += "\"" + strip_deadchannel + "\","
command += "\"" + strip_deadRCC
command += ")\'"

print(command)

os.system(command)

# get the name of the dst
outfiles = os.listdir("{}DST_SVX/".format(outDir))
for file in outfiles:
	if str(runNumber) and str(segNumber) in file:
		dstFile = "{}DST_SVX/{}".format(outDir,file)
		break
print("Found DST: {}".format(dstFile))

##############################################
# run over the data, producing the ntuples
##############################################
print("\n--> Running over the DST to produce cluster ntuple")

os.chdir(outDir)
os.system("mkdir -p ntuple")
os.chdir(condorDir)
print("current directory: " + os.getcwd())

ntupleFile = "{}ntuple/anavtxcluster_{:0>10}-{:0>4}.root".format(outDir,runNumber,segNumber)
print("Output Ntuple: {}".format(ntupleFile))

command = "root -b -q \'run_anavtxcluster.C("
command += str(nevents) + ","
command += "\"" + dstFile + "\","
command += "\"" + ntupleFile
command += ")\'"
print(command)

os.system(command)

##############################################
# DONE!
##############################################

print("\n---- DONE! ----\n")




