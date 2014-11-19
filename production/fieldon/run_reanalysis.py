#!/opt/phenix/bin/python
'''
Python script for re-analyzing the DST's produced
from a previous production. Overwrites previous
ana output.
'''

import os
import sys
import stat

##############################################
# Get the command line arguments
##############################################
runNumber = int(sys.argv[1])
segNumber = int(sys.argv[2])
outDir = sys.argv[3]
vtxalignDir = sys.argv[4]

##############################################
# Set some parameters for running
##############################################
nevents = 0
productionDir = vtxalignDir + "production/fieldon/"



##############################################
# Set environment variables
##############################################
print("\n--> Setting environment variables")

#os.environ["ODBCINI"] = "/opt/phenix/etc/odbc.ini.master"
#print("ODBCINI: {}".format(os.environ["ODBCINI"]))


# get the name of the dst
filenum = "DST_SVX_MB-{:0>10}-{:0>4}.root".format(runNumber,segNumber)
outfiles = os.listdir("{}DST_SVX/".format(outDir))
for file in outfiles:
	if str(filenum) in file:
		dstFile = "{}DST_SVX/{}".format(outDir,file)
		break
print("Found DST: {}".format(dstFile))

##############################################                                                                                                                                                                                                             
# run over the data, producing the ntuples                                                                                                                                                                                                                 
##############################################                                                                                                                                                                                                             
print("\n--> Running over the DST to produce ntuples")

os.chdir(outDir)
os.system("mkdir -p ntuple")
os.chdir(productionDir)
print("current directory: " + os.getcwd())

ntupleFile = "{}ntuple/testvtxproduction_{:0>10}-{:0>4}.root".format(outDir,runNumber,segNumber)
print("Output Ntuple: {}".format(ntupleFile))

command = "root -b -q \'analyze_run14_test_prod_output.C("
command += str(nevents) + ","
command += "\"" + dstFile + "\","
command += "\"" + ntupleFile
command += ")\'"
print(command)

os.system(command)

##############################################
# cleanup
##############################################
print("\n--> Cleaning up.")

#make all the directories we were in group writable
os.system("chmod a+w {}".format(outDir))
for root, dirs, files in os.walk(outDir, topdown=False):
        for dir in dirs:
            os.system("chmod a+w {}".format(outDir+dir))

##############################################
# DONE!
##############################################

print("\n---- DONE! ----\n")




