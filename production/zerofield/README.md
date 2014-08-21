VTX Field-on Production for Testing Millepede Alignment
-------------------------------------------------------

## Setup ##

All macros to run the production should be present in the current directory. 

Running the production assumes you have the following libraries built:

1) `libProdEventCutter.so` - build from `vtx-align/ProdEventCutter/` (git)

2) `libAnaVTXCluster.so` - build from `vtx-align/AnaVTXCluster/` (git)

3) `libTestVTXProduction.so` - build from `offline/AnalysisTrain/mcglinchey_vtxprodcheck/` (cvs)


## Running a production ##

The condor job file necessary for running the production can be produced using the `makeCondorSubmit.py` module. It takes a series of command line arguments to setup the job file. Use `makeCondorSubmit.py -h` for an explanation of expected arguments.

Each condor job is controlled by `run_alignment_production.py`, which does all the necessary heavy lifting, like setting up the production environemnt, copying the PRDFF from dcache and running the production & analysis macros.

## Notes ##

`makeCondorSubmit.py` will create the output directory specified with the `-o` option if it doesn't already exist. It will also create the `condor_logs/` directory expected by condor. 

The `standalone` flag in `run_alignment_production.py` to set which production macro is called
 - `standalone = False` grabs `Fun4All_VTX_ZeroField.C` which runs the full production (standalone + svxcentral).
 - `standalone = True` grabs `Fun4All_VTX_ZeroField_Standalone.C` which runs the standalone tracking only.
