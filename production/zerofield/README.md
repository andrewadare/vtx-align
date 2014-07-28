VTX Field-on Production for Testing Millepede Alignment
-------------------------------------------------------

## Setup ##

All macros to run the production should be present in the current directory. For running the analysis module on the produced DST's, you need a compiled copy of the `AnaVTXCluster` module found checked in this git repository.

## Running a production ##

The condor job file necessary for running the production can be produced using the `makeCondorSubmit.py` module. It takes a series of command line arguments to setup the job file. Use `makeCondorSubmit.py -h` for an explanation of expected arguments.

Each condor job is controlled by `run_alignment_production.py`, which does all the necessary heavy lifting, like setting up the production environemnt, copying the PRDFF from dcache and running the production & analysis macros.

Before submitting jobs, the variable `vtxalignDir` (line 24) should be set to a path pointing to your checked out vtx-align directory.

## Notes ##

`makeCondorSubmit.py` will create the output directory specified with the `-o` option if it doesn't already exist. It will also create the `condor_logs/` directory expected by condor. 

When using the `-c` option to specify the configuration file, you must give it the full path.