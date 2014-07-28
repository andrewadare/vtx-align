VTX Field-on Production for Testing Millepede Alignment
-------------------------------------------------------

## Setup ##

All macros to run the production should be present in the current directory. For running the analysis module on the produced DST's, you need a compiled copy of the `TestVTXProduction` macro found checked in to CVS in `offline/AnalysisTrain/mcglinchey_vtxprodcheck/`

## Running a production ##

The condor job file necessary for running the production can be produced using the `makeCondorSubmit.py` module. It takes a series of command line arguments to setup the job file. Use `python makeCondorSubmit.py -h` for an explanation of expected arguments.

Each condor job is controlled by `run_alignment_production_fieldon.py`, which does all the necessary heavy lifting, like setting up the production environemnt, copying the PRDFF from dcache and running the production & analysis macros.

Before submitting jobs, the variable `vtxalignDir` (line 24) should be set to a path pointing to your checked out vtx-align directory.

## Notes ##