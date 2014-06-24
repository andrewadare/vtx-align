% README file for vtx-align repository

# General info
This repository contains the code used to generate input for Millepede II, run the alignment fitter, visualize the results, and get new PISA geometry parameters in text format. 

The version of Millepede II used for VTX alignment is different than the original version wrapped for PHENIX (for that code see `/offline/packages/mutoo/alignment`). This version (currently v4-02-00) consists of (1) the `Mille` `C++` class (included in `mille/`) and (2) the `pede` executable built from FORTRAN 90 source. Through calls to `Mille` member functions, one or more binary files containing parameter and residual information are generated for consumption by `pede`. One or more text files containing constraint information are also passed to `pede`. A steering file lists the files to be used for input.

# Dependencies
This code requires the `svxgeo` library from `offline/packages/svxgeo`. The rootlogon.C script provides the option to either load the version installed at RCF via `gSystem->Load("libsvxgeo")` or to build the needed libraries from source (in which case the source path is needed).

The `DrawResults.C` script requires `UtilFns.h` which is available from `https://github.com/andrewadare/utils.git`.

# Geometry files in `geom/`
The initial par file is from `svxDetectorGeo::Write_svxPISApar()` as run by the macro `GetParFile.C` in `offline/packages/svx/wrk/`.

The `reformat.py` script uses `SvxTGeo` to read in and write a tidier version that is the same (it is numerically diffed to be sure.)

