TODO items for VTX alignment
----------------------------
# Parse the "sim" option in ConstraintBuilder.h to keep all ladders.
*10/7/2014*
We do not want to exclude any ladders in simulations. Currently I hack this by editing BadLadders.h directly (commenting the ladder list or returning immediately).

Now VtxAlign takes the option "sim" as part of alignMode.

I need to pass the option "sim" through to LadderLabels(), or do something with equivalent effect. LadderLabels() should skip checking the Dead() function from BadLadders.h when using simulated data.

# Modify FilterData() in DataRejector.h to use input beam center #
*10/2/2014*

Now that we are trying to fix the beamcenter, modify to use input beamcenter rather than calculating it from the tracks.

# Modify CalcBeamCenter.C to use input beamcenter #
*10/2/2014*

# Add beamcenter and vertex to tree #
*9/29/2014*

Add fields for the beam center and primary vertex to output tree (see `CreateTree()` in `VtxIO.h`) so that the calculations can be carried around rather than re-calculated each time.




