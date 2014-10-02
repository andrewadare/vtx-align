TODO items for VTX alignment
----------------------------

# Modify FilterData() in DataRejector.h to use input beam center #
*10/2/2014*

Now that we are trying to fix the beamcenter, modify to use input beamcenter rather than calculating it from the tracks.

# Modify CalcBeamCenter.C to use input beamcenter #
*10/2/2014*

# Add beamcenter and vertex to tree #
*9/29/2014*

Add fields for the beam center and primary vertex to output tree (see `CreateTree()` in `VtxIO.h`) so that the calculations can be carried around rather than re-calculated each time.




