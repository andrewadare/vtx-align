TODO items for VTX alignment
----------------------------


# Add beamcenter and vertex to tree #
** 9/29/2014 **

Add fields for the beam center and primary vertex to output tree (see `CreateTree()` in `VtxIO.h`) so that the calculations can be carried around rather than re-calculated each time.

# Possibly create SvxGeoEvent class #
** 9/29/2014 **

Create an `SvxGeoEvent` class with member variables:

1) array of `SvxGeoTrack` objects
2) beamcenter (x,y)
3) primary vertex (x,y,z)

This would allow us to carry around the calculated parameters. 

*Needs further thought on whether this is truly desireable and worth the time*



