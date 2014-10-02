Procedure for running millipede VTX alignment
=============================================


1) Run a production with the ideal geometry and zeroes for all offsets.

2) Run a production with the ideal geometry and zero offsets for vtx-to-cnt and east-to-west. Include a beamcenter determined from step 1. Use this beamcenter for all further alignments & productions.

3) Run millipede on the output of #2 in halflayer mode with [x,y,s,z] parameters. Use two internal iterations.

4) Run production using output geometry from #3. Use the same beam center as in #2. Again, 0's for all offsets.

5) Repeat #3-4

6) Run millipede on production output from #5 in ladder mode with [x,y,z] parameters.

7) Run production using output geometry from #6. Use the same beam center as in #2. Again, 0's for all offsets.