from ROOT import gSystem, gROOT
import os

geo = "/Users/adare/phenix/svxgeo"
geolib = ''.join((os.environ['TMPDIR'], geo, '/SvxTGeo_C.so'))
gSystem.Load(geolib)
from ROOT import SvxTGeo


# infile  = 'geom/svxPISA-ideal-sdg.par'
# outfile = 'geom/svxPISA-ideal.par'
infile  = 'geom/svxPISA-411768-sdg.par'
outfile = 'geom/svxPISA-411768.par'

tgeo = SvxTGeo()
tgeo.ReadParFile(infile)
tgeo.MakeTopVolume(100, 100, 100)
tgeo.AddSensors()
tgeo.GeoManager().CloseGeometry()
tgeo.WriteParFile(outfile)

# Numerically diff text files using ndiff:
# http://www.math.utah.edu/~beebe/software/ndiff/
cmd = 'ndiff -abserr 1.0e-4 {} {}'.format(infile, outfile)

print('Diffing output:')
os.system(cmd)
