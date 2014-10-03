import sys
import os

on = True # enable/disable execution
run = 411768
p = 0  # production step
itersteps = [0,1]

pat = '{}-{}-{}'.format(run, p, 0)
g = 'geom/{}.par'.format(pat)
r = 'rootfiles/{}.root'.format(pat)
root = 'time root -b -q '

if os.path.isfile(g) == False:
    print('Error: {} not found.'.format(g))
    sys.exit()

if os.path.isfile(r) == False:
    print('Error: {} not found.'.format(g))
    sys.exit()

for i in itersteps:

    # Run CalcBeamCenter.C
    m = "'CalcBeamCenter.C+({},{},{})'".format(run, p, i)
    print('Executing ' + root + m)
    if on:
        os.system(root + m)
        os.system("open pdfs/beam-center-run{}-pro{}-sub{}.pdf".format(run,p,i))

    # Run DrawResults.C
    m = "'DiffGeometry.C(\"geom/svxPISA-ideal.par\", \"geom/{}-{}-{}.par\")'".format(run, p, i)
    print('Executing ' + root + m)
    if on:
        os.system(root + m)
        os.system("open pdfs/svxPISA-ideal-vs-{}-{}-{}.pdf".format(run,p,i))

    if i != itersteps[-1]:
        
        # Run VtxAlign.C
        m = "'VtxAlign.C+({},{},{})' >& logs/{}-{}-{}.log".format(run, p, i,
                                                                  run, p, i)
        print('Executing ' + root + m)
        if on:
            os.system(root + m)
            os.system("open pdfs/run{}-pro{}sub{}-vs-pro{}sub{}-vtxtrks.pdf".\
                      format(run,p,i,p,i+1))

