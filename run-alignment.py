import sys
import os

on = True # enable/disable execution
run = 411768
p = 6  # production step
itersteps = [20, 21, 22]

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

    if i != itersteps[-1]:
        
        # Run VtxAlign.C
        m = "'VtxAlign.C+({},{},{})' >& logs/{}-{}-{}.log".format(run, p, i,
                                                                  run, p, i)
        print('Executing ' + root + m)
        if on:
            os.system(root + m)

        # Run DrawResults.C
        m = "'DrawResults.C+({},{},{},{},{})'".format(run, p, i, p, i + 1)
        print('Executing ' + root + m)
        if on:
            os.system(root + m)
