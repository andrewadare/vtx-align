import sys
import os

on = True # enable/disable execution
log = False # Redirect stdout from alignment to logfile
run = 123456
p = 0  # production step
itersteps = [0]

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
    if log:
        m = "'VtxAlign.C+({},{},{})' >& logs/{}-{}-{}.log".format(run, p, i,
                                                                  run, p, i)
    else:
        m = "'VtxAlign.C+({},{},{})'".format(run, p, i)
    print('Executing ' + root + m)
    if on:
        os.system(root + m)
