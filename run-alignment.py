import sys, os

narg = len(sys.argv)
run = 411768
iter = 0

if narg==1:
  pass # Use defaults above
elif narg==2 and sys.argv[1].isdigit():
  iter = sys.argv[1]
elif narg==3 and sys.argv[1].isdigit() and sys.argv[2].isdigit():
  run = sys.argv[1]
  iter = sys.argv[2]
else:
  print("Usage: python run-alignment.py <run> <iter> "),
  print("(run number is optional and iter = 0,1,2,...)")
  sys.exit(0)

root = 'time root -b -q'
cmd1 = "{} 'VtxAlign.C+({},{})' >& logs/iter{}.log".format(root,run,iter,iter)
cmd2 = "{} 'DrawResults.C({},{})'".format(root,run,iter)

print('Executing ' + cmd1)
os.system(cmd1)
print('Executing ' + cmd2)
os.system(cmd2)
