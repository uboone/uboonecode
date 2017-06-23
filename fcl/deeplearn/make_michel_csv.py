import sys
from ROOT import larutil

larp = larutil.LArProperties.GetME()
fout=open('michel_constraint.csv','w')
fout.write('run,subrun,event,x,y,z\n')

for l in open(sys.argv[1],'r').read().split('\n'):
    words = l.split()
    if not len(words) == 5: continue
    run    = int(words[0])
    subrun = int(words[1])
    event  = int(words[2])
    z = float(words[3]) * 0.3
    x = (float(words[4]) - 800) / 2. * larp.DriftVelocity()
    fout.write('%d,%d,%d,%g,0.0,%g\n' % (run,subrun,event,x,z))
fout.close()
