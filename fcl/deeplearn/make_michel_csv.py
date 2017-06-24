import sys
from ROOT import larutil
padding=110

larp = larutil.LArProperties.GetME()
data = {}
for l in open(sys.argv[1],'r').read().split('\n'):
    words = l.split()
    if not len(words) == 5: continue
    run    = int(words[0])
    subrun = int(words[1])
    event  = int(words[2])
    
    z = int(words[3])
    x = int(words[4])
    if x/6 < padding or x > (6400 - 6*padding):
        continue
    if z < padding or z > (3456 - padding):
        continue
    
    z = float(z) * 0.3
    x = (float(x) - 800) / 2. * larp.DriftVelocity()
    event_id = (run,subrun,event)
    if not event_id in data.keys():
        data[event_id] = (x,0.,z)
        continue

    print 'Duplicate event:',event_id
    print 'Previous:',data[event_id]
    print 'Now:',(x,0.,z)
    dist = abs(128. - x)
    if z < 100 : dist += z
    if z > 900 : dist += 1030 - z
    pre_dist = abs(128. - data[event_id][0])
    if data[event_id][2] < 100: pre_dist += data[event_id][2]
    if data[event_id][2] > 900: pre_dist += 1030. - data[event_id][2]
    if dist < pre_dist:
        print '\033[93mUpdating\033[00m'
        data[event_id] = (x,0.,z)
    print
fout=open('michel_constraint.csv','w')
fout.write('run,subrun,event,x,y,z\n')
keys = data.keys()
keys.sort()
for key in keys:
    run,subrun,event = key
    x = data[key][0]
    z = data[key][2]
    fout.write('%d,%d,%d,%g,0.0,%g\n' % (run,subrun,event,x,z))
fout.close()
