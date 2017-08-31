import sys
from ROOT import larutil
padding=110

data = {}
for l in open(sys.argv[1],'r').read().split('\n'):
    words = l.split()
    if not len(words) == 6: continue
    run    = int(words[0])
    subrun = int(words[1])
    event  = int(words[2])
    
    x = float(words[3])
    y = float(words[4])
    z = float(words[5])

    event_id = (run,subrun,event)
    if not event_id in data.keys():
        data[event_id] = (x,y,z)
        continue

    print 'Duplicate event:',event_id
    print 'Previous:',data[event_id]
    print 'Now:',(x,y,z)
    dist = abs(128. - x) + abs(y)
    if z < 100 : dist += z
    if z > 900 : dist += (1030 - z)
    pre_dist = abs(128. - data[event_id][0]) + abs(data[event_id][1])
    if data[event_id][2] < 100: pre_dist += data[event_id][2]
    if data[event_id][2] > 900: pre_dist += 1030. - data[event_id][2]
    if dist < pre_dist:
        print '\033[93mUpdating\033[00m'
        data[event_id] = (x,y,z)
    print
fout=open('pi0_constraint.csv','w')
fout.write('run,subrun,event,x,y,z\n')
keys = data.keys()
keys.sort()
for key in keys:
    run,subrun,event = key
    x = data[key][0]
    z = data[key][2]
    fout.write('%d,%d,%d,%g,0.0,%g\n' % (run,subrun,event,x,z))
fout.close()
