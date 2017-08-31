import sys, os

runids={}
for l in open(sys.argv[2],'r').read().split('\n'):
    words=l.split()
    if not len(words)==3 or not words[0].isdigit(): continue
    runids[((int(words[0]),int(words[1]),int(words[2])))] = False


fout=open('ccpi0_constraint.csv.updated','w')
for l in open(sys.argv[1],'r').read().split('\n'):
    words=l.split(',')
    if not len(words)==6 or not words[0].isdigit(): continue
    runid = (int(words[0]),int(words[1]),int(words[2]))
    if not runid in runids:
        fout.write('%s\n' % l)
    else:
        runids[runid] = True

fout.close()
for runid,used in runids.iteritems():
    if used: continue
    print 'Unused ID:',runid


    
