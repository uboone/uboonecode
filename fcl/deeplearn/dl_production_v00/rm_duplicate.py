import shutil,commands,os

DIR='/pnfs/uboone/scratch/users/vgenty/wirecalib_v00_p00_log/v05_08_00/'

for job in xrange(389):
    dirs=['%s/%s' % (DIR,f) for f in commands.getoutput('ls %s | grep "_%d"' % (DIR,job)).split() if f.endswith('_%d' % job)]
    if len(dirs) <= 1: continue
    dirm={}
    for d in dirs:
        print os.path.getmtime(d),d
        dirm[os.path.getmtime(d)] = d
    keys=dirm.keys()
    keys.sort()
    for x in xrange(len(keys)-1):
        key = keys[x]
        print 'removing',dirm[key]
        shutil.rmtree(dirm[key])
    print job

