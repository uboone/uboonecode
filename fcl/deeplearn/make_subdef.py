import sys,os

parent_def = sys.argv[1]
new_def    = sys.argv[2]
csv_file   = sys.argv[3]

runs=[]
for l in open(csv_file,'r').read().split('\n'):
    words = l.split(',')
    if not words[0].isdigit(): continue
    run = int(words[0])
    subrun = int(words[1])
    runid = (run,subrun)
    if runid in runs: continue
    runs.append(runid)

samdef = 'defname: %s and run_number ' % parent_def
for idx,runid in enumerate(runs):
    if idx: samdef += ', %d.%d' % runid
    else: samdef += ' %d.%d' % runid

samdef = 'samweb create-definition %s "%s"' % (new_def,samdef)
print
print samdef
print
os.system(samdef)



