import os

files = [x for x in os.listdir('.') if x.startswith('lite') and x.endswith('.fcl')]

for f in files:
    print f

    contents = open(f,'r').read().split('\n')

    fout = open(f,'w')
    for line in contents:
        if line.find('PROLOG') >=0: continue
        fout.write('%s\n' % line)
    fout.close()
