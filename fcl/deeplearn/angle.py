import sys
from ROOT import TChain
import numpy as np
ch=TChain("mctruth_generator_tree")
ch.AddFile(sys.argv[1])
fout=open('ana.csv','w')
fout.write('entry,idx,pdg,dcosx,dcosy,dcosz\n')

num_entries = ch.GetEntries()
if len(sys.argv)>2:
    num_entries = int(sys.argv[2])
    if num_entries > ch.GetEntries():
        num_entries = ch.GetEntries()

for entry in xrange(num_entries):

    ch.GetEntry(entry)
    mct = ch.mctruth_generator_branch
    part_v = mct.GetParticles()
    for idx in xrange(part_v.size()):
        part = part_v[idx]
        if not part.StatusCode() == 1: continue
        pdg = int(np.abs(part.PdgCode()))
        if not pdg in [11,13,22,2212,211]: continue

        px = part.Momentum().X()
        py = part.Momentum().Y()
        pz = part.Momentum().Z()

        mom_mag = np.sqrt(np.power(px,2)+np.power(py,2)+np.power(pz,2))
        dcosx = px / mom_mag
        dcosy = py / mom_mag
        dcosz = pz / mom_mag
        
        fout.write('%d,%d,%d,%g,%g,%g\n' % (entry,idx,pdg,dcosx,dcosy,dcosz))
fout.close()

