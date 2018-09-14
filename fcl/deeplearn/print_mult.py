from ROOT import TChain
import sys
ch=TChain("partroi_segment_tree")
ch.AddFile(sys.argv[1])

def mean(ar):
    ctr=0
    for v in ar: ctr+=v
    return float(ctr)/len(ar)

mult={}

for entry in xrange(ch.GetEntries()):

    ch.GetEntry(entry)
    roi_v = ch.partroi_segment_branch.ROIArray()
    if not roi_v.size() in mult:
        mult[roi_v.size()] = [[],[]]
    sctr=0
    tctr=0
    for roi in roi_v:
        if roi.PdgCode() in [11,22]: sctr+=1
        else: tctr+=1
    mult[roi_v.size()][0].append(sctr)
    mult[roi_v.size()][1].append(tctr)

keys = mult.keys()
keys.sort()
for key in keys:
    print key, len(mult[key][0]),mean(mult[key][0]),mean(mult[key][1])



