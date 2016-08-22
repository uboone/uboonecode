from ROOT import TChain
import sys
ch0=TChain("opdigit_pmtreadout_tree")
ch1=TChain("trigger_triggersim_tree")

ch0.AddFile(sys.argv[1])
ch1.AddFile(sys.argv[1])

ch0.GetEntry(0)
ch1.GetEntry(0)

br = ch0.opdigit_pmtreadout_branch
tr = ch1.trigger_triggersim_branch

ctr=0
ttime=[]

for x in br:
    if x.size()<1000: continue
    if not x.TimeStamp() in ttime: ttime.append(x.TimeStamp())
    ctr+=1
    
print tr.TriggerTime()
print ttime
print ctr
