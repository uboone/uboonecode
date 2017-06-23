from ROOT import TChain
import sys

producer=sys.argv[2]

ch=TChain("opdigit_%s_tree" % producer)
ch.AddFile(sys.argv[1])
print ch.GetEntries()

ch.GetEntry(0)

wf_v = None
exec('wf_v = ch.opdigit_%s_branch' % producer)

for entry in xrange(ch.GetEntries()):
    ch.GetEntry(entry)
    wf_info_m={}
    for wf in wf_v:
        channel = wf.ChannelNumber()
        if not channel in wf_info_m:
            wf_info_m[channel]=[]
        wf_info_m[channel].append(wf.size())

    keys=wf_info_m.keys()
    keys.sort()

    print 'entry',entry
    for k in keys:
        print k,wf_info_m[k]
    print
