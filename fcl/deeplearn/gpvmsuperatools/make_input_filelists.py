import os,sys
import ROOT as rt

project = "mcc7bnbcosmic"

if project=="mcc7cosmic":
    original_procs_file = "orig_mcc7cosmic_procs.txt"
    out_procfile = "procs.txt"
    tmp_inputlist_dir = "jobfilelists"
    larlite_simlink_dir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/simlinks/mcc7_cosmic_v00_p00/"
    supera_outdir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/supera/mcc7_cosmic_v00_p00"
    ismc = True
elif project=="mcc7bnbcosmic":
    original_procs_file = "orig_procs_mcc7bnbcosmics.txt"
    out_procfile = "procs_mcc7_bnb_cosmic_v00_p00.txt"
    tmp_inputlist_dir = "jobfilelists"
    larlite_simlink_dir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/simlinks/mcc7_bnb_cosmic_v00_p00/"
    supera_outdir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/supera/mcc7_bnb_cosmic_v00_p00"
    ismc = True
elif project=="extbnb":
    original_procs_file = "orig_procs_extbnb_v00_p00.txt"
    out_procfile = "procs_extbnb.txt"
    tmp_inputlist_dir = "jobfilelists"
    larlite_simlink_dir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/simlinks/data_extbnb_v00_p00"
    supera_outdir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/supera/data_extbnb_v00_p00"
    ismc = False
elif project=="bnb":
    original_procs_file = "orig_bnb_procs.txt"
    out_procfile = "procs_bnb.txt"
    tmp_inputlist_dir = "jobfilelists"
    larlite_simlink_dir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/simlinks/data_bnb_v00_p00"
    supera_outdir = "/pnfs/uboone/persistent/users/tmw/dl_thrumu/supera/data_bnb_v00_p00"
    ismc = False

os.system("rm %s/*"%(tmp_inputlist_dir))

original = open(original_procs_file,'r')

lines = original.readlines()

original_procs = []
for l in lines:
    original_procs.append( int(l) )

output_files = os.listdir( supera_outdir )

processed = []
badlist = []
for output in output_files:
    if ".root" not in output:
        continue
    rfile = rt.TFile( supera_outdir+"/"+output )
    tree = rfile.Get("partroi_tpc_tree")
    try:
        nentries = tree.GetEntries()
    except:
        print supera_outdir+"/"+output," failed test"
        nentries = 0

    if nentries==0:
        badlist.append( supera_outdir+"/"+output )
        continue
    proc = int( output.split(".")[-2].split("_")[-1] )
    processed.append( proc )

out_proc = open( out_procfile, 'w' )

ftypes = ["chstatus","opreco","opdigit","wire"]
if ismc:
    ftypes.append("mcinfo")
    ftypes.append("simch")
remaining = []
for original in original_procs:
    if original not in processed:
        out_proc.write("%d\n"%(original))
        remaining.append( original )
        flist = open("%s/flist_%04d.txt"%(tmp_inputlist_dir,original),'w')
        for ftype in ftypes:
            print >> flist,os.path.realpath( "%s/larlite_%s_%04d.root"%(larlite_simlink_dir,ftype,original)  ),"larlite_%s_%04d.root"%(ftype,original)
        flist.close()

out_proc.close()

fbadlist = open("badlist.txt",'w')
for bad in badlist:
    print >> fbadlist,bad
fbadlist.close()

print "proc file made: ",out_procfile
print "badlist made: badlist.txt"
print "processes finished: ",len(processed)
print "processes remaining: ",len(remaining)
print "Remeber: copy your proc file to /pnfs, remove files from badlist and run checkana again so that project.py removes them when running makeup."
