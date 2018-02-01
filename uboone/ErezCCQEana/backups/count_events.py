# combine *vertices.csv files from a grid job
default_name = "prod_reco2_extbnb_v8_mcc8"
name = raw_input("enter list name:...<"+default_name+">") or default_name

defauld_indirname = "/uboone/data/users/ecohen/book/mcc8/v06_42_00/ccqe_ana/"+name+"_CCQE/"
indirname = raw_input("enter indir:...<"+defauld_indirname+">") or defauld_indirname

print 'counting events from',name
import os
#import ROOT
list_dir = os.listdir(indirname)
print "list of directories to copy from:",list_dir

file_ctr,events_ctr = 0,0
for dir_name in list_dir:#{
    if os.path.isdir(indirname+'/'+dir_name): #{
        print 'in directory',dir_name
        list_files = os.listdir(indirname+'/'+dir_name)
        for file in list_files:#{
            if "larStage0.out" in file:#{
                print 'in file',indirname+'/'+dir_name+'/'+file
                with open(indirname+'/'+dir_name+'/'+file) as searchfile:
                    for line in searchfile:
                        left,sep,right = line.partition("TrigReport Events total = ")
                        if sep: # True iff 'Figure' in line
                            n_events_str = right.split(' ',1)[0]
                            n_events = int(n_events_str)
                            print "int(n_events_str):",int(n_events_str)
                print 'adding %d events'%n_events
                events_ctr = events_ctr + n_events
                file_ctr = file_ctr + 1
            #}
        #}
    #}
#}
print 'done counting events from',file_ctr,'files',', the total number of processed events is:\n',events_ctr
print
