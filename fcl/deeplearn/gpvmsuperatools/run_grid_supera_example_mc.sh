#!/bin/bash

# process ID file
procfile=/pnfs/uboone/persistent/users/tmw/dl_thrumu/simlinks/procs_mcc7_cosmic_v00_p00.txt

# supera config file path
superacfgpath=/pnfs/uboone/persistent/users/tmw/supera_mcc7_bnb_cosmic.cfg

# supera cfg file
superacfgname=supera_mcc7_bnb_cosmic.cfg

# input file list directory
inputlistdir=/pnfs/uboone/persistent/users/tmw/dl_thrumu/jobfilelists/mcc7_cosmic_v00_p00

# copy over the procs file
ifdh cp $procfile procs.txt

# copy over config
ifdh cp $superacfgpath .

# get the process ID
export LARLITE_FILEID=`sed -n ${PROCESS}p procs.txt`
echo "Running job with processid=${PROCESS} and fileid=${LARLITE_FILEID}"

chstatus_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_chstatus_%04d.root",$1) }'`
wire_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_wire_%04d.root",$1) }'`
opdigit_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_opdigit_%04d.root",$1) }'`
opreco_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_opreco_%04d.root",$1) }'`
mcinfo_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_mcinfo_%04d.root",$1) }'`
simch_file=`echo ${LARLITE_FILEID} | awk '{ printf("larlite_simch_%04d.root",$1) }'`
supera_file=`echo ${LARLITE_FILEID} | awk '{ printf("supera_%04d.root",$1) }'`

# get flist
flist=`echo ${LARLITE_FILEID} | awk '{ printf("${inputlistdir}/flist_%04d.txt",$1)}'`

echo "Copying the following files:"
echo $chstatus_file
echo $wire_file
echo $mcinfo_file
echo $simch_file
echo $opdigit_file
echo $opreco_file
echo $supera_file
echo $flist

# copy the filelist
ifdh cp $flist flist.txt

# copy the larlite files
ifdh cp -f flist.txt

supera $superacfgname $supera_file $wire_file $opdigit_file $opreco_file $mcinfo_file $simch_file

rm $chstatus_file 
rm $wire_file
rm $opdigit_file
rm $opreco_file
rm $mcinfo_file

echo "FINISHED SUPERA"
ls -lh