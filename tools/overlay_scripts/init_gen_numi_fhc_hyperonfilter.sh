#!/bin/sh
init_gen_common.sh prodgenie_numi_nu_uboone_fhc.fcl 5000
lar -c wirecell_g4_uboone.fcl -s genfile.root.local -o genfile_g4.root.local -n -1
lar -c run_HyperonFilter.fcl -s genfile_g4.root.local -o genfile_filtered.root.local -n -1
mv genfile_filtered.root.local genfile.root.local 
rm genfile_g4.root.local
