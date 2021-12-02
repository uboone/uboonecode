#!/bin/sh
source nuwro_script.sh
init_gen_common.sh prodhepmc_uboone_10000.fcl  10000
lar -c run_PrimaryNeutronFilter.fcl -s genfile.root.local -o genfile_gen_filtered.root.local -n -1
rm genfile.root.local
lar -c wirecell_g4_uboone.fcl -s genfile_gen_filtered.root.local -o genfile_g4.root.local -n -1
rm genfile_gen_filtered.root.local
lar -c run_NeutronScatterFilter.fcl -s genfile_g4.root.local -o genfile_filtered.root.local -n -1
rm genfile_g4.root.local
