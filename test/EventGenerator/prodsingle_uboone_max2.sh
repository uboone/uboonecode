#! /bin/bash

source /cvmfs/uboone.opensciencegrid.org/bin/mpdsetenv.sh
env > env.txt
lar -c prodsingle_uboone_max2.fcl > prodsingle_uboone_max2.out 2> prodsingle_uboone_max2.err

