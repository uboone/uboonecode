#! /bin/bash

source /cvmfs/uboone.opensciencegrid.org/bin/mpdsetenv.sh
env > env.txt
lar -c test_geometry_uboone.fcl > test_geometry_uboone.out 2> test_geometry_uboone.err
