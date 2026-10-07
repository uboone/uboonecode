#! /bin/bash

source /cvmfs/uboone.opensciencegrid.org/bin/mpdsetenv.sh
env > env.txt
lar -c InteractTest.fcl > InteractTest.out 2> InteractTest.err
