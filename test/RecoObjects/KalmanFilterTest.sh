#! /bin/bash

source /cvmfs/uboone.opensciencegrid.org/bin/mpdsetenv.sh
env > env.txt
lar -c KalmanFilterTest.fcl > KalmanFilterTest.out 2> KalmanFilterTest.err
