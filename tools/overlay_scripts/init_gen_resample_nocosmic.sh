#!/bin/bash

############################################################################################
# Generate events by selecting a set of files from a list of generator-level files, 
# passing them through g4 stage and applying a filter. This is a variation on the script
# designed for generation of no-cosmic samples
# Author: C Thorpe (Manchester U)
#
# Useage:
# init_gen_resample_nocosmic.sh <number of files to sample> <intput filelist> <optional filter fhicl>
############################################################################################

# Parse input parameters
filestouse=$1
inputfilelist=$2
filterfhicl=$3

# Count the number of files in file list to sample from 
GenFiles=$(wc -l  $inputfilelist | awk '{ print $1 }')

touch genfiles.list

i=1
while [ "$i" -le "$filestouse" ]; do

  # Pick a file from list
  LINE=$((1 + $RANDOM % $GenFiles))
  FILE=$(sed "${LINE}q;d" $inputfilelist)

  #Copy that file onto the node and run g4 over it
  ifdh cp $FILE gen_$i.root

  # Make list of all gen files, remove gen file from node
  echo gen_$i.root >> genfiles.list
  #rm gen_$i.root
  i=$(($i + 1))

done

# Run G4 over the list of random files
lar -c wirecell_g4_uboone.fcl -S genfiles.list -o g4.root -n -1

# If post-g4 filter added - run this here
if [ ! -z "$filterfhicl" ]; then
    echo "Running filter ${filterfhicl}"
    lar -c ${filterfhicl} -s g4.root -o filter.root -n -1
    mv filter.root g4.root
fi

i=1
while [ "$i" -le "$filestouse" ]; do
  rm gen_$i.root
  i=$(($i + 1))
done
