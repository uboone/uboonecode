#!/bin/bash

############################################################################################
# Generate events by selecting a set of files from a list of generator-level files, 
# passing them through g4 stage and applying a filter
# Author: C Thorpe (Manchester U)
#
# Useage:
# init_gen_resample.sh <number of files to sample> <intput filelist> <optional filter fhicl>
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
  lar -c wirecell_g4_uboone.fcl -s gen_$i.root -o g4_$i.root -n -1

  # If post-g4 filter added - run this here
  if [ ! -z "$filterfhicl" ]; then
      echo "Running filter ${filterfhicl}"
      lar -c ${filterfhicl} -s g4_$i.root -o filter_$i.root -n -1
      mv filter_$i.root g4_$i.root
  fi

  # Make list of all g4 files, remove gen file from node
  echo g4_$i.root >> genfiles.list
  rm gen_$i.root
  i=$(($i + 1))

done

lar -c standard_overlay_gen_SubRunPOTInEvent.fcl -S genfiles.list -o genfile_filtered.root.local -n -1

i=1
while [ "$i" -le "$filestouse" ]; do
  rm g4_$i.root
  i=$(($i + 1))
done
