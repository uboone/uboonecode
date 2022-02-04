#!/bin/bash

# Get the file(s) with the neutrons
filestouse=5

# Count the number of files 
GenFiles=$(wc -l  gen_neutrons.list | awk '{ print $1 }')

touch genfiles.list

i=1
while [ "$i" -le "$filestouse" ]; do
LINE=$((1 + $RANDOM % $GenFiles))
FILE=$(sed "${LINE}q;d" gen_neutrons.list)
#Copy that file onto the node
ifdh cp $FILE neutrons_gen_$i.root
#echo neutrons_gen_$i.root >> genfiles.list
lar -c wirecell_g4_uboone.fcl -s neutrons_gen_$i.root -o neutrons_g4_$i.root -n -1
lar -c run_NeutronScatterFilter.fcl -s  neutrons_g4_$i.root -o neutrons_filtered_$i.root -n -1
echo neutrons_filtered_$i.root >> genfiles.list
rm neutrons_gen_$i.root
rm neutrons_g4_$i.root
i=$(($i + 1))
done

lar -c standard_overlay_gen_SubRunPOTInEvent.fcl -S genfiles.list -o genfile_filtered.root.local -n -1

i=1
while [ "$i" -le "$filestouse" ]; do
echo REMOVING neutrons_filtered_$i.root
rm neutrons_filtered_$i.root
i=$(($i + 1))
done
