#!/bin/bash

# Get the file(s) with the hyperons
filestouse=10

# Count the number of files 
GenFiles=$(wc -l  gen_hyperons.list | awk '{ print $1 }')

touch genfiles.list

i=1
while [ "$i" -le "$filestouse" ]; do
LINE=$((1 + $RANDOM % $GenFiles))
FILE=$(sed "${LINE}q;d" gen_hyperons.list)
#Copy that file onto the node
ifdh cp $FILE hyperons_gen_$i.root
#echo hyperons_gen_$i.root >> genfiles.list
lar -c wirecell_g4_uboone.fcl -s hyperons_gen_$i.root -o hyperons_g4_$i.root -n -1
#lar -c run_HyperonFilterG4.fcl -s  hyperons_g4_$i.root -o hyperons_filtered_$i.root -n -1
echo hyperons_g4_$i.root >> genfiles.list
rm hyperons_gen_$i.root
#rm hyperons_g4_$i.root
i=$(($i + 1))
done

lar -c standard_overlay_gen_SubRunPOTInEvent.fcl -S genfiles.list -o genfile_filtered.root.local -n -1

i=1
while [ "$i" -le "$filestouse" ]; do
echo REMOVING hyperons_filtered_$i.root
rm hyperons_g4_$i.root
i=$(($i + 1))
done
