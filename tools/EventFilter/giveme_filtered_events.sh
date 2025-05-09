#!/bin/sh
# giveme_filtered_events.sh
# Main Author: Marina Reggiani-Guzzo 
# Tweaked by C Thorpe :) 

# USAGE: giveme_filtered_events.sh <eventlist> <samdef> 
# Make sure to update the input fields below

# Description:
# Purpose of the script is to convert a list of run subrun event list and filter
# these events into one or many root files which you can easily use for scanning 
# events in the event display. 

# How it works:
# The script takes the run subrun event list and if it is a large number of events,
# then it splits it into n pieces where n is the number of events you specify by
# the nmaxevents variable. It then runs a sam query to find the parent files that
# are associated with those events. Be sure to configure the sam definition to be
# the one where the run subrun event list will have come from and the type which 
# tunes the sam query slightly. It then makes a sam defintion with all these files,
# prestages these files so we can access them. Next a filelist with the full path to 
# the files is created. We then run the event filter module which sives through all the
# files and grabs out the events you are interesed in and puts them into a root file.
# enjoy!

# Useful note: Run this with nohup/screen/tmux since it prestages the files you give it

# ---------------------------------------------------------------------------------

# input: please check if variables correspond to your case

# format: "run subrun event"
runsubrunlist=$1

# The name of the sam def which you made the run subrun event list from
originalsamdef=$2

# choose a name for your sam def, any name as long as it doesnt already exist
newsamdef=$2_tmp

# 0=data, 1=offbeam, 2=overlay
#type=2

# max number of events per new samdef created
nmaxevents=200

# keep the auxiliary files created during the process? 0=keep, 1=delete
deletefiles=1

# ---------------------------------------------------------------------------------







# ---------------------------------------------------------------------------------
BLUE="\033[1;34m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
DEFAULT="\033[0m"
# ---------------------------------------------------------------------------------
# because there is a limit of arguments to create a new samdef, it is necessary to 
# split the runsubrunfile into smaller files
# ---------------------------------------------------------------------------------

echo -e "${BLUE}Splitting up run subrun event lists...${DEFAULT}"

nlines=`wc -l < ${runsubrunlist}`
if [ nlines > $nmaxevents ] ; then # split file

    runsubrunlist_new="$(basename ${runsubrunlist} .txt)_split"

    split -l $nmaxevents ${runsubrunlist} ${runsubrunlist_new} # split input file

    for i in $runsubrunlist_new* ; do # loop to add extention to the splitted files
    mv $i $i.txt
    done

    # List the split the file names in the current directory
    ls $runsubrunlist_new*

    runsubrunlist_arr=$( ls $runsubrunlist_new* ) # array of splitted files

else

    runsubrunlist_arr=$runsubrunlist

fi

echo 

# --------------------------------------------------------------------------------- 
# crete a file list by querying the input sam def
# --------------------------------------------------------------------------------- 

declare -a query_array=("and data_stream beam_good" "and data_stream all" "") # each type of data require a different query

# Create a master filelist to catch any duplicate files
touch master_filelist.txt

echo -e "${BLUE}Creating the list of art-root files and checking for duplicate files from the sam query...${DEFAULT}"

for f in ${runsubrunlist_arr[@]}; do

    counter=0

    for line in $(cat $f); do

    if [ $counter -eq 0 ]; then
            run=$line
    elif [ $counter -eq 1  ]; then
            subrun=$line
            runsubrun="$run.$subrun"
            
            # Call sam and get the file. This can give multiple files e.g. beam good and bad files
            #tempfiles=$(samweb list-files "defname:$originalsamdef and run_number $runsubrun ${query_array[${type}]}")
            tempfiles=$(samweb list-files "defname:$originalsamdef and run_number $runsubrun")
            
            # Loop over the files 
            for tempfile in $(echo "$tempfiles"); do
                
                # Add the file to the master list (this is for identifying duplicate files)
                echo "$tempfile" >> master_filelist.txt
		echo "$tempfile"

                # Only add the file if it hasnt come up before
                checkduplicate=$(cat master_filelist.txt | grep $tempfile)
		echo $($check_duplicate)
                if [ $(echo $check_duplicate | wc -l) -eq 1 ]; then
                    echo "$tempfile" >> files_$(basename $f)
                else 
                    echo "Found a duplicate file (likely has multiple run subruns in the list specified) so not adding this to the list."
                fi
            
            done
           
    fi

    counter=$(($counter+1))
    if [ $counter -eq 3 ]; then counter=0; fi

    done

done
    
fileslist_arr=$( ls files_* ) # array with the name of the files that contain the artroot file name of the selected events

echo

# ---------------------------------------------------------------------------------
# create a samdef from the files and prestage the files
# ---------------------------------------------------------------------------------

declare -a newsamdef_arr=()
counter=0

for i in ${fileslist_arr[@]}; do

    newsamdeftemp=${newsamdef}_${counter}
    newsamdef_arr+=( $newsamdeftemp )

    files=`awk '{print $newsamdef}' $i | paste -s -d, -`
    files=$files

    samweb create-definition $newsamdeftemp "file_name $files"

    echo
    echo -e "${BLUE}samweb prestage-dataset --defname="$newsamdeftemp" --parallel=4${DEFAULT}";
    samweb prestage-dataset --defname="$newsamdeftemp" --parallel=4
    echo 
    echo "Prestage complete!"
    echo

    counter=$(($counter+1))
    
done

# ---------------------------------------------------------------------------------
# create a filelist with the full path
# ---------------------------------------------------------------------------------

declare -a filelist_fullpath=()
counter=0

echo 
echo -e "${BLUE}Creating the filelists with full paths...${DEFAULT}"

for def in ${newsamdef_arr[@]}; do

    for files in `samweb list-files defname: $def`; do

    # locate the file
    newfile=`samweb locate-file "${files}"`
    #echo "Newfile:"
    #echo "$newfile"

    good_location=""

    # If file has more than one location, check which is online
    while IFS= read -r line; do

      # get rid of enstore/dcache in the name
      line=${line#"enstore:"}
      line=${line#"dcache:"}

      # get rid of anything after the bracket
      line=${line%(*}
     
      # Check which location is staged 
      location=$(cat $line/."(get)(${files})(locality)")
      location=$(echo "$location" | head -n 1)
      if [ "ONLINE" = "$location" ] || [ "$location" = "ONLINE_AND_NEARLINE" ]; then
        good_location=$line
      fi

    done <<< "$newfile" 

    # Put into the new file                                                     
    echo ${good_location}/${files} >> ${def}_file.list

    done

    echo
    echo -e "${BLUE}${def}_file.list CREATED${DEFAULT}"

    filelist_fullpath+=( ${def}_file.list )

done

# ---------------------------------------------------------------------------------
# running fhicl file
# ---------------------------------------------------------------------------------

counter=0

# copy the event fcl file over (hopefully this fcl wont ever get deleted)
#cp /exp/uboone/app/users/cthorpe/throwaway/test/run_EventFilter.fcl  .
fhicl_loc=$(find_fhicl.sh run_EventFilter.fcl | tail -n 1)
cp $fhicl_loc . 

echo
echo
echo -e "${BLUE}Running the Event Filter...${DEFAULT}"

for s in ${runsubrunlist_arr[@]}; do

    # first delete the following lines
    sed -i '/^physics.filters.filter.EventList/d' run_EventFilter.fcl
    sed -i '/^physics.filters.filter.Selection/d' run_EventFilter.fcl

    # adding those lines again with new file
    echo "physics.filters.filter.EventList: \"${s}\"" >> run_EventFilter.fcl
    echo "physics.filters.filter.Selection: 1" >> run_EventFilter.fcl
    
    echo -e "${RED}lar -c run_EventFilter.fcl -S ${filelist_fullpath[$counter]} -o filtered_${newsamdef_arr[$counter]}.root${DEFAULT}"
    eval "lar -c run_EventFilter.fcl -S ${filelist_fullpath[$counter]} -o filtered_${newsamdef_arr[$counter]}.root"
    echo

    counter=$(($counter+1))

done


# ---------------------------------------------------------------------------------
# print a summary of the events you have just filtered
# ---------------------------------------------------------------------------------

echo -e "${YELLOW}-------------SUMMARY"
echo "samdef:                          ${originalsamdef}"
echo "number of events:                ${nlines}"
echo "new samdef:                      ${newsamdef_arr[@]}"
echo -e "--------------------${DEFAULT}"

# ---------------------------------------------------------------------------------
# delete auxiliary files?
# ---------------------------------------------------------------------------------

# echo ${#runsubrunlist_arr[@]}

if [ $deletefiles -eq 1 ]; then
    echo
    echo "Removing file lists..."
    rm $newsamdef*.list 
    rm *split*.txt
    rm run_EventFilter.fcl
    rm master_filelist.txt
    rm $nmaxevents # have no idea why a file of this name gets made... but lets delete it!
fi

# Cleanup the sam definitions made
for def in ${newsamdef_arr[@]}; do
    echo "samweb delete-definition $def"
    samweb delete-definition $def
done

echo 
echo "FINSHED!"
