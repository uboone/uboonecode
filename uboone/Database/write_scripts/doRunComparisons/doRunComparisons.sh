#!/bin/bash

usage() { printf "\n>> Flags:\n>> Single run: ./doRunComparisons.sh -s run_no
>> Compare selected runs: ./doRunComparisons.sh -c run_no1,run_no2,run_no3
>> Compare a run range: ./doRunComparusins -r run_lowerlimit,run_higherlimit\n\n"; exit 1; }

while getopts "r:c:s:h:" opt; do
  case $opt in
    s)
      set -f # disable glob
      si=($OPTARG)
      it=0
      array[$it]=${si[$it]}
      ;;
    r)
      set -f # disable glob
      IFS=','
      ra=($OPTARG)
      it=0
      for (( i=${ra[0]}; i<=${ra[1]}; i++)) 
      do
        array[$it]=$i
        it+=1
      done
      ;;
    c)
      set -f # disable glob
      IFS=","
      co=($OPTARG)
      for (( it=0; it<${#co[@]}; it++)) 
      do
        array[$it]=${co[$it]}
      done
      ;;
    h)
      usage
      ;;
    *) 
      usage     
      ;;
  esac
done
shift $(( OPTIND - 1 ))

export elecDbVr="v1r2"
export chanStartDbVr="v1r4"

## clean up after possible previous running
rm -rf ./electronicsDbFiles
rm -rf ./chanStatDbFiles
rm -rf *.list

mkdir ./electronicsDbFiles
mkdir ./chanStatDbFiles

# for each run in array, reqeust date/time from SAM, convert to unix time, 
# and pull down database information
for runit in "${array[@]}" 
do
  export fileName=`samweb list-files "run_number=${runit}.0 and data_tier=raw and ub_project.stage=assembler"`
  export basicTime=`samweb get-metadata $fileName | grep "Start Time" | sed 's/.* //' | sed 's/.\{6\}$//'`
  export unixTime=`date -d $basicTime +"%s"`
  echo ">> Getting database information for run ${runit} which corresponds to unix time $unixTime"

  export elecWebString="http://dbdata0vm.fnal.gov:8086/uboonecon_dev/app/data?f=electronicscalib_data&t=${unixTime}&tag=${elecDbVr}"
  export chanStatWebString="http://dbdata0vm.fnal.gov:8186/uboonecon_prod/app/data?f=channelstatus_data&t=${unixTime}&tag=${chanStatDbVr}"
  export singleQt="'"

  export newElecWebString="${singleqt}${elecWebString}${singleqt}"
  export newChanStatWebString="${singleqt}${chanStatWebString}${singleqt}"

  wget -O ./electronicsDbFiles/run_${runit}.csv $newElecWebString
  wget -O ./chanStatDbFiles/run_${runit}.csv $newChanStatWebString
done

# delete files from runs that are empty

export emptyTester=`find ./electronicsDbFiles -size 0`
if [ -n "$emptyTester" ] 
then
  echo ">> Deleting empty files..."
  echo ">> If all selected runs are empty, you're gonna have a bad time."
  find ./electronicsDbFiles -size 0 -print0 |xargs -0 rm
  find ./chanStatDbFiles -size 0 -print0 |xargs -0 rm
fi

sh -c 'ls -d electronicsDbFiles/* > electronicsDbFiles.list'
sh -c 'ls -d chanStatDbFiles/* > chanStatDbFiles.list'

python compareRuns.py

mv electronicsDbFiles.list electronicsDbFiles/
mv chanStatDbFiles.list chanStatDbFiles/
