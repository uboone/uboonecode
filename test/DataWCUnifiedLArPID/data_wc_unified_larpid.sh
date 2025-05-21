#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Exit if stash cache isn't mounted.

UBOONE_EXAMPLE_DATA_DIR=/pnfs/uboone/persistent/stash/uboone_example_data
if [ ! -d $UBOONE_EXAMPLE_DATA_DIR ]; then
  echo "Quittig because stash cache isn't available."
  exit
fi

# Exit if pnfs persistent isn't mounted.

if [ ! -d /pnfs/uboone/persistent ]; then
  exit
fi

# Exit if ups product wcp isn't set up.

if [ x$WCP_DIR = x ]; then
  exit
fi

# LD_LIBRARY_PATH kluge

export LD_LIBRARY_PATH=`dropit -p $LD_LIBRARY_PATH ublarcvapp`

# Set wire cell path.

export WIRECELL_PATH=${UBOONEDATA_DIR}/WireCellData:${WIRECELL_FQ_DIR}/share/wirecell

# Maybe update execute path.

if [ -d $UBANA_DIR/ubana/MicroBooNEWireCell/utils ]; then
  export PATH=$UBANA_DIR/ubana/MicroBooNEWireCell/utils:$PATH
fi
if [ -d $UBCV_DIR/ubcv/LArCVImageMaker/mcc10_prod_scripts ]; then
  export PATH=$UBCV_DIR/ubcv/LArCVImageMaker/mcc10_prod_scripts:$PATH
fi

# Set up python path.

export PYTHONPATH=`pwd`:$UBUTIL_DIR/python:$LARBATCH_DIR/python:$PYTHONPATH
rm -rf project_modules
cp -r $LARBATCH_DIR/python project_modules
touch project_modules/__init__.py

# Set experiment environment variables (not set by mrbsetenv, but needed by IFDH).

export EXPERIMENT=uboone
export SAM_EXPERIMENT=uboone

rm -f *.root
rm -f *.json
rm -f *.rootstat
rm -f *.fcl
rm -f *.out
rm -f *.err
rm -f *.db
rm -rf WCPwork

# This script runs the WC unified chain (includes reco1.5 and reco2, starting from reco1).

input=$UBOONE_EXAMPLE_DATA_DIR/reco1/PhysicsRun-2018_7_2_21_23_1-0017505-00141_20180714T132533_ext_bnb_3_20181207T001756_optfilter_20181224T035203_reco1_postwcct_postdl_20181224T035451.root
for FCL in run_celltreeub_port_prod.fcl \
  run_slimmed_port_data_sp.fcl \
  mcc10_dlreco_w_wirecell_driver_data_lantern.fcl \
  run_wcpplus_port.fcl \
  run_wcpf_port.fcl
do

  # Make a copy of the fcl file in the current directory.
  # Needed by some init source scripts.

  echo
  echo "Running stage $FCL"
  p=`find \`echo $FHICL_FILE_PATH | tr : ' '\` -name $FCL 2> /dev/null | head -1`
  echo $p
  if [ x$p = x ]; then
    echo "Fcl file $FCL not found."
    exit 1
  fi
  cp $p $FCL
  fhicl-dump $FCL > cfg-$FCL

  out=`basename $FCL .fcl`.out
  err=`basename $FCL .fcl`.err
  end_out=end_`basename $FCL .fcl`.out
  end_err=end_`basename $FCL .fcl`.err

  # Source init scripts.

  inits=''
  if [ $FCL = run_slimmed_port_data_sp.fcl ]; then
    inits=update_slimmed_port_data.sh
  elif [ $FCL = run_wcpf_port.fcl ]; then
    inits=update_wcpf_port.sh
  fi
  if [ x$inits != x ]; then
    inits_p=`which $inits`
    echo "Sourcing $inits_p"
    . $inits_p
    stat=$?
    echo "Source script finished with status $stat"
    if [ $stat -ne 0 ]; then
      exit $stat
    fi
  fi

  # Run lar.

  if [ x$input = x ]; then
    cmd="lar --rethrow-all -c $FCL -n 2"
  else
    cmd="lar --rethrow-all -c $FCL -s $input -n 2"
  fi
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi

  # Run end scripts.

  ends=''
  if [ $FCL = run_celltreeub_port_prod.fcl ]; then
    ends=fully_unified_reco2_wirecell.sh

    # Generate fake metadata.

    cat <<EOF > celltreeDATA.root.json
{
  "event_count": 2,
}
EOF

  elif [ $FCL = mcc10_dlreco_w_wirecell_driver_data_lantern.fcl ]; then
    ends=merge_dlreco_fnal_data_lantern.sh
  fi
  if [ x$ends != x ]; then
    echo "Running end script $ends"
    $ends > $end_out 2> $end_err
    stat=$?
    echo "End script finished with status $stat"
    if [ $stat -ne 0 ]; then
      exit $stat
    fi
  fi

  #if [ -f $output ]; then
  #  input=$output
  #fi
  input=`ls -t1 *.root | egrep -v 'celltree|hist|larlite|larcv|Supplemental|TGraphs|merged_dlreco' | head -n1`
done

# Done (success).

#for root in *.root
#do
  #rootstat.py $root > ${root}.rootstat
  #rm $root   # Clean up.
#done

exit 0


