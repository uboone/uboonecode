#!/bin/sh

MY_FRST_FCL_FILE=$1 #prodgenie_bnb_nu_uboone.fcl
MY_SCND_FCL_FILE=$2 #prodgenie_bnb_nu_wirecell_g4_uboone_filtered_CCmuNoPi.fcl
MY_N_EVENTS=$3 #2500
MY_OUTPUT_FILE=${4:-genfile.root.local}
MY_OUT_LOG=larInitGen.out
MY_ERR_LOG=larInitGen.err

echo "#include \"$MY_FRST_FCL_FILE\"" > local_gen.fcl
echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen.fcl
echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen.fcl
echo "services.IFDH: {}" >> local_gen.fcl

echo "#include \"$MY_SCND_FCL_FILE\"" > local_gen_include.fcl
echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen_include.fcl
echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen_include.fcl
echo "services.IFDH: {}" >> local_gen_include.fcl
echo "gen_detail: { physics: {@table::physics} services: {@table::services} outputs: {@table::outputs} source: {@table::source} process_name: @local::process_name }" >> local_gen_include.fcl

if [ -f $MY_OUTPUT_FILE ]; then
    echo "File $MY_OUTPUT_FILE exists, so not generating again."
else
    echo "File $MY_OUTOUT_FILE does not exist, so running generation of it..."
    lar -c local_gen.fcl -T ./first_hist.root -o first.root.temp -n $MY_N_EVENTS
    lar -c $MY_SCND_FCL_FILE -T ./second_hist.root -s first.root.temp -o second.root.temp -n -1
    lar -c standard_overlay_gen_SubRunPOTInEvent.fcl -s second.root.temp -T ./genfile_pot_hist.root -o $MY_OUTPUT_FILE -n -1
    rm first.root.temp second.root.temp
fi

