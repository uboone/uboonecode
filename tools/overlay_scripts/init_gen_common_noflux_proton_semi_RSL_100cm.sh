#!/bin/sh

MY_FCL_FILE=prod_proton_overlay_semi_ana_RSL100cm.fcl #prodgenie_bnb_nu_filtered_NCPiZero_uboone.fcl
MY_N_EVENTS=50
MY_OUTPUT_FILE=${3:-genfile.root.local}
MY_OUT_LOG=larInitGen.out
MY_ERR_LOG=larInitGen.err

echo "#include \"$MY_FCL_FILE\"" > local_gen.fcl
#echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen.fcl
#echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen.fcl
echo "services.IFDH: {}" >> local_gen.fcl

echo "#include \"$MY_FCL_FILE\"" > local_gen_include.fcl
#echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen_include.fcl
#echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen_include.fcl
echo "services.IFDH: {}" >> local_gen_include.fcl
echo "gen_detail: { physics: {@table::physics} services: {@table::services} outputs: {@table::outputs} source: {@table::source} process_name: @local::process_name }" >> local_gen_include.fcl

if [ -f $MY_OUTPUT_FILE ]; then
    echo "File $MY_OUTPUT_FILE exists, so not generating again."
else
    echo "File $MY_OUTOUT_FILE does not exist, so running generation of it..."
    echo "Running command 'lar -c local_gen.fcl -T ./genfile_hist.root -o $MY_OUTPUT_FILE -n $MY_N_EVENTS > $MY_OUT_LOG 2> $MY_ERR_LOG' "
    lar -c local_gen.fcl -T ./genfile_hist.root -o genfile.root.temp -n $MY_N_EVENTS
    lar -c standard_overlay_gen_SubRunPOTInEvent.fcl -s genfile.root.temp -T ./genfile_pot_hist.root -o $MY_OUTPUT_FILE -n -1
    rm genfile.root.temp
fi

