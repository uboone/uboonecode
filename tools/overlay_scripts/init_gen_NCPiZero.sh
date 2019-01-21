MY_OUTPUT_FILE=genfile.root.local
MY_FCL_FILE=prodgenie_bnb_nu_filtered_NCPiZero_uboone.fcl
MY_N_EVENTS=100000
MY_OUT_LOG=larInitGen.out
MY_ERR_LOG=larInitGen.err

echo "#include \"$MY_FCL_FILE\"" > local_gen.fcl
echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen.fcl
echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen.fcl

if [ -f $MY_OUTPUT_FILE ]; then
    echo "File $MY_OUTPUT_FILE exists, so not generating again."
else
    echo "File $MY_OUTOUT_FILE does not exist, so running generation of it..."
    echo "Running command 'lar -c local_gen.fcl -T ./genfile_hist.root -o $MY_OUTPUT_FILE -n $MY_N_EVENTS > $MY_OUT_LOG 2> $MY_ERR_LOG' "
    lar -c local_gen.fcl -T ./genfile_hist.root -o $MY_OUTPUT_FILE -n $MY_N_EVENTS
fi
