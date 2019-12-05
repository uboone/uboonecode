# Production Fhicl files for DL MCC9 Production

The DL reconstruction chain is run in one job, but proceeds within that job, in two stages.

The first stage is the execution of the art producer module, `DLIntegration_module`,
located in the `ubcv` respository.  
This module is responsible for 
  * making `larcv` and `larlite` files used by the DL reconstruction chain, and
  * running convolution neural networks and storing the output into both `larcv` format and into the art file.

The second stage of the job is to execute the DL reconstruction algorithms.
The algorithms use the `larcv` and `larlite` files produced in the first stage.
The results of these algorithms are several different root files, which are merged together to ease 
file book-keeping. 

Note: the DL algorithms are not run in the larsoft environment.
Instead, the algorithms need to be executed using an `endscript` for the job.

Therefore, configuring a job involves specifying a larsoft fcl file for the first stage 
and then the bash script for running the second stage.

## Configurations

We specify the driver fcl files and DL reco. script for each type of data.
Contact Taritree (twongj01@tufts.edu or better on slack), if a data type is missing.

| Data Type    | input stage required     | larsoft driver fcl file       | DL reco. script                           |
| ------------ | ------------------------ | ----------------------------- | ----------------------------------------- |
| Overlay MC   | reco1 (needs simchannel) |  mcc9_dlreco_driver_mc.fcl    |  rundlreco_mcoverlay_ssnetvertexonly.sh   |
| Full MC      | reco1 (needs simchannel) |  mcc9_dlreco_driver_mc.fcl    |  rundlreco_fullmc_ssnetvertexonly.sh      |
| BNB data     | reco1 or reco2           |  mcc9_dlreco_driver_data.fcl  |  rundlreco_bnb_ssnetvertexonly.sh         |
| EXT-BNB data | reco1 or reco2           |  mcc9_dlreco_driver_data.fcl  |  rundlreco_extbnb_ssnetvertexonly.sh      |
| EXT-unbiased | reco1 or reco2           |  mcc9_dlreco_driver_data.fcl  |  rundlreco_extunbiased_ssnetvertexonly.sh |

**WARNING** scripts still need to be finished.

## fcl chain

For data jobs `EXTBNB` and `BNB` data, one can run on Reco 2 files. 
The fcl chain is:

```
  <stage name="ssnet">
    <fcl>mcc9_dlreco_driver_data.fcl</fcl>
    <fcl>standard_dlreco_uboone_metadata.fcl</fcl>
    <endscript>rundlreco_bnb_ssnetvertexonly.sh</endscript>
```

For mc jobs, we want access to simchannels so we can have truth images in the same file.
However, this means we have to run some `reco2` processes.

The fcl chain for `OVERLAY` and `PURE-MC` samples are:

```
  <fcl>run_eventweight_microboone_justSplines.fcl</fcl>
  <fcl>wirecell_detsim_optical_overlay_uboone.fcl</fcl>
  <fcl>standard_overlay_optical_uboone.fcl</fcl>
  <fcl>reco_uboone_mcc9_8_driver_overlay_optical.fcl</fcl>
  <fcl>dlreco2.fcl</fcl>
  <fcl>dl_driver_overlay.fcl</fcl>
  <fcl>standard_dlreco_uboone_metadata.fcl</fcl>
```

One needs to select a driver fcl file and endscript using the table above.  
One also needs to include `standard_dlreco_uboone_metadata.fcl`. Note, this lives in `ubcv/ubcv/ubdlintegration/`.