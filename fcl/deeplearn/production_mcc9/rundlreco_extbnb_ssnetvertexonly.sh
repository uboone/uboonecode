#!/bin/bash


echo "<<<<<< RUN DL RECO SCRIPT >>>>>>"

# OUTPUT FILES FROM PREVIOUS STAGE
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh

echo "<< FILES available >> "
ls -lh

SUPERA=out_larcv_test.root  # has adc image, chstatus, ssnet output, mrcnn
OPRECO=larlite_opreco.root
RECO2D=larlite_reco2d.root

# HERE's OUR HACK: bring down ubdl, bring up dllee_unified
unsetup ubdl

echo "<<< SETUP DLLEE_UNIFIED >>>"
setup dllee_unified v1_0_1 -q e17:prof

# SETUP ENV FOR TAGGER BIN
export PATH=$LARLITECV_BASEDIR/app/TaggerCROI/bin:$PATH

echo "<<< CHECK ENV AFTER DLLEE_UNIFIED >>>"
export 

echo "<<< PRIMARY CHAIN >>>"
echo "< RUN TAGGER >"
TAGGER_CONFIG=$DLLEE_UNIFIED_DIR/dlreco_scripts/tagger_configs/tagger_bnbdata_v2_splity_mcc9.cfg
echo "TAGGER CONFIG: ${TAGGER_CONFIG}"
ls out_larcv_test.root > input_larcv.txt
ls larlite_opreco.root > input_larlite.txt
run_tagger $TAGGER_CONFIG

TAGGER_LARCV=tagger_anaout_larcv.root
TAGGER_LARLITE=tagger_anaout_larlite.root


echo "<<< RUN VERTEXER >>>"
VERTEX_CONFIG=$DLLEE_UNIFIED_DIR/dlreco_scripts/vertex_configs/prod_fullchain_mcc9ssnet_combined_newtag_extbnb_c10_union.cfg
#VERTEX_CONFIG=prod_fullchain_mcc9ssnet_combined_newtag_extbnb_c10_union.cfg # for debug
python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/run_vertexer.py -c $VERTEX_CONFIG -a vertexana.root -o vertexout.root -d ./ $SUPERA $TAGGER_LARCV
VERTEXOUT=vertexout.root
VERTEXANA=vertexana.root

echo "<<< RUN TRACKER >>>"
TRACKER_CONFIG=$DLLEE_UNIFIED_DIR/dlreco_scripts/tracker_configs/tracker_read_cosmo.cfg
python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/run_tracker_reco3d.py -c $TRACKER_CONFIG -i $SUPERA -t $TAGGER_LARCV -p $VERTEXOUT -d ./ 
TRACKEROUT=tracker_reco.root
TRACKERANA=tracker_anaout.root
mv -f tracker_reco_0.root $TRACKEROUT
mv -f tracker_anaout_0.root  $TRACKERANA

#echo "<<< TRACKONLY CHAIN >>>"
#echo "< RUN TRACKONLY VERTEXER >"
#TRKONLY_VERTEX_CONFIG=$DLLEE_UNIFIED_DIR/dlreco_scripts/vertex_configs/prod_fullchain_alltracklabel_combined_newtag_extbnb_c10_union.cfg
##TRKONLY_VERTEX_CONFIG=prod_fullchain_alltracklabel_combined_newtag_extbnb_c10_union.cfg # for debug
##python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/run_vertexer.py -c $TRKONLY_VERTEX_CONFIG -a vertexana_trackonly_temp.root -o vertexout_trackonly.root -d ./ $SUPERA $TAGGER_LARCV
#TRKONLY_VERTEXOUT=vertexout_trackonly.root
#TRKONLY_VERTEXANA=vertexana_trackonly.root
##python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/rename_vertexana.py vertexana_trackonly_temp.root $TRKONLY_VERTEXANA

#echo "<<< RUN TRACKONLY TRACKER >>>"
#TRKONLY_TRACKER_CONFIG=$DLLEE_UNIFIED_DIR/dlreco_scripts/tracker_configs/tracker_read_cosmo_trackonlyvertexer.cfg
##TRKONLY_TRACKER_CONFIG=tracker_read_cosmo_trackonlyvertexer.cfg
##python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/run_tracker_reco3d.py -c $TRKONLY_TRACKER_CONFIG -i $SUPERA -t $TAGGER_LARCV -p $TRKONLY_VERTEXOUT -d ./ 
#TRKONLY_TRACKEROUT=tracker_reco_trackonly.root
#TRKONLY_TRACKERANA=tracker_anaout_trackonly.root
##mv -f tracker_reco_0.root    $TRKONLY_TRACKEROUT
##mv -f tracker_anaout_0.root  $TRKONLY_TRACKERANA

echo "<<< RUN SHOWER RECO >>>"
echo "  < make inter file > "
NUEID_INTER_DIR=${LARLITECV_BASEDIR}/app/LLCVProcessor/InterTool/Sel/NueID/mac/ # using ups
#NUEID_INTER_DIR=./intertool_configs/ # using local scripts for debug
NUEID_INTER_CONFIG=${NUEID_INTER_DIR}/inter_nue_mc_pl2.cfg
echo "python ${NUEID_INTER_DIR}/inter_ana_nue_server.py -c ${NUEID_INTER_CONFIG} -mc -d -id 0 -od ./ -re larlite_reco2d.root vertexout.root"
python ${NUEID_INTER_DIR}/inter_ana_nue_server.py -c ${NUEID_INTER_CONFIG} -mc -d -id 0 -od ./ -re larlite_reco2d.root vertexout.root

echo "  < shower file >"
SHOWER_MAC_DIR=${LARLITECV_BASEDIR}/app/LLCVProcessor/DLHandshake/mac/ # using ups
#SHOWER_MAC_DIR=./shower # using local copy for debug
SHOWER_RECO_CONFIG=$SHOWER_MAC_DIR/config_nueid.cfg
SHOWER_RECO_DQDS=$SHOWER_MAC_DIR/dqds_mc_xyz.txt
python ${SHOWER_MAC_DIR}/reco_recluster_shower.py -c $SHOWER_RECO_CONFIG -mc -id 0 -od ./ --reco2d larlite_reco2d.root -dqds $SHOWER_RECO_DQDS nueid_lcv_out_0.root

echo "<< combine larlite files >>"
python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/combine_larlite.py -o larlite_dlmerged.root larlite_opreco.root larlite_reco2d.root tagger_anaout_larlite.root tracker_reco.root nueid_ll_out_0.root shower_reco_out_0.root
echo "<<< HADD ROOT FILES >>>"
hadd -f merged_dlreco.root $VERTEXOUT $VERTEXANA $TRACKERANA nueid_lcv_out_0.root larlite_dlmerged.root
echo "<<< Append UBDL Products >>>"
python $DLLEE_UNIFIED_DIR/dlreco_scripts/bin/append_ubdlproducts.py merged_dlreco.root out_larcv_test.root

echo "<<< cleanup excess root files >>>"
# we expect the following ROOT files to have been produced
# -rw-rw-r-- 1 tmw microboone  10K Oct 15 17:03 larcv_hist.root
# -rw-rw-r-- 1 tmw microboone 4.3M Oct 15 17:03 larlite_dlmerged.root
# -rw-rw-r-- 1 tmw microboone  16K Oct 15 17:03 larlite_larflow.root
# -rw-rw-r-- 1 tmw microboone 495K Oct 15 17:03 larlite_opreco.root
# -rw-rw-r-- 1 tmw microboone 2.9M Oct 15 17:03 larlite_reco2d.root
# -rw-rw-r-- 1 tmw microboone  23M Oct 15 17:03 merged_dlreco_b9cb82f3-e42a-4dde-ac96-ad4c6d10f466.root
# -rw-rw-r-- 1 tmw microboone  13M Oct 15 17:03 out_larcv_test.root
# -rw-rw-r-- 1 tmw microboone 9.9M Oct 15 17:03 tagger_anaout_larcv.root
# -rw-rw-r-- 1 tmw microboone 1.3M Oct 15 17:03 tagger_anaout_larlite.root
# -rw-rw-r-- 1 tmw microboone  28K Oct 15 17:03 tracker_anaout.root
# -rw-rw-r-- 1 tmw microboone  44K Oct 15 17:03 tracker_reco.root
# -rw-rw-r-- 1 tmw microboone  62K Oct 15 17:03 vertexana.root
# -rw-rw-r-- 1 tmw microboone  14M Oct 15 17:03 vertexout.root

#-rw-r--r-- 1 tmw microboone  10K Oct 24 14:49 larcv_hist.root
#-rw-r--r-- 1 tmw microboone  11K Oct 24 14:37 lcv_trash.root
#-rw-r--r-- 1 tmw microboone 5.4M Oct 24 14:37 merged_dlreco.root
#-rw-r--r-- 1 tmw microboone  72K Oct 24 14:36 nueid_ana_0.root
#-rw-r--r-- 1 tmw microboone 190K Oct 24 14:36 nueid_lcv_out_0.root
#-rw-r--r-- 1 tmw microboone  53K Oct 24 14:36 nueid_ll_out_0.root
#-rw-r--r-- 1 tmw microboone  13M Oct 24 14:49 out_larsoft_copy.root
#-rw-r--r-- 1 tmw microboone  14M Oct 24 14:33 out_larsoft.root
#-rw-r--r-- 1 tmw microboone 173K Oct 24 14:37 shower_reco_out_0.root

rm -f larlite_dlmerged.root larlite_larflow.root larlite_opreco.root larlite_reco2d.root out_larcv_test.root 
rm -f tagger_anaout_larcv.root tagger_anaout_larlite.root tracker_anaout.root tracker_reco.root vertexana.root vertexout.root
rm -f shower_reco_out_0.root nueid_lcv_out_0.root nueid_ll_out_0.root lcv_trash.root nueid_ana_0.root

