#!/bin/sh

#  spectrassembler_pipeline.sh
# Example of shell script to run full OLC pipeline,
# with minimap to compute the overlaps,
# python script spectral_layout_from_minimap to compute layout,
# and POA to compute consensus sequences.
#
#  Created by Antoine Recanati on 23/09/2016.
#

# DEFINE PATHS. PUT YOUR OWN HERE IF NEEDED.
CURR_DIR=`pwd`
SRC_DIR=`pwd` # DIRECTORY WITH THE SOURCE FILES (spectrassembler). By default we suppose we are in there.
WRK_DIR=`pwd` # PREFERABLE TO PUT another DIRECTORY WHERE TO WORK HERE (unless you are not in SRC_DIR).
POA=`which poa` # OR, IF poa IS NOT IN YOUR PATH, SPECIFY IT MANUALLY
MINIMAP=`which minimap` # OR, IF minimap IS NOT IN YOUR PATH, SPECIFY IT MANUALLY

N_PROCS=`getconf _NPROCESSORS_ONLN`
NUM_NODES=$((N_PROCS - 1)) #SET THE NUMBER OF NODES TO COMPUTE CONSENSUS IN PARALLEL WITH POA
#OPTIONS. MODIFY ACCORDING TO YOUR PREFERENCES
DO_CONSENSUS=true
CLEAN_TEMP_FILES=false
W_LEN=2500
W_OVL_LEN=1250
SIM_THR=850
LEN_THR=3500
VERB_OPT="-vv"

# READS (FASTA FILE)
if [$# != 1]; then
  echo "PLEASE PASS ARGUMENT TO THE SCRIPT :
  PATH/TO/FASTAFILE.FASTA (file containing the raw reads)."
  exit
fi
READS_FILE=$1

# Sanity checks
if [${#POA} == 0]; then
  echo "poa not in PATH. Please specify manually the PATH
  to the executable poa inside the shell script (l.15, instead of POA=`which poa`)!"
fi
if[${#MINIMAP} == 0]; then
  echo "minimap not in PATH. Please specify manually the PATH
  to the executable minimap inside the shell script (l.16, instead of POA=`which poa`)!"
fi
# option to be passed to python script
WRITE_POA_FILES_OPT=""
if ["$DO_CONSENSUS" = true]; then
  WRITE_POA_FILES_OPT="-w"
fi


# Compute alignments with minimap
minimap -S $READS_FILE $READS_FILE > overlaps.mini

# Compute consensus
python $SRC_DIR/spectral_layout_from_minimap.py -f $READS_FILE -m overlaps.mini -r $WRK_DIR $WRITE_POA_FILES_OPT $VERB_OPT --w_len $W_LEN --w_ovl_len $W_OVL_LEN --sim_thr $SIM_THR --len_thr $LEN_THR

# COMPUTE CONSENSUS sequences
if ["$DO_CONSENSUS" == true]; then
  NUM_CC=`find $WRK_DIR/cc_* -type d | wc -l`
  SCORE_MAT="$SRC_DIR/poa-score.mat"
  for cc_idx in $(seq 0 $NUM_CC);
  do
    cd $WRK_DIR/cc_$cc_idx
    for file in poa_in_cc_*.fasta
    do
      while [`jobs | wc -l` -ge $NUM_NODES ]
      do
        sleep 2
      done
      if [! test -f $file.clustal.out] #Do not recompute the same thing twice in case you stopped a computation earlier.
      then
        poa -read_fasta $file -clustal $file.clustal.out -hb $SCORE_MAT &
      fi
    done
    cd ../
    python $SRC_DIR/gen_cons_from_poa.py -cc $cc_idx --poa_mat_path $SCORE_MAT $VERB_OPT
  done
fi

# Clean up if asked
if ["$CLEAN_TEMP_FILES" == true]; then
  for cc_idx in $(seq 0 $NUM_CC);
  do
    rm -r $WRK_DIR/cc_$cc_idx
  done
  echo "MAIN : Cleaned up temporary (poa input and output) files"
fi

# Get back to directory where job was launched
cd $CURR_DIR







#
