#!/bin/bash

RE=1300
PR=7
RI=0.035
S=1
ASPECT=2
LEN=14
RES=192
DT=0.01
SIM_TIME=85

FILEPATH="$(printf "%.0f_aspect%.0d_re%.0e_pr%.0e_Ri%.0e_S%.0d" $RES $ASPECT $RE $PR $RI $S)"
FILEPATH="${FILEPATH//+}"

mpiexec -n 3 --map-by core python3 K-H-Instability.py --Re=$RE --Pr=$PR --Ri=$RI --S=$S --aspect=$ASPECT --len=$LEN \
--res=$RES --dt=$DT --sim_time=$SIM_TIME --basepath=$PWD/$FILEPATH
python3 -m dedalus merge_procs ${FILEPATH}_snapshots --cleanup
python3 interpolate_img_sets.py ${FILEPATH}_snapshots/*.h5 --output=${FILEPATH}_frames
python3 make_movie.py ${FILEPATH}_frames --name=${FILEPATH}_movie.mp4
