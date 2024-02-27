#!/bin/sh

# Run from this directory
cd ${0%/*} || exit 1

# source executable
FINCH_DIR=`pwd`/../..
application=$FINCH_DIR/build/install/bin/finch

# run application
mpirun -np 1 $application -i inputs_small.json

# create solidification data file
postprocess=$FINCH_DIR/utilities/combine_solidification_data.sh
$postprocess -i single_line_small -o single_line_small_solidification.csv
