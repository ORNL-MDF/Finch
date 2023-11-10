#!/bin/sh

# Run from this directory
cd ${0%/*} || exit 1

# source executable
FINCH_DIR=`pwd`/../..
application=$FINCH_DIR/build/install/bin/finch

# run application
mpirun -np 8 $application -i inputs.yaml

# create solidification data file
postprocess=$FINCH_DIR/utilities/combine_solidification_data.sh
$postprocess -i single_line -o single_line_solidification.csv
