#!/bin/sh

# Run from this directory
cd ${0%/*} || exit 1

# source executable
FINCH_DIR=`pwd`/../..
application=$FINCH_DIR/build/install/bin/create_scan_paths

# run application
$application -i inputs.json
