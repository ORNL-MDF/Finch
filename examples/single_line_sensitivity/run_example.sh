#!/bin/sh

# Run from this directory
cd ${0%/*} || exit 1

# source executables
FINCH_DIR=`pwd`/../..
application=$FINCH_DIR/build/install/bin/finch
sensitivity=$FINCH_DIR/build/install/bin/finch_sensitivity

# reference double solve, writing the temperature field and solidification data
echo "### reference double solve"
$application -i inputs.json

# one OTI solve gives every parameter derivative, checked against central finite
# differences. Requires Sparrow at configure time; without it the
# finch_sensitivity target is not built and this example cannot run.
echo "### sensitivity, default finite difference step"
$sensitivity -i inputs.json

# The default step straddles the latent heat branch, so density, specific_heat
# and absorption disagree by ~28% above. The derivative is not wrong: a smaller
# step recovers agreement, which is the opposite of how finite difference error
# normally behaves and is the signature of a discontinuity rather than a bug.
echo "### sensitivity, refined finite difference step"
FINCH_FD_STEP=1e-8 $sensitivity -i inputs.json

# Removing the latent heat branch entirely makes every parameter agree at the
# default step, confirming the diagnosis above.
echo "### sensitivity, latent heat removed"
$sensitivity -i inputs_nolatent.json

# Same case and same input file on four ranks: ranks_per_dim is ignored when it
# does not match the communicator size, leaving the decomposition to
# MPI_Dims_create, which gives 2x2x1 here. T_sum and its derivatives match the
# single rank values; T_probe does not, since the probe is the centre of the
# owned index space and so is a different node under domain decomposition.
echo "### sensitivity, four ranks"
mpirun -np 4 $sensitivity -i inputs.json
