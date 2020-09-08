#!/bin/bash

#
# This script automates the setup of the OpenFOAM case for usage
# in the uncertainty analysis.
# The mesh is converted into the OpenFOAM format and the midlayer
# is copied into the case directory.
#
# Note: Specify the path (SAMPLE_DIR) to a template/dummy testcase
#       that is used as a skelleton for the new case directory.
#

MESHFILE=$1 # teatrahedral volume head mesh in GMSH ASCII file format
MIDLAYER=$2 # surface file representing the midlayer surface in STL file format

if [[ -z "$MESHFILE" ]]; then
    echo "Parameter 1 (meshfile) missing!"
elif [[ -z "$MIDLAYER" ]]; then
    echo "Parameter 2 (midlayer) missing!"
fi

SAMPLE_DIR=/path/to/dummy/OF_case

# setup case dir structure
cp -r ${SAMPLE_DIR}/system .
cp -r ${SAMPLE_DIR}/0_orig .
cp ${SAMPLE_DIR}/Field_template .
cp ${SAMPLE_DIR}/runSim.sh .
cp ${SAMPLE_DIR}/interpolate_to_midlayer_and_retrieve_values.sh .
mkdir constant

sleep 10

# setup mesh
gmshToFoam $MESHFILE
sleep 10
checkMesh -allTopology -allGeometry | tee checkmesh.log
sleep 5
transformPoints -scale '(0.001 0.001 0.001)'
surfaceTransformPoints -scale '(0.001 0.001 0.001)' $MIDLAYER midlayer.stl
renumberMesh -overwrite | tee renumbermesh.log
sleep 10

# setup fields
cp -r 0_orig 0
setFields
sleep 5
sed -i 's/1111/%LESION_VAL%/g' 0/sigma
sleep 5
sed -i 's/2222/%CSF_VAL%/g' 0/sigma
sleep 5
sed -i 's/3333/%AIR_VAL%/g' 0/sigma
sleep 5
sed -i 's/4444/%SKIN_VAL%/g' 0/sigma
sleep 5
sed -i 's/5555/%SKULL_VAL%/g' 0/sigma
sleep 5
sed -i 's/6666/%CSF_VAL%/g' 0/sigma
sleep 5
sed -i 's/7777/%WM_VAL%/g' 0/sigma
sleep 5
sed -i 's/8888/%ELECTRODE_VAL%/g' 0/sigma
sleep 5
sed -i 's/9999/%GM_VAL%/g' 0/sigma
sleep 5
mv 0 0_clean

sleep 60

# create process-directories
CURRENT_DIR=$(pwd)
CASE_NAME=$(basename $CURRENT_DIR)
for n in {0,1,2,3,4,5}; do
    cp -r $CURRENT_DIR $(dirname $CURRENT_DIR)/${CASE_NAME}_$n
done

echo "Case setup done."