#!/bin/bash

#### INPUTS ####
SURFACE_FILE=$1
FOAM_CASE=$2
OUTFILE=$3

#### CONSTANTS #####
INTERP_SURFACE_PLUGIN=/home/benny/OpenFOAM/benny-5.0/additional_stuff/PVInterpolateToSurface/build/libInterpolateToSurface.so
FOAM_READER_PLUGIN=/home/benny/OpenFOAM/OpenFOAM-7/platforms/linux64GccDPInt32Opt/lib/paraview-5.6/libPVFoamReader_SM.so
LIBVTKPQCORE=/home/benny/OpenFOAM/ThirdParty-7/build/linux64Gcc/ParaView-5.6.0/lib/libvtkpqCore-pv5.6.so.1
PVPYTHON=/home/benny/OpenFOAM/ThirdParty-7/platforms/linux64Gcc/ParaView-5.6.0/bin/pvpython
TEMP_FILE=/tmp/script$(date +"%N").py

#### PVPYTHON SCRIPT ####
PVPYTHON_SCRIPT_TEMPLATE="from paraview.simple import *;
LoadPlugin('_INTERP_SURFACE_PLUGIN_', remote=False, ns=globals());
LoadPlugin('_FOAM_READER_PLUGIN_', remote=False, ns=globals());
midlayer_surface = STLReader(FileNames=['_SURFACE_FILE_']);
openfoam_case = PVFoamReader(FileName='_FOAM_CASE_');
openfoam_case.IncludeZones = 1;
openfoam_case.MeshParts = ['cellZone_5 - cellZone', 'cellZone_6 - cellZone'];
openfoam_case.VolumeFields = ['magE'];
UpdatePipeline( time=max(openfoam_case.TimestepValues) );
mergeBlocks1 = MergeBlocks(Input=openfoam_case);
print 'Interpolating array: ' + mergeBlocks1.CellData.GetArray(0).GetName()
interpolatetosurface1 = Interpolatetosurface(SurfaceInput=midlayer_surface, FieldData=mergeBlocks1);
SaveData('_OUTFILE_', proxy=interpolatetosurface1, Precision=9);"

#### MAIN ####
# check validity of input files
if [ ! -e $SURFACE_FILE ]; then
    echo "Surface file does not exist! Aborting..."
    exit 1
fi

if [ ! -e $FOAM_CASE ]; then
    echo "FOAM case directory does not exist! Aborting..."
    exit 1
fi

if [ ! -e $(dirname $OUTFILE) ]; then
    echo "Cannot acces output directory! Aborting..."
    exit 1
fi


# setup pvpython script
PVPYTHON_SCRIPT=${PVPYTHON_SCRIPT_TEMPLATE//_INTERP_SURFACE_PLUGIN_/$INTERP_SURFACE_PLUGIN}
PVPYTHON_SCRIPT=${PVPYTHON_SCRIPT//_FOAM_READER_PLUGIN_/$FOAM_READER_PLUGIN}
PVPYTHON_SCRIPT=${PVPYTHON_SCRIPT//_FOAM_CASE_/$FOAM_CASE}
PVPYTHON_SCRIPT=${PVPYTHON_SCRIPT//_SURFACE_FILE_/$SURFACE_FILE}
PVPYTHON_SCRIPT=${PVPYTHON_SCRIPT//_OUTFILE_/$OUTFILE}

# setup envirnment necessary to execute the script
echo "$PVPYTHON_SCRIPT" > $TEMP_FILE

set -- # unset positional commandline arguments as they would interfere with 
       # the OpenFOAM initialization-script

source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc
export LD_PRELOAD=$LIBVTKPQCORE


# execute and cleanup
$PVPYTHON $TEMP_FILE

rm -f $TMP_FILE
