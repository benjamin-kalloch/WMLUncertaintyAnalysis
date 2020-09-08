#!/bin/bash

source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc

TDCSSolver -case $1 | tee -a $1/solver.log
