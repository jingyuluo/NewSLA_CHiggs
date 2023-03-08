#!/bin/bash

variable=${1}
year=${2}
lepton=${3}
nhot=${4}
nT=${5}
nW=${6}
nB=${7}
nJ=${8}
exeDir=${9}
subDir=${10}
condorDir=$PWD

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd $exeDir
eval `scramv1 runtime -sh`

python -u hists.py \
  -v $variable \
  -y $year \
  -l $lepton \
  -nh $nhot \
  -nt $nT \
  -nw $nW \
  -nb $nB \
  -nj $nJ \
  -sd $subDir
