#!/bin/bash 

cd /sps/lsst/groups/FocalPlane/SLAC/run5/butler

echo "LSST setup"
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2022_01/loadLSST.zsh 
setup lsst_distrib 

echo "Starting butler ingestion"
mkdir 13141
ingestImages.py 13141 @13141_list.txt --output 13141

echo "DONE"
