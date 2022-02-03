#!/usr/bin/zsh

cd /sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/eochar/notebooks/CTE

#WARNING: LSST must not be already set up!
echo "LSST load"
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2022_01/loadLSST.zsh 
echo "LSST setup"
setup lsst_distrib 

echo "Starting to run the eochar job"
source export.sh
echo "Job configuration:"
echo ${run}
echo ${raft}
echo ${ccd}
echo ${input_path}
echo ${output_path}
python CTE_diagnostic_BOT.py ${run} ${raft} ${ccd} ${input_path} ${output_path}
echo "DONE"

#####debug lines
#python CTE_diagnostic_BOT.py 13141 R14 S22 /sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_PTC /sps/lsst/users/tguillem/web/PTC/
#python CTE_diagnostic_BOT.py 13141 R14 S22 /sps/lsst/groups/FocalPlane/SLAC/run5/butler/13141 /sps/lsst/users/tguillem/web/batch/PTC/
