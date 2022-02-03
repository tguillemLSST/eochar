#!/usr/bin/zsh

#working
#export input_path_batch=/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/data_PTC
#not working: no plots!
export input_path_batch=/sps/lsst/groups/FocalPlane/SLAC/run5/butler/13141
export output_path_batch=/sps/lsst/users/tguillem/web/batch/debug4/PTC/

#one particular set
##runs=('13141')
##rafts=('R11' 'R14')
##ccds=('S00' 'S22')

#full focal plane
runs=('13141') 
rafts=('R00' 'R01' 'R02' 'R03' 'R10' 'R11' 'R12' 'R13' 'R14' 'R20' 'R21' 'R22' 'R23' 'R24' 'R30' 'R31' 'R32' 'R33' 'R34' 'R41' 'R42' 'R43' 'R04' 'R40' 'R44')
ccds=('S00' 'S01' 'S02' 'S10' 'S11' 'S12' 'S20' 'S21' 'S22')

for run in "${runs[@]}" ;do
    for raft in "${rafts[@]}" ;do
	for ccd in "${ccds[@]}" ;do
	    #echo "---Run: $run"
	    #echo "Raft: $raft"
	    #echo "CCD: $ccd"
	    qsub -P P_lsst -e /sps/lsst/users/tguillem/batch_logs -v run=${run} -v raft=${raft} -v ccd=${ccd} -v input_path=${input_path_batch} -v output_path=${output_path_batch} run_batch.sh
	done
    done	
done

#qsub -P P_lsst -v run=${run_batch} -v raft='R14' -v ccd='S22' -v input_path=${input_path_batch} -v output_path=${output_path_batch} run_batch.sh
