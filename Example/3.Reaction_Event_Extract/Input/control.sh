#!/usr/bin/bash

for((i=293;i<=300;i++))
do
	if [ -d "$i/" ]; then
		
		cd $i
		# # Step 1
		# mv ../$i.xyz ./
		# mv stru_xyz $i.xyz
		
		pwd	
		# # Step 2
		# cp ../md.inp md.inp
		# cp ../qsub.pbs qsub.pbs
		# qsub qsub.pbs
		# sleep 15
		# xtb --md --input md.inp *.xyz > log

		# # Step 3
		pwd >> log
		cp ../qsub.sh .
		qsub qsub.sh
		# dynReacExtr.py -i xtb.trj --refine --ts --xtbopt >> log
		
		cd ..
	fi
done

