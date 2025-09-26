Input:
sampling      # Directory
├─mtd.inp       # mtd job file
├─qsub.pbs     # pbs script
├─control.sh
└─[number]     # Sub-directory of sampling, eg: 1, 2, ...
    └─stru.xyz    # File of Sub-directory


Step 1:
    Rename stru.xyz as $i.

    control.sh:
    ```
    for((i=1;i<=100;i++))
    do
        if [ -d "$i/" ]; then
            
            cd $i
            pwd	
            
            # Step 1
            mv stru.xyz $i.xyz
            
            cd ..
        fi
    done
    ```

    Output:
    sampling      # Directory
    └─[number]     # Sub-directory of sampling, number: 1, 2, ...
        └─[number].xyz    # File of Sub-directory, number: 1, 2, ...


Step 2:
    Send md.inp and qsub.pbs into every sub-directory. Then run the job file.

    control.sh:
    ```
    for((i=1;i<=100;i++))
    do
        if [ -d "$i/" ]; then
            
            cd $i
            pwd	
            
            cp ../md.inp md.inp
            cp ../qsub.pbs qsub.pbs
            qsub qsub.pbs
            sleep 3
            
            cd ..
        fi
    done
    ```

    Output:
    sampling      # Directory
    └─[number]     # Sub-directory of sampling, number: 1, 2, ...
        ├─[number].xyz    # File of Sub-directory, number: 1, 2, ...
        ├─xtb.trj
        └─other temp file.    