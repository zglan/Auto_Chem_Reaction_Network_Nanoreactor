reactPathRefine.py is needed.(Optional)

Input:
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  
  └─[1]        # Sub-directory of ts
    ├─reactPathRefine.py
    └─[numbers]         # Sub-directory of directory 1
      ├─[id]        
      | └─[ids].gjf
      └─info

Usage(For TS Optimization):

  Modify the end of the program and call `G16_TS`:

  ```python
  if __name__ == '__main__':
      G16_TS(10,5).main()
  ```
  Then:
  ```
  reactPathRefine.py
  ```

Output:
eg: 
```
 Moving gaussion_gjf files: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 70/70 [00:00<00:00, 427.74it/s]
[0/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:36:44 2025
[0/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:38:45 2025
[3/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:40:45 2025
[6/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:42:45 2025
[13/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:44:46 2025
[14/14]  Start Time: Thu Mar  6 16:36:44 2025; Current Time: Thu Mar  6 16:46:46 2025       # [Completed / Total]
  Moving gaussion_log files: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 70/70 [00:00<00:00, 250.08it/s]
  Loading reaction-1:   0%|                                                                                                                                                                                | 0/70 [00:00<?, ?it/s]
  Loading reaction-1: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 70/70 [00:00<00:00, 92.08it/s]
```

sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  
  └─[1]        # Sub-directory of ts
    ├─reactPathRefine.py
    ├─[numbers]         # Sub-directory of directory 1
    | ├─[id]     
    | |  └─[ids].log
    | └─workspace       # Sub-directory include log file which successfully complete ts search.
    |     └─[ts_search success].log     
    ├─ts_workspace
    |  ├─[numbers]
    |  | ├─error         # Error information of gaussian.
    |  | ├─ts.gjf
    |  | └─ts.log
    |  └─ts-[n].pbs     
    ├─job_list          # PBS job IDs.
    ├─ts_job.csv        
    └─ts_state.csv      # Include pbs job execution path, pbs job running status, pbs job IDs

Note: If there are too many initial guess file needed to send, you can create a temperory directory such as temp-1.