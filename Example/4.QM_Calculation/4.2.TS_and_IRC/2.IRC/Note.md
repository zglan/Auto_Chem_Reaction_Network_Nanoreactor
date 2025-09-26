reactPathRefine.py is needed.(Optional)

Input:
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  ├─reactPathRefine.py
  └─[1]        # Sub-directory of ts
    └─[numbers]         # Sub-directory of directory 1
      ├─[id]        
      | └─[ids].gjf
      └─info

Output:
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  ├─reactPathRefine.py
  └─[1]        # Sub-directory of ts
    ├─irc_workspace
    |  ├─[numbers]
    |  | ├─error         # Error information of gaussian.
    |  | ├─irc.gjf
    |  | └─irc.log
    |  └─irc-[n].pbs     
    ├─irc       # Include successful finish log file and gjf file.
    ├─irc_temp
    ├─job_list          # PBS job IDs.
    ├─irc_job.csv        
    └─irc_state.csv      # Include pbs job execution path, pbs job running status, pbs job IDs