from pathlib import Path

comp_config = {
    "nproc": 8, 
    "core": 4000
}
qm_config = {
    "g16_base": "SVP",
    "orca_base": "DEF2-SVP",
    "method": "B3LYP"
}
disp_config = {
    "g16": "em=GD3BJ",
    "orca": "GD3"
}
opt_config = {
    "maxcyc": 100,
    "maxstep": 10,
    "g16_addit1": "cartesian,",
    "g16_addit2": "tight,"
}
ts_config = {
    "recalc": 5,
    "maxcyc": 200,
    "g16_addit": "noeigen,calcfc,"
}
irc_config = {
    "stepsize": 10,
    "maxpoints": 300,
    "g16_addit": "noeigen,calcfc,LQA,"
}
loc_config = {
    'ini_read_data': Path('.ini_data.npy'),
    "reaction_smi_folder": Path('smi/reaction/'),
    "network_smi_folder": Path('smi/network/'),
    "job_file": Path('_job.csv'),
    "relation_file": Path('_relation.csv'),
    "ref_relation_file": Path('_ref_relation.csv'),
    "reaction_pic": Path('_reaction.png'),
    "network_pic": Path('_network.png'),
    "reaction_network_pic": Path('_reaction_network.png')
}


default_setting = {
    "comp_config": comp_config, 
    "qm_config": qm_config,
    "disp_config": disp_config,
    "opt_config": opt_config,
    "ts_config": ts_config,
    "irc_config": irc_config,
    "loc_config": loc_config
}
