loc_config = {
    "orca_loc": '/data/home/pengt/orca/orca'
}
user_config = {
    'user': 'pengt'
}
server_config = {
    "nproc": 8,
    "sleep_t": 300, # seconds
    "max_nproc": 120,
    "qsub_num": 20
}
comp_config = {
    "nproc": 8, 
    "core": 4000,
}
qm_config = {
    "g16_base": "SVP",
    "orca_base": "DEF2-SVP",
    "method": "UB3LYP"
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

pbs_ini = f'''
#PBS -N %s-%s
#PBS -l walltime=999:00:00
#PBS -j oe
#PBS -l nodes=1:ppn=%s
#PBS -V

user="{user_config["user"]}"

WORK_DIR=$PBS_O_WORKDIR/%s
TMP_DIR="/state/partition1/scratch/$user/$PBS_JOBID"

cd $PBS_O_WORKDIR
mkdir -p $TMP_DIR
'''

pbs_job = '''
istart=%s
iend=%s

for ((i=${istart}; i<=${iend}; i++))
do
  cp -r $WORK_DIR/${i} $TMP_DIR/${i}
done

cd $TMP_DIR

for ((i=${istart}; i<=${iend}; i++))
do
  cd ${i}
  %s
  cd ..
  cp -rf $TMP_DIR/${i} $WORK_DIR/
done

rm -rf $TMP_DIR
'''

default_setting = {
    "server_config": server_config,
    "comp_config": comp_config, 
    "qm_config": qm_config,
    "disp_config": disp_config,
    "opt_config": opt_config,
    "ts_config": ts_config,
    "irc_config": irc_config,
    "loc_config": loc_config,
    "pbs_content": pbs_ini + pbs_job,
}
