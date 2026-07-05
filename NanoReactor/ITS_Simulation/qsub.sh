#PBS -N mtd 
#PBS -j oe
#PBS -q amd
#PBS -l nodes=1:ppn=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_STACKSIZE=1000m
ulimit -s unlimited

user="lijie"

source /share/home/lijie/anaconda3/etc/profile.d/conda.sh
conda activate xtb-py
export PYTHONPATH=/share/home/lijie/pro_test/ITS_Simulation/:$PYTHONPATH

CURR=$PBS_O_WORKDIR
WORK_DIR=$CURR
TMP_DIR="/state/partition1/scratch/$user/$PBS_JOBID"

cd $PBS_O_WORKDIR
mkdir -p  $TMP_DIR
cp -r $WORK_DIR/* $TMP_DIR
cd $TMP_DIR

pwd > log
#/share/home/lijie/python-debug/DynReacExtr/src/dynReacExtr.py -i md.xyz --ts --xtbopt --refine >> log
# dynReacExtr.py -i md.xyz --ts --xtbopt --refine >> log
# python2 trjSplit.py
python3 run.py

cp -rf $TMP_DIR/*   $WORK_DIR
rm -rf $TMP_DIR
