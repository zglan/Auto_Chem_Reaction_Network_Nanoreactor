#PBS -N mtd 
#PBS -j oe
#PBS -q amd
#PBS -l nodes=1:ppn=8

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OMP_STACKSIZE=1000m
ulimit -s unlimited

user="lijie"

CURR=$PBS_O_WORKDIR
WORK_DIR=$CURR
TMP_DIR="/state/partition1/scratch/$user/$PBS_JOBID"

cd $PBS_O_WORKDIR
mkdir -p  $TMP_DIR
cp -r $WORK_DIR/* $TMP_DIR 
cd $TMP_DIR

pwd > log
dynReacExtr.py -i xtb.trj --ts --xtbopt --refine >> log

cp -rf $TMP_DIR/*   $WORK_DIR
rm -rf $TMP_DIR
