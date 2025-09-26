#!/usr/bin/env python3

# modify part of ' xtb --md --input mtd.inp *.xyz > log'
pbsFile = '''
#PBS -N mtd-%s
#PBS -l walltime=999:00:00
#PBS -j oe
#PBS -l nodes=1:ppn=8
#PBS -V

user="pengt"
istart=%s
iend=%s
WORK_DIR=$PBS_O_WORKDIR
TMP_DIR="/state/partition1/scratch/$user/$PBS_JOBID"

cd $PBS_O_WORKDIR
mkdir -p $TMP_DIR

for ((i=${istart}; i<=${iend}; i++))
do
  cp -r $WORK_DIR/${i} $TMP_DIR/
done

cd $TMP_DIR

for ((i=${istart}; i<=${iend}; i++))
do
  cd ${i}

  xtb --md --input mtd.inp *.xyz > log
  
  # g16 *.gjf
  
  cp -rf $TMP_DIR/${i} $WORK_DIR/
  cd ..
done

rm -rf $TMP_DIR
'''
start = 201
tot = 100
#num = 50
step = 10 #int(tot / num)
l = range(start, start + tot)
groups = [l[i: i + step] for i in range(0, len(l), step)]

for idx, group in enumerate(groups):
    with open('%s.pbs'%(idx + 1), 'w') as pbs:
            pbs.write(
                pbsFile % (idx + 1, group[0], group[-1]))
