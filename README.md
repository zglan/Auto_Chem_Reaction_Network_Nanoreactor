<<<<<<< HEAD
# Auto_Chem_Reaction_Network_Nanoreactor
=======
<h1 align="center"> Automated Exploration of Reaction Network and Mechanism



## Summary

This software is developed in the group of Prof. Zhenggang Lan at South China Normal University, Guangzhou, China. The purpose of this code is the automatic exploration of reaction network and mechanism. The introduction of this software is given in the below reference. Please cite this work if you wish to use the current code.

***Cite this: [Yutai Zhang](https://orcid.org/0009-0001-2327-0623), [Chao Xu](https://orcid.org/0000-0002-4043-2954), [Zhenggang Lan](https://orcid.org/0000-0002-8509-0388), Automated Exploration of Reaction Networks and Mechanisms Based on Metadynamics Nanoreactor Simulations. [J. Chem. Theory Comput. 2023, 19, 23, 8718-8731.](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752)***


The working procedure is given in the below figure.
<div style="text-align: center">

![Image](https://pubs.acs.org/cms/10.1021/acs.jctc.3c00752/asset/images/medium/ct3c00752_0001.gif)
</div>


Next, based on the flowchart in the above figure, we will introduce a step-by-step tutorial on how this software is used in practical calculations. Before entering the specific process, let's briefly introduce the environment, libraries, and packages that this software depends on.

## 1. Initial Settings in Linux:

All codes are written with Python. It is recommended to use the whole code in the operator system of Linux. We did not test these scripts in Windows.

### 1.1 Required Packages:

Several packages are needed for the installation/use of the current codes, which includes 
[Anaconda](https://www.anaconda.com/download/), [xtb](https://xtb-docs.readthedocs.io/en/latest/setup.html), [packmol](https://m3g.github.io/packmol/download.shtml),[VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) and Gaussian 16.

* `Anaconda`:  It is better to set up the Python enrivonments in Anaconda for convenience. 
* `xtb`: This software is used for the metadynamics simulation by the semiempirical DFTB quantum mechanical method GFNn-xTB.
* `Packmol`: This software is used to generate the initial configurations of molecular systems.
* `VMD`: This software is used to view molecular structure and dynamic trajectory.
* `Gaussian 16`: This software is used to perform the TS optimization and IRC.   
Please notice that we also support ORCA code for TS optimizations. If you wish to use ORCA, please contact with us.


### 1.2 Python Library:

After Anaconda is installed, the python running environment is created by:
```commandline
conda create -n my_env python=3.10
```

Then  the Python environment is activated by: 
```commandline
conda activate my_env
```

Serval Python Libraries should be installed:

#### 1.2.1 `Rdkit`:  

  Recommend using the following command for installation
```commandline
conda install -c conda-forge rdkit
```
#### 1.2.2 `Openbabel`:  

  Recommend using the following command for installation
```commandline
conda install openbabel -c conda-forge
```
#### 1.2.3 `Hmmlearn`:  

  The library hmmlearn includes a set of algorithms for unsupervised learning and the inference to Hidden Markov Models. Recommend using the following command for installation
```commandline
conda install conda-forge::hmmlearn=0.2.7
```
#### 1.2.4 `Scikit-learn`:  

  The library Scikit-learn is an open-source Python library for machine learning algorithms. Recommend using the following command for installation
```commandline
conda install anaconda::scikit-learn
```
#### 1.2.5 `NetworkX`:  

  NetworkX is a Python package for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks. Recommend using the following command for installation
```commandline
conda install anaconda::networkx
```
#### 1.2.6 `Tqdm`:  

  The library tqdm is mainly used to provide the labels for Python progress, which can add a progress information in Python long loops. Recommend using the following command for installation

```commandline
conda install conda-forge::tqdm
```

#### 1.2.7 `Matplotlib`:

```commandline
conda install conda-forge::matplotlib
```

#### 1.2.8 `Cmocean`: 

  This library is used to create graphs with the good view

```commandline
 conda install conda-forge::cmocean=3.0.3
```

### 1.3 .bashrc Settings:
It is recommended to add the absolute path of software into the `.bashrc` file. eg: `packmol`:
```commandline
vi ~/.bashrc
PATH=$PATH:/path/to/packmol
source ~/.bashrc
```

And add the absolute path of the following folders and python program into the `.bashrc` file：which include `DynReacExtr`、`ReacNetDraw`、`QmJob`、`sphereMaker.py`、`trj_split.py`.
```
export PATH=.../python:$PATH
export PATH=.../python/QmJob:$PATH
export PATH=.../python/ReacNetDraw:$PATH
export PATH=.../python/DynReacExtr/src:$PATH
export PATH=.../python/sphereMaker.py:$PATH
export PATH=.../python/trj_split.py:$PATH
```

Notes:
1. Using the latest version of `hmmlearn` and `Cmocean` may result in errors. If these problems happen, try to switch to the old versions, i.e. `hmmlearn=0.2.7` and `cmocean=3.0.3`.
2. In principle, the code is written based on the Python 3.10, and we did not perform tests for normal operation under other Python versions.


After completing all the above settings, you can start the actual simulation. Next, we will introduce the specific steps of each section based on the flowchart in the summary. Firstly, there is the initial sampling section. 

---

## 2. Initial Sampling:

The purpose of initial sampling is to obtain a series of reasonable initial molecular configurations under set conditions (such as T=300K). In this section, two installed packages, `packmol` and `xtb`, are required.  In addition, the script `sphereMaker.py` is needed.

According to the program operation process, we will divide the initial sampling section into two parts: initial conformation construction and molecular dynamics (MD) sampling, and introduce them in sequence.

### 2.1 Initial Conformation Construction:

In this step, the whole system is generated initially, which may contain many molecules (same and different). All of compounds are packed in a reactor with the sphere shape or others. Normally, the sphere shape is recommended. 

Here the working procedure is given as follows.

1. A few `.xyz` files are needed as the input, and each file contains a single molecule. If different compounds are involved, each `.xyz` file contains a single compound. For a group of identical molecules, only a single `.xyz` file is needed.

2. Run the script of `sphereMaker.py`, which automatically calls the `packmol` and generate the initial conformation file. 

Please specify the numbers of each molecule, as well as the size of the reactor.

Please check the parameters of `sphereMaker.py` by:

```
sphereMaker.py -h
```
```
------------------------------------------------------------
usage: sphereMaker.py [-h] -i [INPUTFILE ...] -n [NUMS ...] 
-d DENSITY --names [NAMES ...] [--center [CENTER ...]] [--dis DIS]

Process some integers.

options:
  -h, --help            show this help message and exit
  -i [INPUTFILE ...], --inputfile [INPUTFILE ...]
                        Input mole file, e.g. water.xyz OH.xyz
  -n [NUMS ...], --nums [NUMS ...]
                        Input mole nums, e.g. 10 1
  -d DENSITY, --density DENSITY
                        Set system density, unit in g/ml, e.g. 10
  --names [NAMES ...]   Set Output file name
  --center [CENTER ...]
                        Set Central Molecules
  --dis DIS             Set Mol distance, unit is Å
------------------------------------------------------------
```
Notes: When you use `--center`, the number of molecules should be ONE. 

Here is an example to illustrate the usage of `sphereMaker.py`. 

### 2.1.1 Example:

Build a system containing two compounds (A and B). The number of molecule A is 1, and the number of molecule B is 2. The distance between molecules needs to be larger than 3 A, and the density of the system is 1.5. Assuming that two molecular species are given in "1. xyz" and "2. xyz". Please run the following comments:

```python
sphereMaker.py -i 1.xyz 2.xyz -n 1 2 -d 1.5 --name test --dis 3
```
This gives the following output.
```python
*** Sphere Radii is 8.30 A.
*** Sphere Radii is 15.68 Bohr.
```
The value of Sphere Radii (in Bohr) is used for setting the sphere of nanoreactor. Please record this value, which will be used for setting the sphere radius of nanoreactor in MTD simulation. (See 3.3 MTD job file section)

After the initial conformation construction is completed, we can proceed to the second step of initial sampling, molecular dynamics (MD) sampling

## 2.2 MD Sampling:

The initial sampling mainly includes two methods: molecular dynamics (MD) sampling and Monte Carlo (MC) sampling. Here the molecular dynamics sampling method is used which is supported by `xtb`.

### 2.2.1 MD job file:

Prepare a job file "md.inp" for the MD run using xtb. An example is given below. Delete comments (the sentence after #) when you use the below file. Please check the [xtb website](https://xtb-docs.readthedocs.io/en/latest/md.html) for more details.

```
$chrg 0
$spin 0
$scc
   maxiterations=800
   temp=1000
$cma
$md
   sccacc=2
   temp=298 
   shake=2
   hmass=1
   time=30
   dump=200
   step=0.5
$wall
   potential=logfermi
   beta=0.5
   temp=6000
   sphere: 10, all # In bohr，Setting it according to the output result of sphereMaker.py.(See 2.1.1)
$end
```
### 2.2.2 Run 

Use the following command to run MD for initial sampling:

```python
xtb --md --input md.inp test.xyz > log
```
The output file is given as follows.
```
         time (ps)    <Epot>      Ekin   <T>   T     Etot
      0    0.00      0.00000   0.0253    0.    0.   -12.21325
est. speed in wall clock h for 100 ps :  0.03
.............
block <Epot> / <T> :     -12.24057  280.     drift:  0.20D-04   Tbath : 280.
  25200   12.60    -12.23947   0.0106  299.  320.   -12.23146
  25400   12.70    -12.23947   0.0096  299.  290.   -12.23186
  25600   12.80    -12.23949   0.0133  299.  401.   -12.23010
  25800   12.90    -12.23951   0.0114  299.  342.   -12.23034
  26000   13.00    -12.23952   0.0096  299.  287.   -12.23216
  26200   13.10    -12.23954   0.0116  299.  349.   -12.23154
  26400   13.20    -12.23956   0.0096  299.  290.   -12.23097
  26600   13.30    -12.23958   0.0093  299.  281.   -12.23126
  26800   13.40    -12.23959   0.0099  299.  297.   -12.23138
  27000   13.50    -12.23960   0.0120  299.  362.   -12.23118
  27200   13.60    -12.23962   0.0083  299.  248.   -12.23311
  27400   13.70    -12.23963   0.0065  299.  194.   -12.23288
  27600   13.80    -12.23962   0.0082  298.  246.   -12.23081
  27800   13.90    -12.23962   0.0139  297.  418.   -12.23192
  28000   14.00    -12.23963   0.0070  298.  211.   -12.23189
  28200   14.10    -12.23964   0.0103  297.  310.   -12.23180
  28400   14.20    -12.23964   0.0087  297.  263.   -12.23035
  28600   14.30    -12.23965   0.0093  297.  278.   -12.23117
  28800   14.40    -12.23966   0.0137  297.  411.   -12.23088
  29000   14.50    -12.23967   0.0093  297.  281.   -12.23102
  29200   14.60    -12.23968   0.0120  297.  362.   -12.23192
  29400   14.70    -12.23969   0.0103  297.  311.   -12.23277
  29600   14.80    -12.23970   0.0079  297.  236.   -12.23293
  29800   14.90    -12.23971   0.0095  297.  285.   -12.23213
```

Then you may obtain many geometries from the MD run, which were used as initial configurations for the next step.

Notes:
1. Check the evolution of `<T>` to make sure that the temperature is balanced.
2. Don't choose the structure in which the reaction takes place.  
3. If too many reactions take place in the MD run, please lower the temperature or run the short time dynamics.
4. For different reactive systems, please be careful on how to run MD, see below.

### 2.2.3 Different Sampling Situation:

* For weak reactivity systems,  long-time MD samplings are needed. This means: total Sampling Time: > 10 ps


* For highly reactive system, eg: free radical, short-time MD samplings should be done. This means: total Sampling Time: < 0.5 ps

Please also be careful on the spin multiplicities. This is a hard problem. Please check the xTB documents for details.

After the MD run, use `trj_split.py` to generate sampling structures that are saved in `stru.xyz`(`\sampling\*\stru.xyz`). Note that this script is implemented using Python2. Please use the command `python2 trj_split.py` if needed.


When the initial sampling is completed, the initial structure of a series of MTD simulations is ready. Next, we will perform a series of MTD simulations with different parameter settings.


## 3. MTD Simulation with Different Settings:

As long as the time for MD simulation is long enough, all information about the chemical transformation from reactants to products (also known as reaction events) can be obtained. However, in actual simulations, this condition is not feasible. Metadynamics (MTD), as an enhanced sampling method, can significantly increase the frequency of reaction events and is currently a common approach in reactive molecular dynamics.

In MTD simulation, the form of bias potential energy is given as the sum of time-dependent Gaussian functions. In `xtb`, the root mean square deviation (RMSD) in Cartesian space is chosen as a measure of collective variables. 

When conducting MTD simulations, we need to run MTD simulations with different parameter settings to determine convergence. Please refer to the [original paper](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752) for more details.

In this step, the structure file(`stru.xyz`) is needed and the package `xtb` is used. If calculations are run with the linux system, a few shell scripts may be useful for the automatic operation of all calculation procedures. We provide some example files for Shell script and PBS files. Please feel free to modify them. 

### 3.1 Shell scripts:

Here is an example of a shell script (control.sh):

```shell
#!/usr/bin/bash

for((i=1;i<=500;i++)) # Number of Folders
do
      ## 1. Construct Basic Sub-folders
      mkdir $i  
      ##
      if [ -d "$i/" ]; then
         cd $i
         pwd

         mv ../$i.xyz ./ # Reactive System
         mv stru.xyz $i.xyz # Weak Reactivity System
         cp ../qsub.pbs qsub.pbs
         cp ../md.inp md.inp
         ##

         ## 2. Submit Job
         # qsub qsub.pbs
         # sleep 5
         ##

         ## 3. Trajectories Analysis       
         # pwd >> ../log
         # dynReacExtr.py -i xtb.trj --refine >> ../log ## For Network Drawing
         ##
         cd ..
   fi
done
```

Please execute the shell script using the following command:
```
bash control.sh
```

### 3.2 PBS Script:

The PBS job script is a small text file that contains information about the resources required for the job, including time, number of nodes, and memory. The PBS script also includes commands for the supercomputer system to execute. Here is an example of PBS script (like qsub.pbs):

```
#PBS -N xxx
#PBS -j oe
#PBS -l nodes=1:ppn=8
#PBS -V
#PBS -l walltime=999:00:00

# Username
user="xxx"

CURR=$PBS_O_WORKDIR
WORK_DIR=$CURR
TMP_DIR="/state/partition1/scratch/$user/$PBS_JOBID"

export OMP_NUM_THREADS=8

cd $PBS_O_WORKDIR
mkdir -p  $TMP_DIR
cp -r $WORK_DIR/* $TMP_DIR 
cd $TMP_DIR

# xtb command
xtb --md --input xxx.inp *.xyz > log

cp -rf $TMP_DIR/*   $WORK_DIR
rm -rf $TMP_DIR
```

### 3.3 MTD job file:

This file contains the running setting for MTD. For detailed parameter settings, please refer to the documentation on the xtb website. 
Here is an example of a mtd running settings file (mtd.inp):
```
$chrg 0
$spin 0
$scc
   maxiterations=800
   temp=1000
$cma
$md
   shake=0
   hmass=1
   sccacc=2
   temp=1500
   time=30
   dump=2
   step=0.5
$metadyn
   save=1
   kpush=0.8
   alp=0.7
$wall
   potential=logfermi
   beta=0.5
   temp=6000
   sphere: 10, all
$end
```

Please notice that the MTD may give too drastic reactions for highly active systems. In this case, the nomral MD (even at low temperature) is enough for initial samplings.  

As a summary, this step is used to generate trajectories containing reaction events (or so called reactive trajectories). Whether using MTD or MD is dependent on the system reactivity.

Next, we will introduce the reaction event identification.

## 4. Reaction Event Identification:

The purpose of this step is to extract reaction events from last MTD/MD trajectories, and identify all of the reactants and products for each reaction event. For detailed parameter settings and principles, please refer to this article: ***J. Chem. Theory Comput. 2023, 19, 23, 8718-8731.***

In this step, the script `dynReacExtr.py` is needed. This script is mainly based on the hmmlearn library and RDkit library to analyze the MTD/MD trajectories after mtd simulation.

### 4.1 Basic Settings of `dynReacExtr.py`:

The main function of `dynReacExtr.py` is to extract important snapshots of reaction events and generate input files for subsequent quantum chemistry calculations based on these structures. Currently, the quantum chemistry calculation software supported by this script are Gaussian16 and ORCA. The default settings for `dynReacExtr.py` are located in `\DynReacExtr\src\_setting_default.py`. When you use it, please make modifications according to your purpose at any time, such as electronic structure methods, basic sets, TS optimization, IRC, and so on.

```python
# The number of CPU(nproc) and Memory(core):
comp_config = {
    "nproc": 8, 
    "core": 4000
}
# Basis and Functional:
qm_config = {
    "g16_base": "SVP",
    "orca_base": "DEF2-SVP",
    "method": "B3LYP"
}
# Dispersion Correction:
disp_config = {
    "g16": "em=GD3BJ",
    "orca": "GD3"
}
# Structure Optimization Settings:
opt_config = {
    "maxcyc": 100,
    "maxstep": 10,
    "g16_addit1": "cartesian,",
    "g16_addit2": "tight,"
}
# TS Optimization Settings:
ts_config = {
    "recalc": 5,
    "maxcyc": 200,
    "g16_addit": "noeigen,calcfc,"
}
# IRC Settings:
irc_config = {
    "stepsize": 10,
    "maxpoints": 300,
    "g16_addit": "noeigen,calcfc,LQA,"
}
# Path Settings:
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
```

Now let us discuss all step one by one.

### 4.2 Basic Usage of `dynReacExtr.py`:

Try to run 
```
dynReacExtr.py -i trj --xtbopt --ts --refine --tspin 1 --tchrg 0
```
This step will automatically extract the geometries from MTD/MD runs and generate the TS optimization files. 

Besides, the information of the preliminary reaction network based on the MD or MTD is also constructed, in which all possible species are written as the SMILES code. Please notice that here some species may appear due to the high energy MD or MTD, which may not come out in the normal reactive conditions.

This script is written for the analysis of the single trajectory. If there are many trajectories, shell scripts can be used to make analysis for each trajectory in the parallel way (eg: previously mentioned scripts `control.sh`).

```
#!/usr/bin/bash

for((i=1;i<=500;i++))
do
    if [ -d "$i/" ]; then
         cd $i
         pwd

         ## 3. Trajectories Analysis       
         pwd >> ../log
         dynReacExtr.py -i xtb.trj --ts --xtbopt >> ../log
         cd ..
   fi
done
```
Then:
```
bash control.sh
```

Regarding the parameter description and usage of `dynReacExtr.py`, please use the following command to read:

```
dynReacExtr.py -h
```

```
---------------------------------------------------------------------------
usage: dynReacExtr.py [-h] [-i INPUTFILE] [--interval INTERVAL] [--scale SCALE]
                      [--mode MODE] [--hidmatr HIDMATR HIDMATR HIDMATR HIDMATR]
                      [--obvmatr OBVMATR OBVMATR OBVMATR OBVMATR] [--tspin TSPIN]
                      [--tchrg TCHRG] [--jobset JOBSET JOBSET] [--ts] [--opt] [--smooth]
                      [--neb] [--xtbpath] [--xtbopt] [--network] [--draw] [--refine]
                      [--load]

options:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE
                        Input trajectory file
  --interval INTERVAL   Set interval
  --scale SCALE         Set frame in detect reaction, default 30 frames
  --mode MODE           Select mode: reac or prod, default extract from product
  --hidmatr HIDMATR HIDMATR HIDMATR HIDMATR
                        Matrix hidden of HMM parameters
  --obvmatr OBVMATR OBVMATR OBVMATR OBVMATR
                        Matrix obvserve of HMM parameters
  --tspin TSPIN         Set spin for entire system
  --tchrg TCHRG         Set charge for entire system
  --jobset JOBSET JOBSET
                        Set job set, default scale is -30~30, interval is 10
  --ts                  Proceed gaussion ts job
  --opt                 Proceed gaussion opt job
  --smooth              Proceed smooth trajectory job to get double end points. !! Need
                        nebterpolator !!
  --neb                 Proceed orca neb-ts job
  --xtbpath             Proceed xtb path job
  --xtbopt              Proceed xtb opt job
  --network             Proceed network analysis job
  --draw                draw reactions
  --refine              refine reactions
  --load                load reaction pre-analysis
```

Please note that the location for running shell script `control.sh` (or Python script `dynReacExtr.py`) is in the folder of "sampling". The specific directory relationship is as follows:
```
sampling      # Directory
├─control.sh(Or dynReacExtr.py)
└─[number]     # Sub-directory of sampling, eg: 1, 2, ...
    └─xtb.trj    # trajectory file of mtd.
```

### 4.3 Output:

During the process of running script `dynReacExtr.py`, some key information will be extracted and generated, including the SMILES code of involved compounds in the reaction events, the number of occurrences of the reactions, and so on, for example:
```
Filter sigal transition: 100%|███████████████████████████████| 19/19 [00:00<00:00, 311.45it/s]
   *  1 times  |  [H]C1C([H])C(N([H])C(O)C([H])([H])[H])C([H])C([H])C1O+[H]O+[H]O[H]->[H]C1C([H])C(NC(O)C([H])([H])[H])C([H])C([H])C1O+[H]O[H]+[H]O[H]
```
Here we show that a reaction represented by the SMILES codes appears only once.

After running script `dynReacExtr.py`, a series of files and folders will be generated, with the specific directory relationships as follows:

```
[number]     # Sub-directory of sampling, eg: 1, 2, ...
├─sp            # sp directory, for single point calculation. The command "xtb [id].xyz --input [id].inp" have been used.
|  └─[number]   # Sub-directory of sp
|    ├─[id].charges     # charges data
|    ├─[id].xyz         # xyz file
|    ├─[id].inp         # job file
|    ├─[id].log         
|    └─[id].wbo
├─ts            # ts directory, for ts optimization by gaussian.
|  └─[number]   # Sub-directory of ts
|    └─[id].gjf         # TS initial guess structure file.

<!-->
├─xtbopt        # xtbopt directory, for optimization by xtb. The command "xtb [id].xyz --input [id].inp --opt" have been used.
|  └─[number]   # Sub-directory of xtbopt
|    ├─[id]-opt.xyz     # xyz file after optimization.
|    ├─[id].inp         # Optimization Setting file.
|    ├─[id].log         
|    └─[id].xyz         # xyz file before optimization.
<-->

├─_job.csv
├─_opt_xtb_data.csv
├─_ref_relation.csv
├─_relation.csv
├─.ini_data.npy
└─tsjob.info
```
The main function of the document is as follows. For detailed information, please refer to the annotations in the document:
```
# Initial Analysis
.ini_data.npy
_relation.csv
_job.csv
sp

# --refine: This function is mainly to unify the forward and reverse reactions in the file _relation.csv and output the refinement results to the file _ref_relation.csv. Reactant, Products(In SMILES).
_ref_relation.csv

# --xtbopt: This function is mainly used for xtb pre optimization of reactants and products. All Optimized Structures are in xtbopt file.
xtbopt
_opt_xtb_data.csv

# --ts: This function is mainly used to generate initial files for transition state structure optimization (See Folder ts) and to count the types of reactions that occur in all trajectories (See tsjob.info).
tsjob.info
ts

# --load: This function is mainly used to detect whether the result file .ini_data.npy obtained from the initial analysis exists. If it exists, the file can be read directly without the need for analysis again.
.ini_data.npy
```

Notes:
1. If you want to obtain Chemical Transformation Network, the function `--refine` must be selected. Please refer to section 6.1 for specific drawing operations.
2. The existence of folder `sp` is to provide information on charges and spin multiplicities corresponding to molecular fragments of reaction events in the gjf file. For more information, please refer to the [original paper](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752).

If you want to obtain further mechanism information, please continue with the subsequent operations.

## 5. TS Location & MEP Search

The purpose of this step is to perform quantum chemical calculations on the structure of the reaction event obtained in the previous step to obtain the chemical information of transition state (TS) structure and minimum energy path (MEP). 

This step mainly includes two parts: preparation for quantum chemistry calculations and execution of quantum chemistry calculations. The script `qmPreprocess.py` and `reactPathRefine.py` are needed. Next, we will introduce the preparation for quantum chemistry calculations.

### 5.1 Preparation for Quantum Chemistry Calculation

In this step, the script `qmPreprocess.py` is needed. The main function of this script is to classify all the gjf file in `ts` generated in the previous step based on the reaction types in `tsjob.info`, preparing for the subsequent operation of the script `reactPathRefine.py`.

Please run the below script for the preliminary TS optimization and collect all information. 

```
qmPreprocess.py -h
```
```
---------------------------------------------------------------------------
usage: qmPreprocess.py [-h] [--step STEP] [--freq FREQ] [--unite] [--sort] [--ts] [--neb]
                       [--path] [--xtbopt]

options:
  -h, --help   show this help message and exit
  --step STEP   
  --freq FREQ   
  --unite       
  --sort        
  --ts
  --neb         # Not yet developed
  --path        
  --xtbopt      # The reactant and product compounds are also optimized with xTB. 
```

where subsequent TS searches require:
```
  --ts  # Make directory named ts (job type is ts).
  --step STEP # Number of sub-folders(default: NONE, include all the sub-folder)
  --unite   # Categorize the same reactions.(eg: A -> B and B -> A)
  --sort    # Sort according to reaction frequency.
```
The relationship between the execution location of script `qmPreprocess.py` and other files is as follows:

```
sampling      # Directory
├─qmPreprocess.py
└─[number]     # Sub-directory of sampling, eg: 1, 2, ...
    ├─ts
    └─tsjob.info
```

#### 5.1.1 Basic Usage of `qmPreprocess.py`:
The execution command of script `qmPreprocess.py` is very simple. It is recommended to use the following command to execute it:

```
qmPreprocess.py --ts --sort --unite
```

This command means that the script will create a series of folders in descending order of frequency based on the reaction type and frequency information read from `tsjob.info`, and place the corresponding ts.gjf files.

#### 5.1.2 Output:

During the process of running script `qmPreprocess.py`, some information will be output. For example:

```
*** Operating Grouping Job-1 ***
  ts/1/1/293_2
  ts/1/1/293_2
  ts/1/1/293_3
  ts/1/1/294_2
  ts/1/1/295_1
  ts/1/1/295_3
  ts/1/1/295_7
  ts/1/1/296_4
  ts/1/1/296_10
  ts/1/1/297_2
  ts/1/1/297_3
  ts/1/1/298_1
  ts/1/1/298_6
  ts/1/1/299_11
  ts/1/1/300_1
  ts/1/1/300_3
  ts/1/1/300_13
  ts/1/1/298_15
  ts/1/1/298_15
  ts/1/2/293_14
  ts/1/2/293_14
  ts/1/2/293_15
  ts/1/2/294_6
  ts/1/2/295_10
  ts/1/2/295_18
  ts/1/2/296_9
  ts/1/2/296_15
  ts/1/2/297_12
  ts/1/2/298_8
  ts/1/2/298_9
  ts/1/2/298_10
  ts/1/2/298_17
  ts/1/3/296_16
  ts/1/3/296_16
  ts/1/3/296_19
  ts/1/3/296_17
  ts/1/3/296_17
  ts/1/3/296_23
  ts/1/4/296_18
  ts/1/4/296_18
  ts/1/4/296_20
  ts/1/4/296_20
  ts/1/5/296_21
  ts/1/5/296_21
  ts/1/5/296_22
  ts/1/5/296_22
  ts/1/6/299_9
  ts/1/6/299_9
  ts/1/6/299_10
  ts/1/6/299_10
  ts/1/7/299_12
  ts/1/7/299_12
  ts/1/7/299_14
  ts/1/7/299_14
  ......
```

The specific newly generated files and directories are as follows

```
sampling      # Directory
├─qmPreprocess.py
├─[number]
└─ts            # ts directory, for ts optimization by gaussian.
  ├─info       # Record the frequency of reaction occurrence.
  └─[number]   # Sub-directory of ts
    └─info     # The reaction type corresponding to this subfolders
```

The directory `ts` contains many subfolders, which contain initial guess result files for ts optimization (such as. gjf) (If you use the step command, there will be multiple subfolders, and if you don't use it, only subfolder '1' will be generated). The `info` file in the subfolder records the information about the reaction (SMILES code and number of occurrences):

```
[H]O+[H]OC1C([H])C([H])C(N([H])C(O)C([H])([H])[H])C([H])C1[H]->[H]C1C([H])C(N([H])C(O)C([H])([H])[H])C([H])C([H])C1O+[H]O[H]   6
[H]C1C([H])C(N([H])C(O)C([H])([H])[H])C([H])C([H])C1O+[H]O->[H]C1C(N([H])C(O)C([H])([H])[H])C([H])C([H])C(O)C1+[H]O[H]   4
[H]C1CC(O)C([H])C([H])C1N([H])C(O)C([H])([H])[H]+[H]O->[H]C1CC(O)C([H])C([H])C1NC(O)C([H])([H])[H]+[H]O[H]   4
[H]C1CC(O)C([H])C([H])C1NC(O)C([H])([H])[H]+[H]O->[H]C1CC(O)C([H])CC1NC(O)C([H])([H])[H]+[H]O[H]   4
```

After script `qmPreprocess.py` processing is completed, the next step execution of quantum chemistry calculations can be carried out.


### 5.2 Execution of Quantum Chemistry Calculations

In this step, the script `reactPathRefine.py` is needed.
The main function of this script is to generate PBS scripts and submit PBS jobs in batches. According to the different tasks of performing quantum chemistry calculations, it can be divided into two parts: TS structure optimization and IRC calculation

This script is mainly based on library Scikit-learn to further analyzes the results obtained in the previous step. Currently, quantum chemistry calculation software Gaussian 16 and Orca are supported.

Before starting quantum chemistry calculations, we need to modify some basic settings. If you want to use the default settings, you can skip this step and run the script directly.

#### 5.2.1 Basic Settings of `reactPathRefine.py`:

The basic settings of this script mainly include the setting of quantum chemistry calculation tasks (such as IRC) and the setting of PBS script. For example, if you want to modify the content of the PBS script, please go to \QmJob\src\_setting_default.py:

```python
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
```
You can make modifications according to your own needs, such as the execution queue, required number of cores, additional commands that need to be executed, and so on.For more detailed settings, please refer to \QmJob\src\_setting_default.py.

After modifying the settings in _setting_default.py according to your needs, you can run the script `reactPathRefine.py` according to your task requirements.

Next, we will use Gaussian 16 as an example to introduce the usage of this script. Firstly, there is the TS optimization function.

#### 5.2.2 TS optimization

Before using this script, please check if the function corresponding to the TS optimization is turned on. If you have already completed it, you can skip this step (5.2.2.1) and move on to the next step. (5.2.2.2)

##### 5.2.2.1 Settings in `reactPathRefine.py`:

The functions required for TS optimization are included in `class G16_TS`. When you use it, please modify the end of 
the script `reactPathRefine. py` (`if __name__=='main__:`) and call `G16_TS`:

```python
if __name__ == '__main__':
    G16_TS(10,5).main()
```

Comment out other classes which are not needed. `G16_TS(m, n)`:
* m: The number of reaction events of a reaction for transition state search by random selection.
* n: The number of GJF files managed by each PBS script.

You can modify the `main` function of `class G16_TS` in `reactPathRefine.py` if needed. 

```python
def main(self):
    ## 1. Initialize
    self.initialize()
    
    ## 2. Submit Jobs 
    submitJobs(self.type, self.state_csv, self.work_space)
    
    ## 3. Move the Calculation Result
    job_list = csv2List(self.job_csv)
    log_list = getIdxLogList(self.work_space)
    moveLogList(job_list, log_list)
    
    ## 4. TS Structure Cluster Analysis
    clusterEnerFreq_G16()
```

For example, if the task submission fails due to a network disconnection, you can comment out the `self.initialize()` in the `main` function and then re-execute `reactPathRefine.py`.


##### 5.2.2.2 Basic Usages:

After modifying the settings of the corresponding function, please use the following command to execute script `reactPathRefine.py` at /ts/1/ (ts folder is generated by `qmPreprocess.py`) :

```
reactPathRefine.py
```

The relationship between the execution location of script `reactPathRefine.py` and other files is as follows:

```
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  └─[1]        # Sub-directory of ts
    └─[numbers]         # Sub-directory of directory 1
      ├─reactPathRefine.py
      ├─[id]        
      | └─[ids].gjf
      └─info
```

Here, the main function of this script is to create a series of folders and PBS scripts based on the parameters set earlier, and submit TS jobs


##### 5.2.2.3 Output:
Files `ts_job.csv`, `ts_state.csv`, and folders `ts_workspace`,`workspace` will be generated:
- `ts_state.csv`: The result of each PBS script after submission.
- `ts_workspace`: All calculated GJF file and PBS script.
- `job_list`: qsub job id.
- `workspace`: Sub-directory include log file which successfully complete ts search.
 

The specific newly generated files and directories are as follows:

```
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  ├─reactPathRefine.py
  └─[1]        # Sub-directory of ts
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
```

Notes: The function corresponding to Files ` ts_job. csv ` has not been implemented, so the default values displayed are False, which is independent of the actual status of the ts task.

After all the steps are completed, each subfolder will have a `workspace`, and then you can proceed to the `Irc` calculation.


#### 5.2.3 IRC

Before using this script, please check if the function corresponding to the IRC is turned on. If you have already completed it, you can skip this step (5.2.3.1) and move on to the next step. (5.2.3.2)

##### 5.2.3.1 Settings in `reactPathRefine.py`:

The functions required for IRC are included in `class G16_Irc`. When you use it, please modify the end of the script `reactPathRefine. py` (`if __name__=='main__:`) and call `G16_Irc`:

```python
if __name__ == '__main__':
    G16_Irc(5).main()
```

Comment out other classes which are not needed. `G16_Irc(m)`:
m: The number of GJF files managed by each PBS script.

You can modify the `main` function of `class G16_Irc` in `reactPathRefine.py` if needed. 

```
def main(self):
    ## 1. Initialize
    # self.initialize_ref() # Use this item only after you have used G16_TS_Ref
    self.initialize()
    
    ## 2. Submit Jobs 
    submitJobs(self.type, self.state_csv, self.work_space)

    ## 3. Move the Calculation Result
    job_list = csv2List(self.job_csv)
    log_list = getIdxLogList(self.work_space)
    moveLogList(job_list, log_list)
```

For example, if the task submission fails due to a network disconnection, you can comment out the `self.initialize()` in the `main` function and then re-execute `reactPathRefine.py`.


##### 5.2.3.2 Basic Usages:

After modifying the settings of the corresponding function, please use the following command to execute script `reactPathRefine.py`:

```
reactPathRefine.py
```
The relationship between the execution location of script `reactPathRefine.py` and other files is as follows:

```
sampling      # Directory
└─ts                # ts directory, for ts optimization by gaussian.
  ├─reactPathRefine.py
  └─[1]        # Sub-directory of ts  
    └─ts_workspace
```

Here, the main function of this script is to create a series of folders and PBS scripts based on the parameters set earlier, and submit IRC jobs

##### 5.2.3.3 Output: 

Files `irc_job.csv`, `irc_state.csv`, and folders `irc_workspace`, `irc`, `irc_temp` will be generated:
- `irc_state.csv`: The result of each PBS script after submission.
- `irc_workspace`: All calculated GJFs, error and pbs file.
- `irc`: Include successful finish log file and gjf file.

The specific newly generated files and directories are as follows:
```
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
```
Notes: The function corresponding to Files `irc_job. csv` has not been implemented, so the default values displayed are False, which is independent of the actual status of the ts task.

If all the steps are successful, you can view the Irc log in each subfolder and draw the reaction mechanism network.

## 6. Network Construction

This section corresponds to the last part of the flowchart, whose main purpose is to construct chemical transformation network and reaction mechanism network.

In this step, the script `grid_draw.py` is required. This script can automatically constructs chemical transformation network and chemical species index based on the file `_ref_relation.csv` generated in section 4.1 by calling the Networkx library and Cmocean library. 

The Reaction Mechanism Network is recommended to use `chemdraw` to draw based on the log file of IRC generated in section 5.2.3.3.

Image Type | Filename | Example
------ | ------ | ------
Chemical Transformation Network | reac_Event_Network.png |Figure 4.a in the [original paper](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752)
Chemical Species Index | mols-RDkit.svg | Figure 4.b in the [original paper](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752)
Reaction Mechanism Network |  | Figure 3 in the [original paper](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00752)

Before using the script `grid_draw.py`, we suggest modifying some basic settings to achieve the best image display effect. If you do not want to modify the default parameters, please go to the 6.2 Basic Usage section

### 6.1 Basic Settings of `grid_draw.py`:

This script mainly includes three aspects of settings: Network Diagram Parameter Settings, Molecular Structure List Parameter Settings and Filter Condition Parameter Settings. Please refer to the following text for specific setting locations. Next, we will introduce these settings one by one.

#### 6.1.1 Network Diagram Parameter Settings:

This section is mainly used to set the relevant parameters of chemical transformation network. In other words, it is to set the display effect of the `reac_Event_Network.png` file.

You can modify the parameters at the `default_net_param` of the `\ReacNetDraw\_network.py` file.

```python
default_net_param = {
    # Node
    'node_size': 600,
    'r_node_size': 2, # Node Scaling

    # Edge
    'edge_width': 2,
    'r_edge_width': 40, # Edge Scaling

    # Transparency
    'r_alpha': 12, # Transparency Scaling
    'node_alpha': 0.5,
    'edge_alpha': 0.25, 
    'node_alpha_criter': 0.6,
    'edge_alpha_criter': 0.7,

    # Arrow
    'bound': [3, 6, 10, 20, 35, 50], 

    # Label
    'cbar_label_size': 20,
    'cbar_tick_size': 14,
    'cbar_label_font': 'Detetected Times',

    'net_label_size': 16,
    'spec_label': {1: '1'},
    'spec_label_color': 'Gray',
    
    # Figure Settings
    'fig_size': [24, 15],
    'fig_dpi': 400,
    'fig_name': 'reac_Event_Network.png'
    }
```

#### 6.1.2 Molecular Structure List Parameter Settings:

This section is mainly used to set the relevant parameters of the chemical species index. In other words, it is to set the display effect of Skeletal formula of chemical spexies in the `mols-RDkit.svg` file.

You can modify the parameters of the molecular structure list diagram at `judgeImgSize` in the `\ReacNetDraw\grid_draw.py` file.

The library `rdkit` is used for plotting. By adjusting the `subImgSize`, the size of the molecule (max_atom_nums) can be modified for optimal display. You can also modify `judgeImgSize` according to the `Max atom number` if the figure is not well drawn.

```python
def judgeImgSize(max_atom_nums):
    if max_atom_nums <= 15:
        subImgSize = (400, 400)
        legendFontSize = 300
    elif max_atom_nums <= 25:
        subImgSize = (600, 600)
        legendFontSize = 300
    elif max_atom_nums <= 35:
        subImgSize = (1600, 1000)
        legendFontSize = 300
    elif max_atom_nums <= 45:
        subImgSize = (1700, 1000)
        legendFontSize = 300
    elif max_atom_nums <= 60:
        subImgSize = (2000, 1300)
        legendFontSize = 300
    print(f' * Sub-Grid Image Size:  {subImgSize}')
    return subImgSize, legendFontSize
```

#### 6.1.3 Filter Condition Parameter Settings:

This section is mainly used to set the Filter Condition of  unsuitable reaction species. You can modify the relevant parameters of `anylseMD` located in the script `md_analy.py`:

```
smi_dict, smi_list, idx_rela_list, filter_smi_dict \
    = anylseMD(freq_criter=2, NC_range=None, Natom_range=None)
```

The specific functions of the parameters are as follows:
- `freq_criter`: `int` type. When the frequency of the reaction is greater than `freq_criter`, the reaction will be shown in the Chemical Transformation Network and Chemical Species Index; otherwise, the reaction will be eliminated.

- `NC_range`: The display range of the number of C atoms, `list` type, e.g. `[1:5]`. `[1:5]` means that when the number of C in the reactant or product is greater than 1 and less than or equal to 5, the reactant and product are retained. Otherwise, the reactant and product will be eliminated.

- `Natom_range`: The range of atomic numbers is displayed, and `list` is supported, e.g. `[1:20]`. `[1:20]` means that when the number of all atoms in the reactants and products is greater than 1 and less than or equal to 20, the reactants and products retain.Otherwise, the reactant and product will be eliminated.

When all relevant settings have been modified, the next step is to use the script `grid_draw.py` to draw the chemical conversion network and chemical species index.

### 6.2 Basic Usage

The usage of this script is very simple, just run the following command in the `sampling` folder:

```
grid_draw.py
```

The relationship between the execution location of script `grid_draw.py` and other files is as follows:

```
sampling      # Directory
├─grid_draw.py
└─[numbers]
  └─_ref_relation.csv
```

This script constructs Chemical Transformation Network and Chemical Species Index by reading the `_ref_relation.csv` file in each subfolders (such as subfolder 293).

### 6.3 Output:
If running normally, the screen will display the following information:
```
 * Reactions Analysis Report
┌────────────────────────────────────────┐
│ Reactive trajectory number:          50│
│ Reaction number:                    412│
│ Reaction class number:              339│
│ Main Reaction number:                73│
│ Main Reaction class number:          11│
└────────────────────────────────────────┘


Drawing _Grid_Image......
 * Max atom number:  31
 * Sub-Grid Image Size:  (1600, 1000)
```

Files `mols-RDkit.svg`, `reac_Event_Network.png`, and `reac_main_eve.log` will be generated:
- `mols-RDkit.svg`: Chemical Species Index
- `reac_Event_Network.png`: Chemical Transformation Network
- `reac_main_eve.log`: SMLIES Code of Reaction

The specific newly generated files and directories are as follows:
```
sampling      # Directory
├─[numbers]
├─mols-RDkit.svg
├─reac_Event_Network.png
├─reac_main_eve.log
└─grid_draw.py
```
>>>>>>> master
