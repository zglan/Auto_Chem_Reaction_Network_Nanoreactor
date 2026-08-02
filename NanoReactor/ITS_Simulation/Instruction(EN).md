# Summary

This program is developed in the group of Prof.Zhenggan Lan at South China Normal University, Guangzhou, China. At present, this program is only a preliminary version and needs further improvement in the future.(As of July 2026)

# 1. Initial Settings:

Assuming you have installed Anaconda, the Python running environment is created by the following command:

```
conda create -n my_env ase=3.29.0 -c conda-forge  
```

The program relies on Python library ASE. If you already have a conda environment, the installation command is:

```
conda install ase=3.29.0 -c conda-forge  
```

# 2. Exploration of  Chemical Space based on ITS Simulation :

Our main goal is the automatic exploration of reaction network and mechanism. 

Generally speaking, the automatic exploration of reaction networks and mechanisms mainly includes three parts: (1) Exploration of chemical reaction spaces; (2) Identification of chemical reactions; and (3) Construction of the reaction network. This program mainly focuses on the exploration of chemical reaction space.

## 2.1. Introduction of the program

The program mainly includes the following files and folder:

```
\src           # Source Code folder
setting.py     
run.py          
qsub.py        # job submission script
job.json       # job template file
```

job.json file is a job template file whose main function is to set the parameters of molecular dynamics simulation. For example:

```
{
    "md": {
        "integrator": 1, 
        "thermostat": {"type": 1, "gamma": 7}, # 1：Langevin Dynamics; gamma, ps-1
        "Temperature": 300, # K
        "dt": 0.5, # fs
        "time": 20, # ps
        "calculator_type": "GFN2-xTB"
    },

    "nano": {
        "k_wall": 0.019, # a.u.
        "beta_wall": 10.0, # Bohr-1
        "R_wall": 3.59 # Bohr
    },
	
	"its": {
		"bin_temp": 136, # number of temperature bin
		"temp_low": 300, # K
		"temp_high": 3000 # K
	}
    
}
```

The constrained potential energy for nanoreactor is :

$$
V_{sphere} = \sum_A k_{wall} log\{1 + \exp[\beta_{wall} (|R_A - R_0| - R_{wall})]\}
$$

where $k_{wall}$ is the coefficient to scale the wall potential energy intensity, $\beta_{wall}$ controls the steepness of the wall potential energy change , $R_A$ is the Cartesian coordinate of atom $A$, $R_0$ is the center of the reaction system, and $R_{wall}$ is the radius of the confining sphere.
Following the suggestions of previous work, default values of 
$k_{wall} = 0.019 \, a.u.$ and $\beta_{wall} = 10 \, Bohr^{-1}$  are used.

The Integrated Tempering Sampling(ITS) method is also included. In the ITS method, a non-Boltzmann factor is defined:

$$
f_{ITS}(U) = \sum_k^N n_k \exp(-\beta_k U) 
$$

where $n_k$ is weighting factor, $\beta_k=\dfrac{1}{k_bT_k}$.

In ITS-enhanced molecular dynamics simulations or Monte Carlo simulations, the following relationships exist:

$$
f_{ITS}(U) = \exp(-\beta_0 U_{ITS})
$$

where $U_{ITS}$ is the effective potential energy of the effective system:

$$
U_{ITS} = -\dfrac{1}{\beta_0} \ln \sum_k^N n_k \exp(-\beta_k U) 
$$

Simultaneously, the effective force $F_{ITS}$ of the effective system can be exported:

$$
F_{ITS} = -\dfrac{\partial U_{ITS}}{\partial r} = \dfrac{\sum_k^N n_k \cdot \beta_k \cdot \exp(-\beta_k U)}{\beta_0 \sum_k^N n_k \cdot \exp(-\beta_k U)} \cdot F = S \cdot F
$$

where $S$ is so-called enhanced factor.

Before applying ITS-enhanced MD, however, there is no analytical expression for the weight factor $n_k$, and its value is often unknown beforehand. Therefore, an iterative method will be used to obtain the value of the weight factor $n_k$:

* `1. Initialize`: 
  
 Divide the temperature range [$T_{min}$, $T_{max}$] into N intervals, with each interval corresponding to $\beta_k$(k=1,..., N). A set of initial $n_k$ values is pre-set and the ratio between adjacent values is calculated:

 $$
 m_k = \dfrac{n_{k}}{n_{k+1}}
 $$

* `2. ITS-enhanced MD Simulation`: 

Calculate the effective $F_{ITS}$. At the same time, constrained potential energy is also applied to prevent molecules from escaping from the nanoreactor.

* `3. Evoluate $p_k$`:

 During Simulation, the configuration integral $P_k$ is estimated based on the sampling situation. Subsequently, $p_k$, which is normalized $P_k$ coresponed to $\beta_k$, is calculated by:

 $$
    p_k = \frac{P_k}{\sum_l P_l}
$$

where:

$$
P_k = n_k \int \exp(-\beta_k U)dr
$$

* `4. Update $m_k$`:

 In the i-th iteration step, a weighing function $w_k(i)$ is firstly calculated by:

 $$
w_k(i) = \sqrt{ p_k(i-1)p_{k+1}(i-1)},
 $$

And the  cumulative historical sum of $w_k(i)$, 
$W_k(i)$, is also  calculated by

$$
W_k(i)=\sum_{l=0}^{i-1} w_k(l).
$$

Next a bias term $w'_k(i) = c_w w_k(i)$ is introduced with the bias factor $c_w$. The default value of bias factor $c_w$ is set to 1 in practices.
The adjusted cumulative term $W'_k(i)$ is then obtained:

$$
W'_k(i) = W_k(i-1) + w'_k(i).
$$

Finally，$m_k$ will be updated according to:

$$
m_k(i) = \frac{ m_k(i-1) }{  W'_k(i-1) + w'_k(i) } \, 
\left [ W'_k(i-1) + w'_k(i) \dfrac{p_{k+1}(i-1)}{p_k(i-1)} \right ] 
$$

Once a new set of $m_k$ is obtained, a new set of $n_k$ can be obtained.

* `5. Iteration`:
  
  Repeat steps (2) to (4). When $p_{k+1}\approx p_k$ and $m_k(i)\approx m_k(i-1)$, the $n_k$ of syteam reaches convergence for all temperatures $T_k$.

## 2.2. Application of the program

Overall, the process of automated exploration of reaction networks and mechanisms has been introduced in [README.md](https://github.com/zglan/Auto_Chem_Reaction_Network_Nanoreactor/blob/main/README.md), including initial sampling, chemical space exploration, chemical event analysis, refinement and reaction network construction. This work mainly primarily focuses on exploring the chemical reaction space. In this section, we assume that the initial sampling has been completed.

Before the formal simulation, it is necessary to set the relevant parameters:

<table style="margin-left: auto; margin-right: auto;">
  <tr>
    <th>Job Type</th>
    <th>Key</th>
    <th>default</th>
    <th>description</th>
  </tr>
  <!-- Parameters of md -->
  <tr>
    <td rowspan="6">md</td>
    <td>integrator</td>
    <td>1</td>
    <td>Integrator of MD(default: velocity-verlet)</td>
  </tr>
  <tr>
    <td>thermostat</td>
    <td>"type": 1, "gamma": 7</td>
    <td>Thermostat of MD(default: 1. type: Langevin dynamics; 2. gamma: 7 ps-1)</td>
  </tr>
  <tr>
    <td>Temperature</td>
    <td>None</td>
    <td>Target temperature of system; K</td>
  </tr>
  <tr>
    <td>dt</td>
    <td>1.0</td>
    <td>Time step of md(default: 1.0 fs)</td>
  </tr>
  <tr>
    <td>time</td>
    <td>5</td>
    <td>Total time of md(default: 5 ps)</td>
  </tr>
  <tr>
    <td>calculator_type</td>
    <td>"GFN2-xTB"</td>
    <td>Type of calculator(default: "GFN2-xTB")</td>
  </tr>
  <!-- Parameters of nano -->
  <tr>
    <td rowspan="3">nano</td>
    <td>k_wall</td>
    <td>0.019</td>
    <td>The coefficient to scale the wall potential energy intensity(default: 0.019 a.u.)</td>
  </tr>
  
  <tr>
    <td>beta_wall</td>
    <td>10.0</td>
    <td>The steepness of the wall potential energy change(default: 10.0 Bohr<sup>-1</sup>)</td>
  </tr>
  <tr>
    <td>R_wall</td>
    <td>None</td>
    <td>The radius of nanoreactor; Bohr</td>
  </tr>
  <!-- Parameters of its -->
  <tr>
    <td rowspan="3">its</td>
    <td>bin_temp</td>
    <td>10</td>
    <td>Number of temperature bins (default: 10)</td>
  </tr>
  <tr>
    <td>temp_low</td>
    <td>300</td>
    <td>Lower limit of temperature range(default: 300.0 K)</td>
  </tr>
  <tr>
    <td>temp_high</td>
    <td>650</td>
    <td>Upper limit of temperature range(default: 650.0 K)</td>
  </tr>
</table>

It should be noted that in ITS simulation job, in addition to the temperature range and the number of segmented intervals, the initial nk value also needs to be defined. The specific location is in the script `its_main.py`:

```
...

class ITS:

...

def __init__(self,
                 T_target: float,
                 temp_low: float,
                 temp_high: float,
                 bin_temp: int,
                 update_step: int = 500):
... 
self.nk = ...
...
```

The initial value of $n_k$ can be modified according to actual needs. For example, $n_k = 2^{-k+1} (k \in [1,...,N])$.

After all the settings are completed, the qsub.sh script can be used for task submission. The following is an example of a job submission script qsub.sh that can be modified according to the actual situation: 

```
#PBS -N job_name
#PBS -j oe
#PBS -q amd
#PBS -l nodes=1:ppn=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_STACKSIZE=1000m
ulimit -s unlimited

user="XXX"

source /Pathway/To/conda.sh
conda activate my_env
export PYTHONPATH=/Pathway/To/ITS_Simulation/:$PYTHONPATH

CURR=$PBS_O_WORKDIR
WORK_DIR=$CURR
TMP_DIR="XXXXXX"   

cd $PBS_O_WORKDIR
mkdir -p  $TMP_DIR
cp -r $WORK_DIR/* $TMP_DIR
cd $TMP_DIR

pwd > log
python3 run.py

cp -rf $TMP_DIR/*   $WORK_DIR
rm -rf $TMP_DIR
```

After the modification is completed, the following command can be used to execute it:

```
qsub qsub.sh
```

After successful execution, an md.xyz file will be generated, which can be used for subsequent analysis.(Section 4: Reaction Event Identification in [README.md](https://github.com/zglan/Auto_Chem_Reaction_Network_Nanoreactor/blob/main/README.md)) The detail information of subsequent steps can refer to [README.md](https://github.com/zglan/Auto_Chem_Reaction_Network_Nanoreactor/blob/main/README.md).

# 3. Citation of related work

If you have applied the program in this section, it is recommended to cite the following references:

* `ITS Principle`：
***J. Chem. Phys. 128, 064105 (2008)***

* `Iteration Formula of $n_k$`：
***J. Chem. Phys. 128, 134111 (2008)***

* `Langevin Dynamics`：
***Phys. Rev. E 75, 056707 (2007)***

* `ASE`：
***J. Phys.: Condens. Matter 29 273002 (2017)***

* `Constrained Potential Energy`：
***J. Chem. Theory Comput. 2019, 15, 5, 2847–2862*** 