# TDCI-CAP

This code models ionization by performing state-based Time Dependent Configuration Interaction (TDCI) in the presence of a complex absorbing potential. Supports CI Singles (TD-CIS) for modeling first ionization, and TD-CISD-IP for second ionization.

## Compilation

TDCI-CAP can be compiled with nvfortran on unix systems supported by NVIDIA's HPC SDK. We have tested only CentOS, Ubuntu, and Archlinux x64 builds.

The NVIDIA HPC SDK can be downloaded at https://developer.nvidia.com/hpc-sdk-downloads
NVIDIA's install instructions https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html

Once the environment is prepared, you can clone the repo and compile the tdci code.

```
git clone https://github.com/hbschlegel/TDCI-CAP
cd TDCI-CAP
make
```

Once the compilation is complete, you will find the tdci executable in the `bin/` subdirectory.

## Tutorial




### Gaussian Input
Before running TDCI, we must run a preliminary Gaussian Configuration Interaction Singles (CIS) calculation to obtain the one and two electron integrals, dipole, and absorbing potential matrix elements.
We use Gaussian development version gdvj28p.
Below is an example Gaussian input file. Several nonstandard inputs are necessary, and are described in more detail below.


```bash
%chk=h2o.chk
%mem=10GB
%nproc=10
#P ucis(mo,nstates=5)/aug-cc-pVDZ extrabasis gfprint
int(acc2e=12) Pop=Full IOp(3/194=12003501,6/18=1,8/10=91,8/37=2,8/38=-10,9/127=2)
OUTPUT=(faf,i4lab,Files=(619,751,870,1,2,3,4,5,6,7,8,9,10,11,12,13,14))

h2o

0 1
O        0.000000    0.000000    0.110812
H        0.000000   -0.783976   -0.443248
H        0.000000    0.783976   -0.443248

O H 0
 s 1 1.00
 .0128000000D-00 1.0000000000D+00
 s 1 1.00
 .0064000000D-00 1.0000000000D+00
 p 1 1.00
 .0128000000D-00 1.0000000000D+00
****

MatrixElements.faf


```

`IOp(3/194)` can set the absorbing potential parameters proportional to each atom's van der Waals radius in the following format:
`aaabbb01`
Where the outer (inner) boundary is set at aaa/10 (bbb/10) times the atom's van der Waals radius, respectively.
For example, the `3/194=12003501` in the example sets the inner boundary to 3.5 radii, and the outer boundary to 12 radii.

`IOp(3/195)` sets the maximum value of the absorbing potential in atomic units.

`IOp(6/18)=1` and `IOp(9/127)=2` are necessary for writing the MatrixElements.faf

`IOp(8/10)=91` enables `IOp(8/37-38)`, which specify the orbital cutoffs.

`IOp(8/37)=2` selects the first occupied orbital to be included in the active space. Excitations out of orbitals below will not be considered.

`IOp(8/38)=-10` selects the end of the active space. If `IOp(8/38)=N` where N is positive, the last virtual orbital in the active space will be orbital `N-1`. If N is negative, the N highest energy virtual orbitals will be removed from the active space. In practice, we select this to ignore orbitals above about 3.0 Hartree. 

Don't change the line that starts with `OUTPUT=`, it specifies the indices of the various files to be included in the MatrixElements.faf file.

When selecting an active space, one must be careful to ensure degenerate orbitals remain in the same (non)-active set. If one degenerate orbital is left out of the active space while others are included, it will cause unphysical results.

If you're having convergence difficulties, try combinations of these parameters: `int(acc2e=12) SCF(QC,NoVarAcc,NoIncFock,MaxCYC=64,VShift=500)`

For more information on the IOps, consult the [Gaussian overlay documentation](https://gaussian.com/overlay9/). The SCF parameters are documented [here](https://gaussian.com/scf/).

NOTE: If you are using shared computing like WSU's Grid, please remember to follow their rules and policies. This means you should run Gaussian and TDCI from an sbatch script that is scheduled through SLURM. See your computing resource's documentation for more details.

Once your Gaussian calculation successfully finishes, there should be a MatrixElements.faf in the Gaussian job directory.

### TDCI input

TDCI expects both 'MatrixElements.faf' and 'input' to be in the current working directory. Sample tdci input files can be found in the test/tests folder. To perform a tdci simulation, simply execute the tdci binary (or sbatch script) from inside the job directory. We recommend separating the gaussian and tdci job directories, and creating symlink to the MatrixElements.faf in your tdci directory.
Your setup should look something like this:

```bash
~$ cd tdci_h2o
~$ ln -s ../gaussian_h2o/MatrixElements.faf .  # Create symlink.
~/tdci_h2o$ ls -alt
drwxr-xr-x 18 user user     4096 Sep 23 16:47 ..
drwxr-xr-x 18 user user     4096 Sep 23 16:47 .
-rw-r--r--  1 user user     2052 Sep 22 18:42 input
lrwxrwxrwx  1 user user       39 Sep 22 18:42 MatrixElements.faf -> ../gaussian_h2o/MatrixElements.faf

~/tdci_h2o$ ~/TDCI-CAP/bin/tdci

```

TDCI itself exports some data tables in the job directory, but TDCI is also packaged with a python analysis script that generates many tables and plots. The tools/rate\_analyzer.py script must be executed from the job directory, and passed the path of the gaussian output log as an argument.
```
~/tdci_h2o$ python3 ~/TDCI-CAP/tools/rate_analyzer.py ../gaussian_h2o/h2o.log
```

Next, we provide a sample tdci `input` file. A table with descriptions of each parameter is provided below.

```bash
 &DYNAMICS
 init_states(0) = 2
 init_states(1) = 1
 init_states(1) = 2
 init_coeffs(1) = (1.00000,0.00000)
 restart        = .false.
 Qsave          = .false.
 /
 &FIELD_units
 omega_units = 'au'
 /
 &FIELD
 dirform  = 'polar'
 ellipt   =  1.000
 envelope = 'stat'
 ncyc     =     7
 omega    =  0.057
 phase    = 90.000
 euler    =  0.000
 /
&FIELD_strengths 
 nemax =       1
 read_emax(1) = 0.0500
 /
&FIELD_directions
 ndir = 2 
 read_theta( 1) = 0.d0   ;  read_phi( 1) = 0.d0 
 read_theta( 2) = 30.d0  ;  read_phi( 2) = 0.d0 
/
 &SYSTEM_units
 dt_units     = 'au' 
 eigmax_units = 'au' 
 /
 &SYSTEM
 dt          =  0.050
 eigmax      = 10.000
 ionization  = -1.000
 jobtype     =  'cis'
 nstep       =   1000
 outstep     =    100
 nactive     =      0
 nvirtual    =      0
 socfac      = 1.00000000
 socfacZ     = 1.00100000
 IP_alpha_beta = .false.
 QeigenDC = .true.
 /
 &InOutFILES
 Qread_tdcidata  = .false.
 tdcidatfile     = 'MatrixElements.faf'
 outputfile      = 'OUTPUT'
 restartbinfile  = 'OUTPUT_RESTART.bin'
 Qwrite_ion_coeff= .false.
 Qread_ion_coeff = .false.
 Qmo_dens        = .true.
 write_binaries = .false.
 /
&Davidson
 flag_davidson = .False. 
 /
&ReadU_NO
 flag_ReadU_NO = .True.
 /


```

Table of input parameters:


| **Namelist**       | **Parameter**       | **Description**                                                                                                     | **Allowed values**                          | **Notes**                                                                                                                                                        |
|--------------------|---------------------|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DYNAMICS           | Init_states(i)      | Init_state(0) is the number of initial states. Init_state(i) for i>0 is the state index of state i.                 | N/A                                         |                                                                                                                                                                   |
| DYNAMICS           | Init_coeffs(i)      | Initial coefficient of state i                                                                                      | Complex (e.g., (1.00, 0.00))                | Example: Init_coeffs(1) = (1.0, 0.0)                                                                                                                              |
| DYNAMICS           | restart             | Reads RESTART.bin to restart the propagation                                                                       | .true. or .false.                           |                                                                                                                                                                   |
| DYNAMICS           | Qsave               | Creates RESTART.bin                                                                                                | N/A                                         |                                                                                                                                                                   |
| FIELD_units        | Omega_units         | Units for the omega parameter                                                                                       | ‘au’                                        | ‘au’ – Atomic units.                                                                                                                                             |
| FIELD              | Dirform             | Coordinate system used for FIELD_directions namelist.                                                              | ‘polar’, ‘cart’                             | ‘polar’: Polar coordinates<br>‘cart’: Cartesian coordinates                                                                                                       |
| FIELD              | Ellipt              | Ellipticity for circularly polarized fields                                                                        | Float from 0.0–1.0                          | A value of 1.0 corresponds to a circular pulse.                                                                                                                   |
| FIELD              | Envelope            | Envelope that outlines the carrier field.                                                                          | ‘none’, ‘cos2’, ‘trap’, ‘gaus’, ‘stat’, ‘band’, ‘ramp’, ‘sin2’, ‘cirl’, ‘cirr’ | For a full description of the options, see the relevant section in documentation.                                                                                 |
| FIELD              | Ncyc                | Number of carrier field periods that are enclosed in the envelope                                                  | Positive integer                            | Does not apply to fields ‘none’, ‘static’, and ‘ramp’.                                                                                                            |
| FIELD              | Omega               | Frequency of the carrier field                                                                                     | Positive float                              | Units are set by the omega_units parameter.                                                                                                                       |
| FIELD              | Phase               | Carrier wave’s phase factor in degrees                                                                             | Positive float                              |                                                                                                                                                                   |
| FIELD              | Euler               | Sets the Euler angle (gamma) for elliptical pulses.                                                                | Non-negative float                          |                                                                                                                                                                   |
| FIELD_strengths    | nemax               | Propagations to perform with different field strengths                                                             | Positive integer                            | Read_emax(i) must be set for i=1,…,nemax                                                                                                                           |
| FIELD_strengths    | Read_emax(i)        | Maximum field strength of field i.                                                                                 | Non-negative float (i must be ≥ 1)          |                                                                                                                                                                   |
| FIELD_directions   | ndir                | Number of propagations to perform with different field polarization directions.                                    | Positive integer                            | Read_theta(i) and Read_phi(i) must be set for i=1,…,ndir. The total number of propagations will be nemax×ndir.                                                    |
| FIELD_directions   | Read_theta(i)       | Theta coordinate for field polarization direction i                                                                | Float from 0–360                            | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | Read_phi(i)         | Phi coordinate for field polarization direction i                                                                  | Float from 0–360                            | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | Read_x(i)           | x coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | Read_y(i)           | y coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | Read_z(i)           | z coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| SYSTEM_units       | dt_units            | Units for the ‘dt’ parameter.                                                                                      | ‘au’                                        | ‘au’ – Atomic units                                                                                                                                              |
| SYSTEM_units       | Eigmax_units        | Units for the ‘eigmax’ parameter.                                                                                  | ‘au’                                        | ‘au’ – Atomic units                                                                                                                                              |
| SYSTEM             | Dt                  | Size of the propagation timestep                                                                                   | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | Eigmax              | Ignore states above this energy                                                                                    | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | Ionization          | Ionization parameter (details unspecified)                                                                         | Float                                       |                                                                                                                                                                   |
| SYSTEM             | Jobtype             | Selects the type of propagation to be performed                                                                    | ‘cis’, ‘cisip’, ‘soc’, ‘socip’              | ‘cis’: TD-CIS<br>‘cisip’: TD-CISD-ip<br>‘soc’: TD-CIS with spin-orbit coupling<br>‘socip’: TD-CISD-ip with spin-orbit coupling                                   |
| SYSTEM             | Nstep               | Number of timesteps to be propagated                                                                               | Positive integer                            |                                                                                                                                                                   |
| SYSTEM             | Outstep             | Output data will only be written on step numbers divisible by outstep                                              | Positive integer                            | If nstep=50 and outstep=10, data will be written on steps 10, 20, 30, 40, and 50.                                                                                 |
| SYSTEM             | Nactive             | Number of active occupied orbitals (lower orbitals are ignored)                                                    | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | Nvirtual            | Number of active virtual orbitals (higher orbitals are ignored)                                                    | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | socfac              | Spin-orbit coupling factor (details unspecified)                                                                   | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | socfacZ             | Scales z-component of spin-orbit coupling (development test)                                                       | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | IP_alpha_beta       |                                                                                                                     | .true. or .false.                           |                                                                                                                                                                   |
| SYSTEM             | QeigenDC            |                                                                                                                     | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_tdcidata      | Type of Gaussian data file to read                                                                                 | .true. or .false.                           | For legacy TDCI.dat files, use .true.<br>For new MatrixElements.faf files, use .false.                                                                            |
| InOutFILES         | tdcidatfile         | Specifies the filename of the Gaussian data file                                                                   | String                                      |                                                                                                                                                                   |
| InOutFILES         | Outputfile          | Name of output log file                                                                                            | String                                      |                                                                                                                                                                   |
| InOutFILES         | Restartbinfile      | Name of the restart file                                                                                           | String                                      | Setting parameter ‘restart’ to .true. will read from this file.<br>Setting parameter ‘Qsave’ to .true. will write this file.                                     |
| InOutFILES         | Qwrite_ion_coeff    | For sequential double ionization, writes the ion coefficients                                                      | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_ion_coeff     | For sequential double ionization, reads the ion coefficients                                                       | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qmo_dens            | Writes the density matrix to MO_density files                                                                      | .true. or .false.                           |                                                                                                                                                                   |
| Davidson           | Flag_davidson       | Use Davidson diagonalization for CISD (not implemented)                                                            | .true. or .false.                           |                                                                                                                                                                   |





## How to cite








