
## Theoretical Background

The time-dependent configuration interaction (TDCI) code can be used to simulate the electronic structure of a molecule interacting with the electric field of an intense, ultra-short laser pulse.
The time-independent, ground and excited configurations of a molecule in the absence of a field are combined with time-dependent linear coefficients to model the molecular wavefunction in the electric field of the laser pulse (the nuclear positions kept fixed).
The interaction of the molecule with the electric field of the laser pulse is treated in the semiclassical dipole approximation and the pulses can be linear or elliptically polarized.
For strong field ionization, the molecule is surrounded by a complex absorbing potential (CAP) to remove the outgoing electron (see Figure 1).  

![Figure 1](img/Figure1.png)
***Figure 1.** Hydrogen atom wavefunction (light blue), Coulomb potential (dark blue) and complex absorbing potential (dashed dark blue) without a field. The strong field from a laser pulse suppresses the Coulomb potential (red) and the wavefunction (orange) can go over the Coulomb barrier or tunnel through it and be absorbed by the complex potential (dashed red).*

For single ionization, the wavefunction includes the Hartree-Fock ground state and all singly excited configurations (CIS).
Ionization of cations can be modeled by a CISD-IP wavefunction which is built from singly ionized and singly excited, singly ionized configurations of the Hartree-Fock ground state of the neutral molecule.
For molecules with heavy elements, spin-orbit coupling (SOC) can be included with an effective one electron SOC operator.
The wavefunction is propagated with the exponential of the field-dependent Hamiltonian using a Trotter factorization.
The ionization yield for TDCI-CAP simulation is calculated from the decrease in the norm squared of the wavefunction as it interacts with the absorbing potential.
The rate is computed from the expectation value of the absorbing potential using the time-dependent wavefunction.
The integrals needed to construct the Hamiltonian are computed in an electronic structure program such as Gaussian.
For strong field ionization standard molecular basis sets are augmented with several sets of diffuse functions to support the wavefunction as the electron leaves the molecule and interacts with the absorbing potential.
The followings sections provide a brief overview of the theory and computational approaches of the TDCI code.
More details on applications can be found in some of our papers.[^1][^2][^3][^4][^5][^6][^7]

### Hamiltonian 
In the TD-CI approach, the electronic wavefunction is propagated with the time dependent Schrödinger equation.
$$
i\hbar \frac{\partial}{\partial t} \Psi(t) = \hat{H}(t)\Psi(t) = [\hat{H}_{el} + \hat{V}^{SOC} - \hat{\vec{\mu}} \cdot \vec{E}(t) - i \hat{V}^{abs}] \Psi(t) \tag{1}
$$
$\hat{H}_{el}$ is the field-free, time-independent non-relativistic electronic Hamiltonian and the interaction with the electric field of an intense laser pulse is treated in the semiclassical dipole approximation, where $\hat{\vec{\mu}}$ is the dipole operator and $\vec{E}$ is the electric field.
The spin-orbit coupling term is approximated by an effective one electron spin-orbit coupling operator. 
$$
\hat{V}^{SOC} = -\frac{\alpha_0^2}{2} \sum_A Z_A^{eff} \frac{(\vec{r} - \vec{r}_A) \times \nabla}{|\vec{r} - \vec{r}_A|^3} \tag{2}
$$
Suitable values for $Z^{eff}$ have been reported by Koseki, Gordon and co-workers and by Chiodo and Russo.
$Z^{eff}$ can be adjusted to reproduce the experimentally observed spin-orbit splitting.
$\hat{V}^{abs}$  is the absorbing potential.

### Absorbing Potential
The absorbing potential, $\hat{V}^{abs}$, surrounds the molecule starting at a distance of about 10 – 15 bohr and removes the outgoing part of the wavefunction thereby simulating ionization.
The absorbing potential for the molecule is constructed from spherical potentials centered on each atom.
The potential on atom C is written in terms of the distance from the nucleus, $r_c = |\vec{r} - \vec{C}|$.
In the quadratic form, the atomic potential starts at $R_A$ rises quadratically to $(R_A+R_B)/2$, turns over quadratically to $R_B$ and is constant, $V_{max}$, beyond $R_B$. 

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
 init_states(2) = 3
 init_coeffs(1) = (1.00000,0.00000)
 init_coeffs(2) = (1.00000,0.00000)
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

```

Table of input parameters:


| **Namelist**       | **Parameter**       | **Description**                                                                                                     | **Allowed values**                          | **Notes**                                                                                                                                                        |
|--------------------|---------------------|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DYNAMICS           | init_states(i)      | Init_state(0) is the number of initial states. Init_state(i) for i>0 is the state index of state i.                 | N/A                                         |                                                                                                                                                                   |
| DYNAMICS           | init_coeffs(i)      | Initial coefficient of state i                                                                                      | Complex (e.g., (1.00, 0.00))                | Example: Init_coeffs(1) = (1.0, 0.0)                                                                                                                              |
| DYNAMICS           | restart             | Reads RESTART.bin to restart the propagation                                                                       | .true. or .false.                           |                                                                                                                                                                   |
| DYNAMICS           | Qsave               | Creates RESTART.bin                                                                                                | N/A                                         |                                                                                                                                                                   |
| FIELD_units        | omega_units         | Units for the omega parameter                                                                                       | ‘au’                                        | ‘au’ – Atomic units.                                                                                                                                             |
| FIELD              | dirform             | Coordinate system used for FIELD_directions namelist.                                                              | ‘polar’, ‘cart’                             | ‘polar’: Polar coordinates<br>‘cart’: Cartesian coordinates                                                                                                       |
| FIELD              | ellipt              | Ellipticity for circularly polarized fields                                                                        | Float from 0.0–1.0                          | A value of 1.0 corresponds to a circular pulse.                                                                                                                   |
| FIELD              | envelope            | Envelope that outlines the carrier field.                                                                          | ‘none’, ‘cos2’, ‘trap’, ‘gaus’, ‘stat’, ‘band’, ‘ramp’, ‘sin2’, ‘cirl’, ‘cirr’ | For a full description of the options, see the relevant section in documentation.                                              |
| FIELD              | ncyc                | Number of carrier field periods that are enclosed in the envelope                                                  | Positive integer                            | Does not apply to fields ‘none’, ‘static’, and ‘ramp’.                                                                                                            |
| FIELD              | omega               | Frequency of the carrier field                                                                                     | Positive float                              | Units are set by the omega_units parameter.                                                                                                                       |
| FIELD              | phase               | Carrier wave’s phase factor in degrees                                                                             | Positive float                              |                                                                                                                                                                   |
| FIELD              | euler               | Sets the Euler angle (gamma) for elliptical pulses.                                                                | Non-negative float                          |                                                                                                                                                                   |
| FIELD_strengths    | nemax               | Propagations to perform with different field strengths                                                             | Positive integer                            | Read_emax(i) must be set for i=1,…,nemax                                                                                                                           |
| FIELD_strengths    | read_emax(i)        | Maximum field strength of field i.                                                                                 | Non-negative float (i must be ≥ 1)          |                                                                                                                                                                   |
| FIELD_directions   | ndir                | Number of propagations to perform with different field polarization directions.                                    | Positive integer                            | Read_theta(i) and Read_phi(i) must be set for i=1,…,ndir. The total number of propagations will be nemax×ndir.                                                    |
| FIELD_directions   | read_theta(i)       | Theta coordinate for field polarization direction i                                                                | Float from 0–360                            | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | read_phi(i)         | Phi coordinate for field polarization direction i                                                                  | Float from 0–360                            | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | read_x(i)           | x coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | read_y(i)           | y coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | read_z(i)           | z coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| SYSTEM_units       | dt_units            | Units for the ‘dt’ parameter.                                                                                      | ‘au’                                        | ‘au’ – Atomic units                                                                                                                                              |
| SYSTEM_units       | eigmax_units        | Units for the ‘eigmax’ parameter.                                                                                  | ‘au’                                        | ‘au’ – Atomic units                                                                                                                                              |
| SYSTEM             | dt                  | Size of the propagation timestep                                                                                   | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | eigmax              | Ignore states above this energy                                                                                    | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | ionization          | Ionization parameter (details unspecified)                                                                         | Float                                       |                                                                                                                                                                   |
| SYSTEM             | jobtype             | Selects the type of propagation to be performed                                                                    | ‘cis’, ‘cisip’, ‘soc’, ‘socip’              | ‘cis’: TD-CIS<br>‘cisip’: TD-CISD-ip<br>‘soc’: TD-CIS with spin-orbit coupling<br>‘socip’: TD-CISD-ip with spin-orbit coupling                                   |
| SYSTEM             | nstep               | Number of timesteps to be propagated                                                                               | Positive integer                            |                                                                                                                                                                   |
| SYSTEM             | outstep             | Output data will only be written on step numbers divisible by outstep                                              | Positive integer                            | If nstep=50 and outstep=10, data will be written on steps 10, 20, 30, 40, and 50.                                                                                 |
| SYSTEM             | nactive             | Number of active occupied orbitals (lower orbitals are ignored)                                                    | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | nvirtual            | Number of active virtual orbitals (higher orbitals are ignored)                                                    | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | socfac              | Spin-orbit coupling factor (details unspecified)                                                                   | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | socfacZ             | Scales z-component of spin-orbit coupling (development test)                                                       | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | IP_alpha_beta       |                                                                                                                     | .true. or .false.                           |                                                                                                                                                                   |
| SYSTEM             | QeigenDC            |                                                                                                                     | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_tdcidata      | Type of Gaussian data file to read                                                                                 | .true. or .false.                           | For legacy TDCI.dat files, use .true.<br>For new MatrixElements.faf files, use .false.                                                                            |
| InOutFILES         | tdcidatfile         | Specifies the filename of the Gaussian data file                                                                   | String                                      |                                                                                                                                                                   |
| InOutFILES         | outputfile          | Name of output log file                                                                                            | String                                      |                                                                                                                                                                   |
| InOutFILES         | restartbinfile      | Name of the restart file                                                                                           | String                                      | Setting parameter ‘restart’ to .true. will read from this file.<br>Setting parameter ‘Qsave’ to .true. will write this file.                                     |
| InOutFILES         | Qwrite_ion_coeff    | For sequential double ionization, writes the ion coefficients                                                      | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_ion_coeff     | For sequential double ionization, reads the ion coefficients                                                       | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qmo_dens            | Writes the density matrix to MO_density files                                                                      | .true. or .false.                           |                                                                                                                                                                   |
| Davidson           | flag_davidson       | Use Davidson diagonalization for CISD (not implemented)                                                            | .true. or .false.                           |                                                                                                                                                                   |





## How to cite



[^1]: [Krause, P.; Sonk, J. A.; Schlegel, H. B., Strong field ionization rates simulated with time-dependent configuration interaction and an absorbing potential. J. Chem. Phys. 2014, 140, 174113. 10.1063/1.4874156](https://doi.org/10.1063/1.4874156) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/369.pdf))
[^2]: [Krause, P.; Schlegel, H. B., Angle-dependent ionization of small molecules by time-dependent configuration interaction and an absorbing potential. J. Phys. Chem. Lett. 2015, 6, 2140-2146. 10.1021/acs.jpclett.5b00929](https://doi.org/10.1021/acs.jpclett.5b00929) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/381.pdf))
[^3]: [Hoerner, P.; Schlegel, H. B., Angular dependence of strong field ionization of CH3X (X = F, Cl, Br, or I) using time-dependent configuration interaction with an absorbing potential. J. Phys. Chem. A 2017, 121, 5940-5946. 10.1021/acs.jpca.7b06108](https://doi.org/10.1021/acs.jpca.7b06108) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/409.pdf))
[^4]: [Hoerner, P.; Schlegel, H. B., Angular dependence of strong field ionization of haloacetylenes, HCCX (X = F, Cl, Br, I) using time-dependent configuration interaction with an absorbing potential. J. Phys. Chem. C 2018, 122, 13751–13757. 10.1021/acs.jpcc.8b00619](https://doi.org/10.1021/acs.jpcc.8b00619) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/416.pdf))
[^5]: [Lee, M. K.; Hoerner, P.; Li, W.; Schlegel, H. B., Effect of spin-orbit coupling on strong field ionization simulated with time-dependent configuration interaction. J. Chem. Phys. 2020, 153. 10.1063/5.0034807](https://doi.org/10.1063/5.0034807) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/435.pdf))
[^6]: [Lee, M. K.; Li, W.; Schlegel, H. B., Angular dependence of strong field sequential double ionization for neon and acetylene simulated with time-dependent configuration interaction using CIS and CISD-IP. J. Chem. Phys. 2020, 152, 064106. 10.1063/1.5133659](https://doi.org/10.1063/1.5133659) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/429.pdf))
[^7]: [Schlegel, H. B.; Hoerner, P.; Li, W., Ionization of HCCI neutral and cations by strong laser fields simulated with time dependent configuration interaction. Front. Chem. 2022, 10, 866137. 10.3389/fchem.2022.866137](https://doi.org/10.3389/fchem.2022.866137) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/445.pdf))


