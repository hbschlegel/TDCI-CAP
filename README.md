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

To test your installation, you can use the test utility in the test folder
```
python3 test/test.py 
```

## Usage

TDCI uses integrals and other data from a Gaussian MatrixElements.faf array file. Example gaussian inputs can be found in the test/tests folder.

IOp(3/194) can set the absorbing potential parameters proportional to each atom's van der Waals radius in the following format:
aaabbb01
Where the inner (outer) boundary is set at aaa/10 (bbb/10) times the atom's van der Waals radius, respectively.
IOp(3/195) sets the maximum value of the absorbing potential in atomic units.

TDCI expects both 'MatrixElements.faf' and 'input' to be in the current working directory. Sample tdci input files can be found in the test/tests folder. To do a tdci simulation, simply execute the tdci binary from inside the job directory. Your setup should look something like this:

```bash
~$ cd tdci_h2o
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


## How to cite








