---
title: "TDCI-CAP Manual"
author:
  - "Andrew S. Durden"
  - "H. Bernhard Schlegel"
date: "February 2025"
lang: "en"
toc: true
toc-depth: 3
number-sections: true
---



# Theoretical Background


The time-dependent configuration interaction (TDCI) code can be used to simulate the electronic structure of a molecule interacting with the electric field of an intense, ultra-short laser pulse.
The time-independent, ground and excited configurations of a molecule in the absence of a field are combined with time-dependent linear coefficients to model the molecular wavefunction in the electric field of the laser pulse (the nuclear positions kept fixed).
The interaction of the molecule with the electric field of the laser pulse is treated in the semiclassical dipole approximation and the pulses can be linear or elliptically polarized.
For strong field ionization, the molecule is surrounded by a complex absorbing potential (CAP) to remove the outgoing electron (see Figure 1).

![**Figure 1.** Hydrogen atom wavefunction (light blue), Coulomb potential (dark blue) and complex absorbing potential (dashed dark blue) without a field. The strong field from a laser pulse suppresses the Coulomb potential (red) and the wavefunction (orange) can go over the Coulomb barrier or tunnel through it and be absorbed by the complex potential (dashed red).](img/Figure1.png)

For single ionization, the wavefunction includes the Hartree-Fock ground state and all singly excited configurations (CIS).[^1],[^2]
Ionization of cations can be modeled by a CISD-IP wavefunction which is built from singly ionized and singly excited, singly ionized configurations of the Hartree-Fock ground state of the neutral molecule.
For molecules with heavy elements, spin-orbit coupling (SOC) can be included with an effective one electron SOC operator.[^5]
The wavefunction is propagated with the exponential of the field-dependent Hamiltonian using a Trotter factorization.
The ionization yield for TDCI-CAP simulation is calculated from the decrease in the norm squared of the wavefunction as it interacts with the absorbing potential.
The rate is computed from the expectation value of the absorbing potential using the time-dependent wavefunction.
The integrals needed to construct the Hamiltonian are computed in an electronic structure program such as Gaussian.
For strong field ionization standard molecular basis sets are augmented with several sets of diffuse functions to support the wavefunction as the electron leaves the molecule and interacts with the absorbing potential.
The followings sections provide a brief overview of the theory and computational approaches of the TDCI code.
More details on applications can be found in some of our papers.[^3][^4][^6][^7]

## Hamiltonian

In the TD-CI approach, the electronic wavefunction is propagated with the time-dependent Schrödinger equation:

\begin{equation}
i\hbar \frac{\partial}{\partial t} \Psi(t) \;=\; \hat{H}(t)\,\Psi(t) \;=\; 
\Bigl[\hat{H}_{el} + \hat{V}^{SOC} \;-\; \hat{\vec{\mu}} \cdot \vec{E}(t) \;-\; i \,\hat{V}^{abs}\Bigr] \,\Psi(t)
\tag{1}
\end{equation}

$\hat{H}_{el}$ is the field-free, time-independent non-relativistic electronic Hamiltonian and the interaction with the electric field of an intense laser pulse is treated in the semiclassical dipole approximation, where $\hat{\vec{\mu}}$ is the dipole operator and $\vec{E}$ is the electric field.
The spin-orbit coupling term is approximated by an effective one electron spin-orbit coupling operator. 

\begin{equation}
\hat{V}^{SOC} \;=\; -\,\frac{\alpha_0^2}{2}\,\sum_A Z_A^{eff}\,\frac{(\vec{r} - \vec{r}_A)\,\times\,\nabla}{|\vec{r} - \vec{r}_A|^3}
\tag{2}
\end{equation}

Suitable values for $Z^{eff}$ have been reported by Koseki, Gordon and co-workers and by Chiodo and Russo.
$Z^{eff}$ can be adjusted to reproduce the experimentally observed spin-orbit splitting.
$\hat{V}^{abs}$  is the absorbing potential.

## Absorbing Potential

The absorbing potential, $\hat{V}^{abs}$, surrounds the molecule starting at a distance of about 10 – 20 bohr and removes the outgoing part of the wavefunction, thereby simulating ionization.  
The absorbing potential for the molecule is constructed from spherical potentials centered on each atom.  
The potential on atom C is written in terms of the distance from the nucleus, $r_c = |\vec{r} - \vec{C}|$.  
In the quadratic form, the atomic potential starts at $R_A$, rises quadratically to $\frac{R_A+R_B}{2}$, turns over quadratically to $R_B$, and is constant, $V_{max}$, beyond $R_B$:

\begin{equation}
\hat{V}^{abs}_{C}(r_c) \;=\;
\begin{cases} 
0, & r_c \le R_A,\\[6pt]
V_{max}, & r_c \ge R_B,\\[6pt]
2V_{max}\Bigl(\frac{r_c - R_A}{R_B - R_A}\Bigr)^2, & R_A \le r_c \le \tfrac{R_A + R_B}{2},\\[6pt]
V_{max} \;-\; 2V_{max}\Bigl(\frac{R_B - r_c}{R_B - R_A}\Bigr)^2, & \tfrac{R_A + R_B}{2} \le r_c \le R_B.
\end{cases}
\tag{3}
\end{equation}

In the $\sin^2$ form, the atomic potential on an atom rises as $\sin^2$ between $R_A$ and $R_B$,

\begin{equation}
\hat{V}^{abs}_{C}(r_c) \;=\; 
V_{max}\,\sin^2\Bigl(\tfrac{\pi}{2}\,\tfrac{r_c - R_A}{R_B - R_A}\Bigr), 
\quad R_A \,\le\, r_c \,\le\, R_B.
\tag{4}
\end{equation}

The absorbing potential for the molecule is equal to the minimum of the values of the atomic absorbing potentials. This has also been termed a Voronoi absorbing potential. 

\begin{equation}
\hat{V}^{abs}(\vec{r}) \;=\; 
\min\Bigl(\hat{V}^{abs}_1(r_1),\;\hat{V}^{abs}_2(r_2),\;\dots,\;\hat{V}^{abs}_{\text{Natoms}}(r_{\text{Natoms}})\Bigr).
\tag{5}
\end{equation}

Typical values of the parameters for the absorbing potential are $R_A = 3.5$ times the van der Waals radius for each atom (so that the absorbing potential starts far enough beyond the Coulomb well), $R_B = 10$ times the van der Waals radius (so that the rise in the absorbing potential is gradual enough to minimize any reflection) and $V_{\text{max}} = 10$ hartree (so that the potential is strongly absorbing but the integrals remain finite).  
Figure 2 shows examples of atomic and molecular absorbing potentials.

![**Figure 2.** (a) Shape of an atomic absorbing potential ($\sin^2$ form in red, quadratic form in dashed blue), (b) 2D plot of $R_A$ for absorbing potential of carbon (black), oxygen (red) and hydrogen (grey) in $CH_2O$ (c) 3D plot of $R_A$ (gold) and $\tfrac{R_A+R_B}{2}$ (blue) for the molecular absorbing potential for $CH_2O$.](img/Figure2.png)

## Electric Field of Laser Pulse

The electric field for a linearly polarized $n$-cycle pulse with a $\sin^2$ envelope is:

\begin{equation}
E(t) \;=\;
\begin{cases}
E_{max}\,\sin^2\!\Bigl(\tfrac{\omega t}{2n}\Bigr)\,\cos(\omega t), & 0 \le \omega t \le 2\pi n,\\[6pt]
0, & \omega t > 2\pi n.
\end{cases}
\tag{6}
\end{equation}

For an $n$-cycle circularly polarized 800 nm pulse in the $xz$-plane with a $\sin^2$ envelope:

\begin{equation}
E_x(t) \;=\; E_{max}\,\sin^2\!\Bigl(\tfrac{\omega t}{2n}\Bigr)\,\Bigl[-\,\cos(\omega t)\,\cos(\gamma) \;-\; \sin(\omega t)\,\sin(\gamma)\Bigr],
\end{equation}

\begin{equation}
E_z(t) \;=\; E_{max}\,\sin^2\!\Bigl(\tfrac{\omega t}{2n}\Bigr)\,\Bigl[\cos(\omega t)\,\sin(\gamma) \;-\; \sin(\omega t)\,\cos(\gamma)\Bigr],
\end{equation}

\begin{equation}
E_x(t) \;=\; E_z(t) \;=\; 0,\quad \omega t \ge 2\pi n.
\tag{7}
\end{equation}

$E_{max}$ is the maximum value for the electric field and $\gamma$ is the Euler angle that determines the direction of the major axis of an elliptically polarized pulse.  
To obtain directional information for ionization, a static field can be used instead of an oscillating field.  
To avoid non-adiabatic excitations, the electric field slowly ramps up and turns over to a constant value:

\begin{equation}
E(t) \;=\;
\begin{cases}
E_{max}\Bigl(1 - \Bigl(1 - \tfrac{t}{t_{ramp}}\Bigr)^4\Bigr), & 0 \le t \le t_{ramp},\\[6pt]
E_{max}, & t > t_{ramp}.
\end{cases}
\tag{8}
\end{equation}

The shapes of a 7-cycle linear pulse and a static pulse are shown in Figure 3.  
Several other pulse shapes available in the TDCI code are described in the Appendix.

![**Figure 3.** (a) A 7-cycle linearly polarized 800 nm pulse with a $\sin^2$ envelope (Eq. 6), (b) 2-cycle circularly polarized 800 nm pulse with a $\sin^2$ envelope showing $x$ and $z$ components in green and red, respectively (Eq. 7), and (c) a “static” pulse used to probe the angular dependence of strong field ionization (Eq. 8).](img/Figure3.png)

## Wavefunctions

In the TD-CI approach, the wavefunction is written as a linear combination of time-dependent coefficients with time-independent states of the field-free Hamiltonian, $\Psi_s$:

\begin{equation}
\Psi(t) \;=\; \sum_{s}\,C_{s}(t)\,\Psi_s
\tag{9}
\end{equation}

For single ionization of a neutral molecule to cations, a wavefunction consisting of the ground state and all singly excited configurations (CIS) is suitable,

\begin{equation}
\Psi(t) \;=\; c_0(t)\,\psi_0 \;+\; \sum_{i,a}\,c_i^a(t)\,\psi_i^a
\tag{10}
\end{equation}

where $\psi_i^a$ are singly excited determinants and $i, j$ etc. and $a, b$ etc. are indices for the occupied and virtual molecular spin-orbitals, respectively. For ionization of a cation to a dication, a spin-unrestricted wavefunction could be used, but this treats the $\alpha$ and $\beta$ spin-orbitals differently. This problem can be overcome by using a CISD-IP wavefunction. The wavefunction is constructed using the molecular orbitals of the closed shell system and includes singly ionized determinants, $\psi_x$, and singly excited, singly ionized determinants, $\psi_{xi}^a$:

\begin{equation}
\Psi(t) \;=\; \sum_{x}\,c_x(t)\,\psi_x \;+\; \sum_{x,i,a}\,c_{xi}^a\,\psi_{xi}^a
\tag{11}
\end{equation}

where $x$ indicates the ionized molecular spin-orbitals. The ionizations generate the ground and excited states of the cation, and the single excitations serve to improve the energies of these cation states. As a result, the energies of the cation states will depend slightly on the number of virtual orbitals used to construct the singly excited, singly ionized determinants.

## Propagation of the Wavefunction

For typical strong field ionization studies, numerous simulations of 20 – 50 femtoseconds are needed (e.g., different pulse shapes, intensities and directions, different delay time in pump-probe simulations, etc.). Because of the strong fields and the wide range of the energies of the excited configurations, small time steps are required for the propagation. In the TDCI code the wavefunction is propagated with the exponential of the Hamiltonian, since this unitary transformation is very stable and allows for larger time steps. Typically, a time step of 0.05 au = 1.2 attoseconds can be used. Reducing the time step by a factor of 2 changes the ionization rate by less than 0.02% for typical simulations. The exponential of the Hamiltonian can be obtained via a Trotter factorization with components depending on the field-free Hamiltonian, the absorbing potential, and the dipole moment matrix. 

For a linearly polarized pulse, the propagation for a time step $\Delta t$ is

\begin{equation}
\Psi(t + \Delta t) \;=\; \exp\bigl(-\,i\,\hat{H}\,\Delta t\bigr)\,\Psi(t)
\tag{12}
\end{equation}

\begin{equation}
\begin{aligned}
C(t + \Delta t) \;=\; 
&\exp\Bigl(-\,i\,\hat{H}_{el}\,\tfrac{\Delta t}{2}\Bigr)\,
\exp\Bigl(-\,V^{abs}\,\tfrac{\Delta t}{2}\Bigr)\, \\
&\exp\Bigl(i\,E\bigl(t + \tfrac{\Delta t}{2}\bigr)\,D\,\Delta t\Bigr)\, \\
&\exp\Bigl(-\,V^{abs}\,\tfrac{\Delta t}{2}\Bigr)\,
\exp\Bigl(-\,i\,\hat{H}_{el}\,\tfrac{\Delta t}{2}\Bigr)\,C(t)
\end{aligned}
\tag{13}
\end{equation}

\begin{equation}
\begin{aligned}
C(t + \Delta t) \;=\; 
&\exp\Bigl(-\,i\,\hat{H}_{el}\,\tfrac{\Delta t}{2}\Bigr)\,U^{T}\, \\
&\exp\Bigl(i\,E\bigl(t + \tfrac{\Delta t}{2}\bigr)\,d\,\Delta t\Bigr)\, \\
&U\,\exp\Bigl(-\,i\,\hat{H}_{el}\,\tfrac{\Delta t}{2}\Bigr)\,C(t)
\end{aligned}
\tag{14}
\end{equation}

This requires some initial diagonalizations but avoids the exponentiation of a full matrix at every time step. The field-free Hamiltonian is time-independent and can be diagonalized once at the beginning of the simulation. By working in the eigenbasis of the field-free Hamiltonian, $\exp(-i\,\hat{H}_{el}\,\Delta t / 2)$ is a diagonal matrix and is easy to calculate. Because the absorbing potential is time-independent, $\exp(-V^{abs}\,\Delta t / 2)$ needs to be calculated only once. The calculation of $\exp(i\,\vec{E}(t + \Delta t/2)\,\mathbf{D}\,\Delta t)$ would require an exponentiation of a full matrix at each time step. However, by diagonalizing $\mathbf{d} = \mathbf{W}^T\,\mathbf{D}\,\mathbf{W}$ once at the beginning of the simulation and working in the eigenbasis of $\mathbf{D}$, the contribution reduces to an exponential of a **time-dependent diagonal matrix**, $\exp(i\,\vec{E}(t + \Delta t/2)\,\mathbf{d}\,\Delta t)$. The product $\mathbf{U} = \exp(-V^{abs}\,\Delta t/2)\,\mathbf{W}^T$ is formed once at the beginning of the propagation. Thus, all of the $N^3$ steps need to be done only once at the beginning and can be reused for many subsequent simulations. A propagation step for a linearly polarized pulse with fixed nuclear positions scales as $N^2$ and involves two full matrix-vector multiplies ($\mathbf{U}$ and $\mathbf{U}^T$) and three diagonal matrix-vector multiplies ($\exp(-i\,\hat{H}_{el}\,\Delta t/2)$ and $\exp(i\,\vec{E}(t + \Delta t/2)\,\mathbf{d}\,\Delta t)$).

The corresponding Trotter factorization for a circularly polarized pulse involves two oscillating fields:

\begin{equation}
\begin{aligned}
C(t + \Delta t) = & \exp(-i \hat{H}_{el} \Delta t / 2) \exp(-V^{abs} \Delta t / 2) \\
& \mathbf{W}_2^T \exp\left(i E_z(t + \Delta t / 2) d_2 \Delta t / 2\right) \mathbf{W}_2 \\
& \mathbf{W}_1^T \exp\left(i E_x(t + \Delta t / 2) d_1 \Delta t\right) \mathbf{W}_1 \\
& \mathbf{W}_2^T \exp\left(i E_z(t + \Delta t / 2) d_2 \Delta t / 2\right) \mathbf{W}_2 \\
& \exp(-V^{abs} \Delta t / 2) \exp(-i \hat{H}_{el} \Delta t / 2) C(t).
\end{aligned}
\tag{15}
\end{equation}

Here, $\mathbf{W}_1\,\mathbf{D}_1\,\mathbf{W}_1^T = d_1$ and $\mathbf{W}_2\,\mathbf{D}_2\,\mathbf{W}_2^T = d_2$ are the eigenvalues and eigenvectors of the transition dipole matrices $\mathbf{D}_1$ and $\mathbf{D}_2$ in the two orthogonal field directions. The product $\mathbf{U}' = \mathbf{W}_1\,\mathbf{W}_2$ is formed once at the beginning of the propagation. A propagation step for a circularly polarized pulse with fixed nuclei involves four full matrix-vector multiplies ($\mathbf{U}$, $\mathbf{U}^T$, $\mathbf{U}'$, $\mathbf{U}^{\prime T}$) and five diagonal matrix-vector multiplies.

## Ionization Rate

The ionization rate is taken as the rate of decrease in the norm squared of the wavefunction as the wavefunction interacts with the absorbing potential (Figure 4).

![**Figure 4.** Strong field ionization of hydrogen atom (a) electric field of a “static” pulse with three different intensities (b) dipole moment demonstrating hydrogen responds adiabatically to the field (c) decrease in the norm of the wavefunction as hydrogen is ionized by the electric field.](img/Figure4.png)

Using the time-dependent Schrödinger equation,

\begin{equation}
i\,\frac{\partial\,\Psi(t)}{\partial t} \;=\; \Bigl[\hat{H}(t) \;-\; i\,\hat{V}^{\text{abs}}\Bigr]\,\Psi(t),
\tag{16}
\end{equation}

the rate can be related to the expectation value of the absorbing potential:

\begin{equation}
\begin{aligned}
\text{rate}(t)
&=\; -\,\frac{1}{ \langle \Psi(t)\,|\,\Psi(t)\rangle} \; \frac{d\,\langle \Psi(t)\,|\,\Psi(t)\rangle}{dt} \\[6pt]
&=\; -\,\frac{\Bigl\langle \Psi(t)\,\Big|\, \frac{d\,\Psi(t)}{dt} \Bigr\rangle}{ \langle \Psi(t)\,|\,\Psi(t)\rangle} \;+\;\text{(c.c.)} \\[6pt]
&=\; -\,\frac{\Bigl\langle \Psi(t)\,\Big|\,-\,i\,\bigl[\hat{H}(t) \;-\; i\,\hat{V}^{\text{abs}}\bigr]\,\Psi(t)\Bigr\rangle}{ \langle \Psi(t)\,|\,\Psi(t)\rangle} \;+\;\text{(c.c.)} \\[6pt]
&=\; 2\,\frac{\Bigl\langle \Psi(t)\,\Big|\,\hat{V}^{\text{abs}}\,\Big|\,\Psi(t)\Bigr\rangle}{ \langle \Psi(t)\,|\,\Psi(t)\rangle} \\[4pt]
&=\; 2\,\sum_{r,s}\,C_r^*(t)\,C_s(t)\,\langle \,\Psi_r\,|\,\hat{V}^{\text{abs}}\,|\,\Psi_s\,\rangle,\quad \sum_{r,s}\,C_r^*(t)\,C_s(t)=1 .
\end{aligned}
\tag{17}
\end{equation}




where $r$ and $s$ are indices for configurations. For simplicity, the CIS wavefunction is used in the following equations, but these can be easily extended to CISD-IP. In terms of determinants, the ionization rate is

\begin{equation}
\begin{aligned}
\text{rate}(t)
=\;2\Bigl[\; c_0^*(t)\,c_0(t)\,\langle \Psi_0 | \hat{V}^{\text{abs}} | \Psi_0 \rangle
\;+\;\sum_{i,a}\,c_i^{a*}(t)\,c_0(t)\,\langle \Psi_i^a | \hat{V}^{\text{abs}} | \Psi_0 \rangle \\[6pt]
\quad\quad+\;\sum_{j,b}\,c_0^*(t)\,c_j^b(t)\,\langle \Psi_0 | \hat{V}^{\text{abs}} | \Psi_j^b \rangle
\;+\;\sum_{i,j,a,b}\,c_i^{a*}(t)\,c_j^b(t)\,\langle \Psi_i^a | \hat{V}^{\text{abs}} | \Psi_j^b \rangle
\Bigr].
\end{aligned}
\tag{18}
\end{equation}

The rate can also be written in terms of the density matrix and the absorbing potential in the basis of the molecular orbitals, $\phi_p$:

\begin{equation}
\text{rate}(t)
\;=\; 2\,\sum_{p,q}^{\text{all}}\;\rho_{pq}(t)\;\langle \phi_p\,|\,\hat{V}^{\text{abs}}\,|\,\phi_q \rangle,
\tag{19}
\end{equation}

where $\rho_{pq}(t)$ is the time-dependent one-particle density matrix derived from the CI wavefunction, and indices $p$ and $q$ run over all occupied and virtual orbitals.  
This form is also applicable to rt-TD-DFT.

## Basis Sets

For simulations of strong field ionization with TDCI, a standard molecular basis set must be augmented by several sets of diffuse functions to support the wavefunction as it is distorted by the strong field and interacts with the absorbing potential. Various sets of diffuse functions have been developed and evaluated for their ability to treat strong field ionization. The added diffuse functions should represent the unbound electron density optimally as it interacts with the laser field and the absorbing potential but should be limited in number so that the simulations can be carried out efficiently. These sets include diffuse *s*, *p*, *d*, and *f* gaussian functions with selected even-tempered exponents of the form 0.0001×$2^n$ placed on each atom. When combined with the aug-cc-pVTZ molecular basis set and an absorbing potential starting at 3.5 times the van der Waals radius for each atom, the least diffuse *s*, *p* and *f* functions should have exponents of 0.0256 and 0.0512 for *d* functions. The most diffuse *s*, *p*, *d*, and *f* functions should have exponents of 0.0032, 0.0032, 0.0064 and 0.0064, respectively (or possibly smaller). Because diffuse functions on adjacent centers overlap strongly, the exponents should not be too small since this leads to severe linear dependencies and SCF convergence problems.

![**Figure 5.** Example of diffuse s, p, d and f gaussian basis functions added to support the wavefunction in the region between the Coulomb well and the absorbing potential.](img/Figure5.png)



# Tutorial

## Compilation
TDCI-CAP can be compiled with **nvfortran** on Unix systems supported by NVIDIA's HPC SDK. We have tested only CentOS, Ubuntu, and Archlinux x64 builds.

The NVIDIA HPC SDK can be downloaded at  
<https://developer.nvidia.com/hpc-sdk-downloads>  
NVIDIA's install instructions:  
<https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html>

Once the environment is prepared, you can clone the repository and compile the code:

```bash
git clone https://github.com/hbschlegel/TDCI-CAP
cd TDCI-CAP
make
```

After compilation completes, you will find the `tdci` executable in the `bin/` subdirectory.


## Gaussian Input

Before running TDCI, we must run a preliminary Gaussian Configuration Interaction Singles (CIS) calculation to obtain the necessary integrals. The one-electron integrals are written in the atomic orbital basis and include the dipole moment, spin-orbit, and absorbing potential integrals ($Z_{\text{eff}}$ needs to be specified for each atom for spin-orbit calculations, see Nuclei Props in in the [molecule specification](https://gaussian.com/molspec/) section). The transformed two-electron integrals are written in the molecular orbital basis in the Gaussian “bucket” format. Since most strong field simulations use only a limited number of occupied orbitals and eliminate the highest energy virtual orbitals, the Gaussian input needs to specify the number of occupied and virtual orbitals for the integral transformation. These files can be very large but can be used for multiple simulations as long as the molecular geometry and basis set are the same. 
We use Gaussian development version gdvj28p.  
Below is an example Gaussian input file. Several inputs options are necessary and are described in more detail below:

```
%chk=h2o.chk
%mem=10GB
%nproc=10
#P ucis(mo,nstates=5)/aug-cc-pVDZ extrabasis gfprint Pop=Full int(acc2e=12) 
IOp(3/194=10003501,3/195=10,6/18=1,8/10=91,8/37=2,8/38=-10,9/127=2)
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

`IOp(3/194)` sets the absorbing potential parameters proportional to each atom's van der Waals radius in the following format:
`bbbaaa01` Where the inner boundary ($R_{\text{A}}$) is set at aaa/10 times the atom's van der Waals radius, and the outer boundary ($R_{\text{B}}$) is set at bbb/10 times the atom's van der Waals radius.
For example, the `3/194=10003501` in the example sets $R_{\text{A}}$ to 3.5 radii, and $R_{\text{B}}$ to 10 radii.

`IOp(3/195)` sets the maximum value of the absorbing potential ($V_{\text{max}}$) in atomic units. 10 is suitable.

`IOp(6/18)=1` and `IOp(9/127)=2` are necessary for calculating some of the matrix elements that are written to MatrixElements.faf

`IOp(8/10)=91` enables `IOp(8/37)` and `IOp(8/38)`, which specify the active space.

`IOp(8/37)=N` selects orbital N as the first occupied orbital to be included in the active space. Excitations out of orbitals below will not be considered. Note tIn Gaussian the orbitals are indexed staring at 1, not 0.

`IOp(8/38)=N` If N is positive, selects orbital N as the last orbital to be included in the active space. If N is negative, the N highest energy virtual orbitals will be removed from the active space. In practice, we select this to ignore orbitals with energies above about 3.0 Hartree. 

Don't change the line that starts with `OUTPUT=`; it specifies the indices of the various files to be included in the MatrixElements.faf file.

When selecting an active space, one must be careful to ensure degenerate orbitals remain in the same active or non-active set. If one degenerate orbital is left out of the active space while others are included, it will cause unphysical results.

If you're having convergence difficulties, try combinations of these parameters:
 
> `SCF(QC,NoVarAcc,NoIncFock,MaxCYC=64,VShift=500)`

For more information on the IOps, consult the [Gaussian overlay documentation](https://gaussian.com/iops/). The SCF parameters are documented [here](https://gaussian.com/scf/).

NOTE: If you are using shared computing like WSU's Grid, please remember to follow their rules and policies. This means you should run Gaussian and TDCI from an sbatch script that is scheduled through SLURM. See your computing resource's documentation for more details.

Once your Gaussian calculation successfully finishes, there should be a MatrixElements.faf in the Gaussian job directory.


## TDCI input

TDCI expects both `MatrixElements.faf` and `input` to be in the current working directory. Sample TDCI input files can be found in the test/tests folder, or in the section below. To perform a TDCI simulation, simply execute the tdci binary (or sbatch script) from inside the job directory.
We recommend separating the Gaussian and TDCI job directories, and creating a [symbolic link (symlink)](https://en.wikipedia.org/wiki/Symbolic_link) to the MatrixElements.faf in your TDCI job directory.
Here's how to set up a TDCI job directory in the recommended way. We assume that the Gaussian calculation is performed inside a directory called `gaussian_h2o`.

```bash
mkdir tdci_h2o                                # Create TDCI job directory
cd tdci_h2o                                   # Enter TDCI job directory
ln -s ../gaussian_h2o/MatrixElements.faf .    # Create symlink.
ls -alt                                       # Print out a directory listing
```

Once you finish setting up your input file, your directory should look something like this:

```
drwxr-xr-x 18 user user     4096 Sep 23 16:47 ..
drwxr-xr-x 18 user user     4096 Sep 23 16:47 .
-rw-r--r--  1 user user     2052 Sep 22 18:42 input
lrwxrwxrwx  1 user user       39 Sep 22 18:42 MatrixElements.faf -> ../gaussian_h2o/MatrixElements.faf
```

Then you can run TDCI simulation. If you are running on your own hardware, you can directly invoke the executable:
```
~/TDCI-CAP/bin/tdci
```
If you are using shared computing, make sure to wrap this in an sbatch script. 


<!-- Commenting out this section so we don't have to explain write_matrices = .true.

 TDCI itself exports some data tables in the job directory, but TDCI is also packaged with a python analysis script that generates many tables and plots.
The tools/rate\_analyzer.py script must be executed from the job directory, and passed the path of the gaussian output log as an argument.
```
~/tdci_h2o$ python3 ~/TDCI-CAP/tools/rate_analyzer.py ../gaussian_h2o/h2o.log
```
-->

Next, we provide a sample tdci `input` file.
The `input` file is separated into "namelist" sections that start with `&` and end with `/`.

A table with descriptions of each parameter is provided below.

In this example, we set the initial wavefunction as a superposition between the ground and second excited states with the `init_states` and `init_coeffs` parameters.

The `FIELD` section specifies that a "static" field (slowly ramps up to a static value) will be applied to the system.

The `FIELD_strengths` section specifies that we will only propagate one field strength of 0.0500 atomic units.

The `FIELD_directions` section says that we will apply this field in two different directions, the $(\theta,\phi)=(0,0)$ direction, and the $(\theta,\phi)=(30,0)$ direction.

The parameters in the `SYSTEM` section control the propagation scheme, including the duration and number of timesteps. 


```
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
 envelope = 'stat'
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
 Qmo_dens        = .false.
 /

```

Table of input parameters:


| **Namelist**       | **Parameter**       | **Description**                                                                                                     | **Allowed values**                          | **Notes**                                                                                                                                                        |
|--------------------|---------------------|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DYNAMICS           | init_states(i)      | Init_state(0) is the number of initial states. Init_state(i) for i>0 is the state index of state i.                 | N/A                                         |                                                                                                                                                                   |
| DYNAMICS           | init_coeffs(i)      | Initial coefficient of state i                                                                                      | Complex (e.g., (1.00, 0.00))                | Example: Init_coeffs(1) = (1.0, 0.0)                                                                                                                              |
| DYNAMICS           | restart             | Reads RESTART.bin to restart the propagation                                                                       | .true. or .false.                           |                                                                                                                                                                   |
| DYNAMICS           | Qsave               | Writes RESTART.bin                                                                                                | .true. or .false.                                        |                                                                                                                                                                   |
| FIELD_units        | omega_units         | Units for the omega parameter                                                                                       | ‘au’                                        | ‘au’ angular frequency in atomic units \\ 'nm' wavelength in nanometer                                                                                                                                             |
| FIELD              | dirform             | Coordinate system used for FIELD_directions namelist.                                                              | ‘polar’, ‘cart’                             | ‘polar’: Polar coordinates \\ ‘cart’: Cartesian coordinates                                                                                                       |
| FIELD              | ellipt              | Ellipticity for circularly polarized fields                                                                        | Float from 0.0–1.0                          | A value of 1.0 corresponds to a circular pulse.                                                                                                                   |
| FIELD              | envelope            | Envelope that outlines the pulse                                                                          | ‘none’, ‘cos2’, ‘trap’, ‘gaus’, ‘stat’, ‘band’, ‘ramp’, ‘sin2’, ‘cirl’, ‘cirr’ | For a full description of the options, see the relevant section in documentation.                                              |
| FIELD              | ncyc                | Number of optical cycles that are enclosed by the envelope                                                  | Positive integer                            | Does not apply to fields ‘none’, ‘static’, and ‘ramp’.                                                                                                            |
| FIELD              | omega               | Frequency of the pulse                                                                                   | Positive float                              | Units are set by the omega_units parameter.                                                                                                                       |
| FIELD              | phase               | Carrier envelope phase (CEP) in degrees                                                                             | float                              |                                                                                                                                                                   |
| FIELD              | euler               | Sets the Euler angle (gamma) for elliptical pulses.                                                                | float                          |                                                                                                                                                                   |
| FIELD_strengths    | nemax               | Number of propagations to perform with different field strengths                                                             | Positive integer                            | Read_emax(i) must be set for i=1,…,nemax                                                                                                                           |
| FIELD_strengths    | read_emax(i)        | Maximum field strength of field i.                                                                                 | Non-negative float (i=1,…,nemax)          |                                                                                                                                                                   |
| FIELD_directions   | ndir                | Number of propagations to perform with different field polarization directions.                                    | Positive integer                            | Read_theta(i) and Read_phi(i) must be set for i=1,…,ndir. The total number of propagations will be nemax×ndir.                                                    |
| FIELD_directions   | read_theta(i)       | Theta coordinate for field polarization direction i                                                                | Float from 0–180 (i=1,…,ndir)                           | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | read_phi(i)         | Phi coordinate for field polarization direction i                                                                  | Float from 0–360 (i=1,…,ndir)                           | Requires dirform=’polar’                                                                                                                                         |
| FIELD_directions   | read_x(i)           | x coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | read_y(i)           | y coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| FIELD_directions   | read_z(i)           | z coordinate for field polarization direction i                                                                    | Float                                       | Requires dirform=’cart’                                                                                                                                          |
| SYSTEM_units       | dt_units            | Units for the ‘dt’ parameter.                                                                                      | ‘au’, 'as', 'fs', 'ps'                                        | ‘au’ atomic units \\ 'as' -attosec \\ 'fs' femtosec \\ 'ps' picosec                                                                                                                                              |
| SYSTEM_units       | eigmax_units        | Units for the ‘eigmax’ parameter.                                                                                  | ‘au’, 'eV'                                        | ‘au’ atomic units \\ 'eV' electron volts                                                                                                                                              |
| SYSTEM             | dt                  | Size of the propagation timestep                                                                                   | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | eigmax              | Ignore states with energies above eigmax+ionization                                                                                    | Positive float                              |                                                                                                                                                                   |
| SYSTEM             | ionization          | estimated ionization energy                                                                          | Float                                       | if ionization = -1.0, use -HOMO energy                                                                                                                                                                  |
| SYSTEM             | jobtype             | Selects the type of propagation to be performed                                                                    | ‘cis’, ‘cisip’, ‘soc’, ‘socip’              | ‘cis’: TD-CIS \\ ‘cisip’: TD-CISDIP \\ ‘soc’: TD-CIS with spin-orbit coupling \\ ‘socip’: TD-CISDIP with spin-orbit coupling                                   |
| SYSTEM             | nstep               | Number of timesteps to be propagated                                                                               | Positive integer                            |                                                                                                                                                                   |
| SYSTEM             | outstep             | Output data will only be written on step numbers divisible by outstep                                              | Positive integer                            | If nstep=50 and outstep=10, data will be written on steps 10, 20, 30, 40, and 50.                                                                                 |
| SYSTEM             | nactive             | Number of active occupied orbitals (lower orbitals are ignored)                                                    | positive integer                                         | nactive = 0 to use all occupied orbitals                                                                                                                                                                  |
| SYSTEM             | nvirtual            | Number of active virtual orbitals (higher orbitals are ignored)                                                    | positive integer                                         | nvirtual = 0 to use all virtual orbitals                                                                                                                                                                  |
| SYSTEM             | socfac              | Spin-orbit coupling scale factor                                                                  | float                                        |  socfac = 1.0 for full SOC \\ socfac = 0.0 to turn off SOC                                                                                                                                                                 |
| SYSTEM             | socfacZ             | Scales z-component of spin-orbit coupling                                                       | N/A                                         |                                                                                                                                                                   |
| SYSTEM             | IP_alpha_beta       |  Which electrons to ionize                                                                                                                   | .true. or .false.                           | .true. ionize both alpha and beta \\ .false. ionize only beta electrons                                                                                                                                                                  |
| SYSTEM             | QeigenDC            | use divide-and-conquer eigenvalue solver                                                                                                                    | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_tdcidata      | Type of Gaussian data file to read                                                                                 | .true. or .false.                           | For legacy TDCI.dat files, use .true. \\ For new MatrixElements.faf files, use .false.                                                                            |
| InOutFILES         | tdcidatfile         | Specifies the filename of the Gaussian data file                                                                   | String                                      |                                                                                                                                                                   |
| InOutFILES         | outputfile          | Name of output log file                                                                                            | String                                      |                                                                                                                                                                   |
| InOutFILES         | restartbinfile      | Name of the restart file                                                                                           | String                                      | Setting parameter ‘restart’ to .true. will read from this file. \\ Setting parameter ‘Qsave’ to .true. will write this file.                                     |
| InOutFILES         | Qwrite_ion_coeff    | For sequential double ionization, writes the ion coefficients                                                      | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qread_ion_coeff     | For sequential double ionization, reads the ion coefficients                                                       | .true. or .false.                           |                                                                                                                                                                   |
| InOutFILES         | Qmo_dens            | Writes the density matrix to MO_density files                                                                      | .true. or .false.                           |                                                                                                                                                                   |


## Output Files

This section describes the main output files generated by the TDCI code.
Each file (aside from OUTPUT) is named according to the electric-field strength index (X) and direction index (Y), typically appearing in the filename as -eX-dY.
For example, the RESULTS file for the progataion with the first specified field and direction will be named RESULTS-e1-d1.dat


### OUTPUT
This is the main log file. Similar to a Gaussian log file, it contains both data tables and information about the progress of the calculatiAfter running a TDCI calculation, this should be the first place you check to make sure that it completed successfully.
The parameters from 'input' and various control parameters are listed at the beginning of the file.
OUTPUT also lists the molecule specification, dipole moment, number of basis functions, occupied orbitals, vitrual orbitals and states, specifications for the laser pulse, propagation time step,  simulation time, etc.
For 'cis' or 'soc' calculations, energies and coefficients of the ionized state and the neutral state are listed next (for 'ip' or 'socip', the energies and coefficients are for the doubly ionized state and the singly ionized state).
The coefficients of the CI determinants are C(i,a) for 'cis' and 'soc', and C(i) and C(i,j,a) for 'ip' and 'socip'.
The indices are negative for alpha electrons, positive for beta electrons.
Abs(i) runs from 1 to the number of occupied orbitals; abs(a) runs from 1 to the number of virtual orbitals.
A summary is printed for each individual propagation and includes the maximum field strength and direction, the final norm and dipole moment of the wavefunction.
In a good ionization simulation, you should aim for a final norm between 0.3 and 0.6 in the direction with most ionization.
You may have to play around with the Emax to achieve this.

### RESULTS-eX-dY.dat

The RESULTS file tabulates values of interest at each timestep (the first 3 lines describe the contents).
It provides time, field, total norm, ionization rate, and dipole moment components (`mu_x,y,z`).
The MO99 number describes the total number of orbitals needed to contribute 99% of the total ionization rate. NO99 is the same for a user-supplied set of orbitals.
The rate data must be divided by norm2 to get the rate for a normalized wavefunction.

### POP-eX-dY.dat

Tracks the populations of occupied orbitals at each timestep (the first 2 lines describe the contents).
Also contains the ionization rate for the normalized wavefunction partitioned into contributions from each occupied or virtual orbital, as described in equation [^6].

### ION-eX-dY.dat

Tracks the populations of orbitals and states in the ionized system (the first line describe the contents).



[^1]: [Krause, P.; Sonk, J. A.; Schlegel, H. B., Strong field ionization rates simulated with time-dependent configuration interaction and an absorbing potential. J. Chem. Phys. 2014, 140, 174113. 10.1063/1.4874156](https://doi.org/10.1063/1.4874156) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/369.pdf))
[^2]: [Krause, P.; Schlegel, H. B., Angle-dependent ionization of small molecules by time-dependent configuration interaction and an absorbing potential. J. Phys. Chem. Lett. 2015, 6, 2140-2146. 10.1021/acs.jpclett.5b00929](https://doi.org/10.1021/acs.jpclett.5b00929) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/381.pdf))
[^3]: [Hoerner, P.; Schlegel, H. B., Angular dependence of strong field ionization of CH3X (X = F, Cl, Br, or I) using time-dependent configuration interaction with an absorbing potential. J. Phys. Chem. A 2017, 121, 5940-5946. 10.1021/acs.jpca.7b06108](https://doi.org/10.1021/acs.jpca.7b06108) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/409.pdf))
[^4]: [Hoerner, P.; Schlegel, H. B., Angular dependence of strong field ionization of haloacetylenes, HCCX (X = F, Cl, Br, I) using time-dependent configuration interaction with an absorbing potential. J. Phys. Chem. C 2018, 122, 13751–13757. 10.1021/acs.jpcc.8b00619](https://doi.org/10.1021/acs.jpcc.8b00619) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/416.pdf))
[^5]: [Lee, M. K.; Hoerner, P.; Li, W.; Schlegel, H. B., Effect of spin-orbit coupling on strong field ionization simulated with time-dependent configuration interaction. J. Chem. Phys. 2020, 153. 10.1063/5.0034807](https://doi.org/10.1063/5.0034807) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/435.pdf))
[^6]: [Lee, M. K.; Li, W.; Schlegel, H. B., Angular dependence of strong field sequential double ionization for neon and acetylene simulated with time-dependent configuration interaction using CIS and CISD-IP. J. Chem. Phys. 2020, 152, 064106. 10.1063/1.5133659](https://doi.org/10.1063/1.5133659) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/429.pdf))
[^7]: [Schlegel, H. B.; Hoerner, P.; Li, W., Ionization of HCCI neutral and cations by strong laser fields simulated with time dependent configuration interaction. Front. Chem. 2022, 10, 866137. 10.3389/fchem.2022.866137](https://doi.org/10.3389/fchem.2022.866137) -- ([Mirror](https://schlegelgroup.wayne.edu/Pub_folder/445.pdf))

