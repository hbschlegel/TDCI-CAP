#!/usr/bin/env python3


#
# Usage: Run from a completed TDCI directory, argument is relative or absolute path to corresponding gaussian output
# This script partitions the rate by atomic orbital, and can generate plots
#   showing contribution to the rate based on atom/quantum numbers. 
#

import math
import numpy as np
import struct, sys, re, os, csv, random
import scipy
from scipy.optimize import curve_fit
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(formatter={'float': lambda x: "{0:0.4e}".format(x)})


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

import pickle
from collections import defaultdict

#from matplotlib.colors import ListedColormap, BoundaryNorm
#from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# globals
au2fs = 0.0241888432650900
fs2au = (1/au2fs)

# WARNING: BIN FILES ARE NOT PORTABLE BETWEEN COMPUTERS
#          WITH DIFFERENT DOUBLE IMPLEMENTATIONS  
def read_bin_array(filepath, length):
  #print("l: "+str(length))
  if length == 0:
    return np.array([])
  f = open(filepath, 'rb')
  return np.array(struct.unpack('d'*length, f.read()))

def write_bin_array(array, filepath):
  if (len(np.shape(array)) != 1):
    print("WARNING: write_bin_array() is for 1D double arrays! Your array has shape:")
    #print(np.shape(array))
  f = open(filepath, 'wb')
  f.write( struct.pack('d'*len(array), *array ))
  f.close(); del f

# Parses the file 'input' for the theta/phi of each direction
def parse_tdci_input():
  polar_coords = []
  f = open('input', 'r')
  lines = []
  l = f.readline()
  lines.append(l)
  while l.strip() != "&FIELD_directions":
    l = f.readline()
    lines.append(l)
    if l == "":
      print("Can't find &FIELD_directions !!!")
  l = f.readline()
  lines.append(l)
  nemax = int(l.split()[-1])
  l = f.readline()
  lines.append(l)
  while l.strip() != "/":
    if l == "":
      print("Can't find end of directions !!!")
      sys.exit()
    index, theta, phi = parse_tdci_input_line(l)
    polar_coords.append( (theta,phi) )
    l = f.readline()
    lines.append(l)
  #print("Polar coords:")
  #print( (nemax, len(polar_coords)) )
  #print(polar_coords)

  while l != "": # Finish reading file
    l = f.readline()
    lines.append(l)

  nsteps = 0
  nprintstep = 0
  dt = 0.0 # Timestep size in au
  Emax = -99.0 
  for line in lines:
    if len(line.split()) > 1:
      #print(line.split()[0])
      if line.split()[0] == 'nstep':
        nsteps = int(line.split()[-1])
      if line.split()[0] == 'outstep':
        nprintstep = int(line.split()[-1])
      if line.split()[0] == 'dt':
        dt = float(line.split()[-1])
      if line.split()[0] == "read_emax(1)":
        Emax = float(line.split()[-1])
  assert nsteps != 0
  assert nprintstep != 0
  assert dt != 0.0
  assert Emax != -99.0
      

  return polar_coords, nsteps, nprintstep, dt, Emax
  
  
# Get ndirs and the theta,phi for each direction
def parse_tdci_input_line(line):
  # read_theta( 1) = 0.d0    ;  read_phi( 1) = 0.d0 
  # regex magic collects first index and the params, d/D case insensitive, ignore spaces.
  pattern = r'read_theta\((\d+)\)\s*=\s*([\d.]+[dD][\d+-]+)\s*;\s*read_phi\((\d+)\)\s*=\s*([\d.]+[dD][\d+-]+)'
  match = re.search(pattern, line.replace(' ', ''))
  if match:
    index = int(match.group(1))
    theta = float(match.group(2).replace('d', 'e').replace('D', 'e'))
    phi = float(match.group(4).replace('d', 'e').replace('D', 'e'))
    return index, theta, phi
  print("NO MATCH ON LINE!!: "+str(line))
  sys.exit()


def get_word_from_last_line(file_path, number):
  with open(file_path, 'r') as f:
    return f.readlines()[-1].split()[number]


def matrix2csv(mat,n, outname):
  f = open(outname,'w')
  i, j = 0,0
  while i<n:
    j=0
    while j<n:
      f.write(f"{mat[i][j]:.10E},")
      j+=1
    f.write("\n")
    i+=1
  f.close()
  return 0

def matrix2csvRECT(mat,n,m, outname):
  f = open(outname,'w')
  i, j = 0,0
  while i<n:
    j=0
    while j<m:
      f.write(f"{mat[i][j]:.10E},")
      j+=1
    f.write("\n")
    i+=1
  f.close()
  return 0
    

# Accepts OUTPUT_FIELD_SHAPE path and returns Efield(t)
def parse_tdci_field(filename):
  f = open(filename, 'r')
  f.readline()
  l = f.readline()
  Efield = []
  while l != "":
    Efield.append(float(l.split()[-1]))
    l = f.readline()
  return Efield



# Get CMO and some parameters from an fchk file
def parse_fchk(filename):
  nea = 0
  nbasis = 0
  nbsuse = 0
  cmo = []

  f = open(filename,'r')
  keystr = "Number of alpha electrons"
  line = f.readline()
  while line != "":
    ls = line.split()
    if "Number of alpha electrons" in line:
      nea = int(ls[-1])
    if "Number of basis functions" in line:
      nbasis = int(ls[-1])
    if "Number of independent functions" in line:
      nbsuse = int(ls[-1])
    if "Alpha MO coefficients" in line:
      N = int(ls[-1]) # number of elements in fchk section
      assert N == nbasis*nbsuse
      line = f.readline()
      while len(cmo) < N:
        ls = line.split()
        #print(ls)
        for entry in ls:
          #print(entry)
          cmo.append(float(entry))
        line = f.readline()
      assert len(cmo) == N

    line = f.readline() 
  return nea, nbasis, nbsuse, np.array(cmo)
  

# Get AO density from population analysis of gaussian output
# For sanity testing MO2AO transformation
def parse_gdv_AO_density(filename,nao):
  f = open(filename, 'r')
  line = f.readline()
  # This thing is triangular
  density_AO_tri = []
  while ("Density Matrix" not in line):
    line = f.readline()
  while (("Beta Density Matrix" not in line) and ("Mulliken" not in line)):
    ls = line.split()
    for entry in ls:
      if "." in entry: density_AO_tri.append(float(entry))
    line = f.readline()

  f.close()
  density_AO_tri = np.array(density_AO_tri)
  # Unpack triangular matrix, arranged in columns of 5
  density_AO_2D = np.zeros((nao, nao))
  idx = 0
  BLOCK_SIZE = 5
  for block_start in range(0, nao, BLOCK_SIZE):
      for i in range(block_start, nao):
          for j in range(block_start, min(i+1, block_start+BLOCK_SIZE)):
              density_AO_2D[i, j] = density_AO_tri[idx]
              density_AO_2D[j, i] = density_AO_tri[idx]  # Symmetry
              idx += 1
  return density_AO_2D*2

def mm_MO2AO(nbasis, norb, density_MO, CMO):
  #print("start MO2AO")
  density_MO.resize( (norb,norb) )
  CMO.resize((norb,nbasis))
  mm = np.matmul
  density_AO = mm( CMO.T, mm( density_MO, CMO )  )
  return density_AO
 

def MO2AO_sanity(ref_density_AO, CMO, nbasis, norb, nocc):
  # Make ground state MO density matrix
  density_MO = np.zeros((norb,norb))
  for i in range(0,nocc):
    density_MO[i,i] = 2.0

  # Transform to AO
  #density_AO = tdci_MO2AO(nbasis, norb, density_MO, CMO)
  density_AO = mm_MO2AO(nbasis, norb, density_MO, CMO)

  # Compare with reference
  print("Comparing density_AO with ref_density_AO")
  good = True
  density_AO.resize( (nbasis*nbasis) )
  ref_density_AO.resize( (nbasis*nbasis) )
  for i in range(0,nbasis**2):
  #for i in range(0,1000):
    if (np.abs(density_AO[i]-ref_density_AO[i]) > 1e-3):
      good = False
      print(f"BAD ({i//nbasis},{i%nbasis}) : {density_AO[i]:.3E} , {ref_density_AO[i]:.3E}")
  if good: print("Passed!!")
  else: print("Failed! :( ")

  

class orb:
  def __init__(self):
    self.orb_index = 0
    self.atom_index = 0
    self.atom_type = 0
    self.n = 0
    self.l = 0
    self.m = 0
    self.prim_expon = []
    self.prim_coeff = []

class gdvlog_parser:
  def __init__(self, filename):
    self.gaussdir = os.path.dirname(filename)+"/"
    self.filename = filename 
    self.nao, self.nmo, self.natoms, self.nocc, self.Ra, self.Rb = 0, 0, 0, 0, None, None
    self.get_params()
    #self.run_MO2AO_sanity( self.nao, self.nmo)
    self.orbs = []
    # What should I do with 'SP' orbitals? I think just treat them as P?
    self.l_map = {'S': 0, 'P': 1, 'SP': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}
    self.inv_l_map = ['S', 'P', 'D', 'F', 'G', 'H', 'I'] # Only used for labels
    # the 'SP' entry breaks this
    #self.inv_l_map = {value: key for key, value in self.l_map.items()}
    self.atomtypes = {} # map iatom -> atom type (H, He, etc.)


  # Something doesn't work in this test either... not sure
  def run_MO2AO_sanity(self, nao, nmo):
    CMO = read_bin_array( "matrices/CMO.bin", nmo*nao)
    CMO.resize((nmo,nao)) # This is correct.
    # MO2AO density transformation sanity check
    print("getting reference density from gdv log")
    ref_density_AO = parse_gdv_AO_density(self.filename, nao)
    print("doing sanity check")
    MO2AO_sanity(ref_density_AO, CMO, nao, nmo, self.nocc)
    print("finished")
    

  # Parses TDCI output
  def get_params(self):
    nao = 0
    nmo = 0
    natoms = 0
    nocc = 0

    # Parse tdci output for array sizes
    f = open("OUTPUT" , 'r')
    line = f.readline()
    while line != "":
      ls = line.split()
      if len(ls) > 2:
        if ls[0] == "noa":
          nocc = int(ls[2])
        #     charge = 0  multiplicity = 1  natoms = 4
        if ls[0] == "charge":
          natoms = int(ls[8])
        #      nbasis  = 390
        if ls[0] == "nbasis":
          nao = int(ls[2])
        #      nrorb   = 267
        if ls[0] == "nrorb":
          nmo = int(ls[2])
          
      line = f.readline()

    #print(f"nao, nmo : {nao}, {nmo}")
    self.nao, self.nmo, self.natoms, self.nocc = nao, nmo, natoms, nocc
    #self.Ra = Ra
    #self.Rb = Rb
    #assert len(Ra) == natoms
    #return nao, nmo, natoms, nocc, Ra, Rb
    return 0

  # Return tuple (n, l, m)
  # Examples:
  # "1S" -> (1, 0, 0)
  # "9PX" -> (9, 1, -1)
  # "20D+2" -> (20, 2, 2)
  # "29F 0" -> (29, 3, 0) 
  def decipher_orbstring(self,orb_string):
    l_map = self.l_map

    # get rid of space and lowercase
    orb_string = orb_string.replace(" ", "").upper()
    l_letter = ''.join(filter(str.isalpha, orb_string))[0]
    l = l_map[l_letter]
    lstr_index = orb_string.index(l_letter)
    n = int( orb_string[:lstr_index] )

    m=0
    # Decode m
    if 'S' in orb_string:
      m = 0
    elif 'X' in orb_string: # for PX,PY,PZ
      m = -1
    elif 'Y' in orb_string:
      m = 1
    elif 'Z' in orb_string:
      m = 0
    else: # for D, F, and up
      m = int( orb_string[lstr_index+1:] ) 

    return (n,l,m)


  # NOTE atom_index and orb_index are 1-indexed!!
  def parse_vabs(self):
    f = open(self.filename, 'r')
    keystr = "Absorbing boundary integrals in AO basis"
    line = f.readline()
    while (keystr not in line) and (line != ""):
      line = f.readline()
      
    # Next 3 lines are header info
    f.readline() ; f.readline() ; f.readline()
    # Read first data line
    line = f.readline().split()
    # If a line starts a new atom, it should have 6 words
    # 180 3   H  1S          0.00000   0.00000
    # Otherwise, it should have 4 OR 5
    # 154       20D+2        0.12572   1.79229
    # 131       16D 0        0.00000   0.00000
    # (that space in the orbstring is annoying...)
    # Every atom should start with 1S, so we shouldn't have to
    # worry about about strings like "131  5 O 16D 0 0.000 0.000" etc with 7 words
    # Final line
    # RMS (Vabsrob - Sum Vabsorb(i)) =    0.13469D+00

    atom_index = 0
    atom_type = ""
    while line[0][0].isnumeric():
      orb_index = 0
      orbstring = ""
      vabs_diag = 0.0
      vabs_rowsum = 0.0
      neworb = orb()

      orb_index = int(line[0])

      if len(line) == 6: # New atom
        atom_index = int(line[1])
        atom_type = line[2]
        self.atomtypes[atom_index] = atom_type
        orbstring = line[3] # Never will have space like '16D 0'
      elif len(line) == 4: # 154       20D+2        0.12572   1.79229
        orbstring = line[1]
      elif len(line) == 5: # 131       16D 0        0.00000   0.00000
        orbstring = line[1]+line[2]

      vabs_diag = float(line[-2])
      vabs_rowsum = float(line[-1])
        
      # Save info
      neworb.orb_index = orb_index
      neworb.atom_index = atom_index
      neworb.atom_type = atom_type
      neworb.n, neworb.l, neworb.m = self.decipher_orbstring(orbstring)
      neworb.vabs_diag = vabs_diag
      neworb.vabs_rowsum = vabs_rowsum
      self.orbs.append(neworb) # Append to output data
      # End loop, make sure next line doesn't crash the while loop
      line = f.readline().split()
      if len(line)==0:
        break
      elif len(line[0])==0:
        break
    #print("Done Parsing!")

  def parse_gdv_params(self):
    self.nea = 0
    self.nbasis = 0
    self.nbsuse = 0
    Ra = []
    Rb = []
    f = open(self.filename, 'r')
    line = f.readline()
    while line != "":
      ls = line.split()
      if "NBasis=" in line:
        self.nbasis = int( line.split()[1] )
      if "NBsUse=" in line:
        self.nbsuse = int( line.split()[1] )
      if "alpha electrons" in line:
        self.nea = int( line.split()[0] )
      if "IAtom   RA        RB" in line:
        #print("FOUND VABS RA RB")
        line = f.readline()
        ls = line.split()
        # in my jobs the line after this section is:
        # There are a total of     1322692 grid points.
        while (len(ls) == 3):
          Ra.append(float(ls[1]))
          Rb.append(float(ls[2]))
          line = f.readline()
          ls = line.split()
      line = f.readline()
    self.Ra = Ra
    self.Rb = Rb
    #print(Ra)


    
    

  # This should be executed after parse(), and it assumes that the list self.orbs
  #   is properly in order ( orbs[i] == orbs[i].orb_index ), as it should be
  def parse_prims(self):
    #print("Entered parse_prims")
    f = open(self.filename, 'r')
    keystr = "AO basis set (Overlap normalization)"
    line = f.readline()
    while (keystr not in line) and (line != ""):
      line = f.readline()

    # Read first data line -- new shell
    line = f.readline()
    # First line after data looks like 
    #   185 basis functions,   237 primitive gaussians,   216 cartesian basis functions
    # in h2o, in hcci it looks like:
    #   There are   206 symmetry adapted cartesian basis functions of A1  symmetry.
    # So for our keystring we could use 'basis' or 'function'

    # Data lines will look like either:
    #  Atom C3       Shell    62 S   7     bf  212 -   212          0.000000000000          0.000000000000         -4.998683
    # or
    #  0.8236000000D+04  0.5419783203D-03
    # (exponent, coefficient) 

    # Initialize shell variables in outer scope
    atomi_ = atomtype_ = n_ = nprim = bfi_min = bfi_max = None
    xcoord = ycoord = zcoord = None
    prim_i = 0

    while ( 'basis' not in line ):
      ls = line.split()
      #print(ls[0])
      if ls[0] == "Atom": # New shell
        # Mostly extra info used for asserts
        #atomi_ = int(ls[1][1:])
        atomi_ = int(''.join(filter(lambda x: x.isdigit(), ls[1]))) # Remove non-numbers for He31 etc.
        atomtype_ = ''.join(filter(lambda x: not x.isdigit(), ls[1])) # Remove numbers for He31 etc.
        #n_ = int(ls[3]) # What actually is this value? Shell number?
        l_ = self.l_map[ls[4]] 
        nprim = int(ls[5]) # Not needed
        bfi_min = int(ls[7]) # INCLUSIVE RANGE, 1-indexed
        bfi_max = int(ls[9])
        xcoord = float(ls[10]) # I don't think we need these, just putting them here incase we do later
        ycoord = float(ls[11])
        zcoord = float(ls[12])
        #print(ls)
        #print(f"l_, nprim, bfi_min, bfi_max: {(l_, nprim, bfi_min, bfi_max)}")

      if ( (len(ls) == 2) or (len(ls) == 3) ): # New prim
        #print(f"bfi_min, bfi_max+1: {(bfi_min, bfi_max+1)}")
        for i in range(bfi_min, bfi_max+1):
          # Make sure nothing is fishy
          # Remember orbs is 0-indexed while orb_index, bfi are 1-indexed.
          assert i == self.orbs[i-1].orb_index , f"Orb index mismatch, {i} != {self.orbs[i-1].orb_index}" 
          assert( atomtype_ == self.orbs[i-1].atom_type )
          assert( atomi_ == self.orbs[i-1].atom_index )
          # Actually, the following assert will fire for SP orbitals, since l_map['SP'] = 1, 
          #  but they contain an S primitive.
          #assert l_ == self.orbs[i-1].l, f"Orb {i}, l mismatch: {l_} != {self.orbs[i-1].l}" 
          # Add the prims
          #   Python requires scientific notation have an E instead of D
          self.orbs[i-1].prim_expon.append(float(ls[0].replace('D','E'))) 
          coeffs = tuple(float(coeff.replace('D', 'E')) for coeff in ls[1:])
          self.orbs[i-1].prim_coeff.append(coeffs)

      line = f.readline()
    for idx in range(0,len(self.orbs)):
      if (len(self.orbs[idx].prim_coeff)==0):
        print("BAD ORBITAL!!")
        print(str(idx)+"th orbital:")
        print("orb_index  : "+str(self.orbs[idx].orb_index ) )
        print("atom_index : "+str(self.orbs[idx].atom_index ))
        print("atom_type  : "+str(self.orbs[idx].atom_type ) )
        print("n          : "+str(self.orbs[idx].n ) )
        print("l          : "+str(self.orbs[idx].l ) )
        print("m          : "+str(self.orbs[idx].m ) )
        print("prim_expon : "+str(self.orbs[idx].prim_expon ))
        print("prim_coeff : "+str(self.orbs[idx].prim_coeff ))
        

  
    # debug print to test
    if False:
      idx = 78
      for idx in range(0,len(self.orbs)):
        print(str(idx)+"th orbital:")
        print("orb_index  : "+str(self.orbs[idx].orb_index ) )
        print("atom_index : "+str(self.orbs[idx].atom_index ))
        print("atom_type  : "+str(self.orbs[idx].atom_type ) )
        print("n          : "+str(self.orbs[idx].n ) )
        print("l          : "+str(self.orbs[idx].l ) )
        print("m          : "+str(self.orbs[idx].m ) )
        print("prim_expon : "+str(self.orbs[idx].prim_expon ))
        print("prim_coeff : "+str(self.orbs[idx].prim_coeff ))
      return 0

    

##############################
# End gdvlog_parser class
##############################

class ratechecker:
  def __init__(self, gdvparser):
    self.gdvparser = gdvparser
    self.gaussdir = gdvparser.gaussdir
    self.nao, self.nmo = gdvparser.nao, gdvparser.nmo
    self.natoms, self.nocc = gdvparser.natoms, gdvparser.nocc
    self.ntri = self.nao*(self.nao+1)//2
    self.CMO = None
    self.Vabs_AO = None
    self.prep_matrices()


  # Read in and transform 
  def prep_matrices(self):
    mm = np.matmul
    nao, nmo, ntri = self.nao, self.nmo, self.ntri
    
    vabs_tri = read_bin_array( "matrices/Vabs_AO.bin", ntri )


    CMO = read_bin_array( "matrices/CMO.bin", nmo*nao)
    CMO.resize((nmo,nao))
    self.CMO = CMO
    matrix2csvRECT(CMO,nmo,nao, "CMO.csv")

    # Unpack triangular matrix
    vabs_AO = np.zeros((nao,nao))

    for iao in range(0,nao):
      for jao in range(0,nao):
        ij = (iao+1)*(iao)//2 + jao
        if (jao>iao) :
          ij = (jao+1)*(jao)//2 + iao
        vabs_AO[iao][jao] = vabs_tri[ij]

    test_vabs_MO = mm( CMO, mm( vabs_AO, CMO.T))

    self.Vabs_AO = vabs_AO
    self.Vabs_MO = test_vabs_MO
    matrix2csv(vabs_AO, nao, f"Vabs_AO.csv")

    return 0

  # Expect timestep and d to be 1-indexed
  def make_vdens(self,timestep,e,d,writecsv=False,prevtime=0):
    nao, nmo = self.nao, self.nmo
    CMO = self.CMO
    dens_MO = None
    if os.path.isfile(f"matrices/MO_density-e{e}-d{d}.{timestep}.bin"):
      dens_MO = read_bin_array( f"matrices/MO_density-e{e}-d{d}.{timestep}.bin", nmo**2)
    else: # dirty hack for missing files
      dens_MO = read_bin_array( f"matrices/MO_density-e{e}-d{d}.{prevtime}.bin", nmo**2)
    dens_MO.resize((nmo,nmo))
    mm = np.matmul
    dens_AO = mm(CMO.T, mm(dens_MO, CMO))

    # Test trace
    #MOtr = np.trace(dens_MO)
    #AOtr = np.trace(dens_AO)
    #print(f"Dir{d}, Step {timestep}, MOtr={MOtr:.3f}, AOtr={AOtr:.3f}")

    #CMO_pinv = self.CMO_pinv
    #dens_AO = CMO_pinv @ dens_MO @ CMO_pinv.T
    if writecsv:
      matrix2csv(dens_AO, nao, f"DensityAO-d{d}-{timestep}.csv")

    dens_AO.resize((nao**2))
    self.Vabs_AO.resize((nao**2))
    vdens = np.multiply(self.Vabs_AO,dens_AO) # Element-wise multiplication
    vdens = vdens*2 # density is for alpha electrons, so *2 for full rate
    vdens.resize((nao,nao))
    if timestep == 1:
      sumrate = np.sum(vdens)
      tracerate = 0.0
      for i in range(0,nao):
        tracerate += vdens[i][i]
      #print(f"vdens dir {d} from AO : sumrate: {sumrate} | {2*sumrate/au2fs}")
      #print(f"vdens dir {d} from AO : tracerate: {tracerate}")

      # Test MO rate
      vabs_MO = read_bin_array( "matrices/Vabs_MO.bin", nmo**2)
      dens_MO.resize((nmo**2))
      vdens_MO = vabs_MO * dens_MO
      vdens_MO.resize((nmo,nmo))
      sumrate = np.sum(vdens_MO)
      tracerate = 0.0
      for i in range(0,nmo):
        tracerate += vdens_MO[i][i]
      #print(f"vdens dir {d} from MO : sumrate: {sumrate} | {2*sumrate/au2fs}")
      #print(f"vdens dir {d} from MO : tracerate: {tracerate}")

    #print("in make_vdens")
    #print(vdens[:2,:5])
    #print(np.sum(vdens))
    return vdens


class RateTrajectoryData:
  def __init__(self, nsteps, nprintstep, natoms, lmax, norbs, direction):
    self.npts = nsteps//nprintstep
    npts = self.npts
    self.nsteps = nsteps
    self.rate = np.zeros(npts)
    self.norm2 = np.zeros(npts)
    self.atomsum = np.zeros((npts,natoms))
    self.lsum = np.zeros((npts,natoms,lmax))
    self.totsum = np.zeros(npts) # Sanity check for rate
    self.orbrate = np.zeros((npts,norbs))
    self.maxrate = 0
    self.direction = direction
    self.set_norm2()
    self.expdict_final = {} # At final timestep.

  def set_norm2(self):
    f = open(f"RESULTS-e1-d{self.direction+1}.dat", 'r')
    f.readline() ; f.readline() ; f.readline()
    line = f.readline() #first data line, step 1
    for i in range(self.npts):
      self.norm2[i] = float( line.split()[9] ) 
      line = f.readline()
    f.close()
    return 0


class RatePlotter:
  def  __init__(self, gdvlog, plotname=None ):
    if plotname is None:
      try:
        plotname =  os.path.dirname(gdvlog).replace("_"," ") # Attempt at autogenerating name
      except:
        plotname = ""
    self.plotname = plotname
    self.gaussdir = os.path.dirname(gdvlog)+"/"
    self.parser = gdvlog_parser(gdvlog)
    parser = self.parser
    parser.parse_vabs()
    parser.parse_prims()
    parser.parse_gdv_params()

    self.polar_coords, self.nsteps, self.nprintstep, self.dt, self.Emax = parse_tdci_input()
    self.Efield = parse_tdci_field("OUTPUT_FIELD_SHAPE")
    self.dataset = [] # RateTrajectoryData indexed by direction index.
    self.ndir = len(self.polar_coords)
    print("ndir: "+str(self.ndir))

    self.nao, self.nmo, self.natoms, self.nocc = parser.nao, parser.nmo, parser.natoms, parser.nocc
    nao, nmo, natoms = self.nao, self.nmo, self.natoms
    self.lmax = 4
    nsteps = self.nsteps
    lmax = self.lmax
    self.npts = self.nsteps//self.nprintstep
   
    e_idx = 1 # Field number
    d_idx = 1 # Direction

    self.rc = ratechecker(self.parser)
    expdict_max = {}
    # Generate index keys for angular and exponent breakdown
    for orb in parser.orbs:
      if len(orb.prim_expon)>0:
        l_ = orb.l
        exp_ = orb.prim_expon[0]
        # Maximum magnitude across directions
        expdict_max[ (orb.atom_index, l_, exp_) ] = 0.0

    for d in range(0,self.ndir):
      d_data = RateTrajectoryData(self.nsteps, self.nprintstep, self.natoms, self.lmax, len(parser.orbs), d)
      rate = d_data.rate
      atomsum = d_data.atomsum
      lsum = d_data.lsum
      totsum = d_data.totsum
      orbrate = d_data.orbrate
      norm2 = d_data.norm2
      expdict_final = d_data.expdict_final
      expdict = {}


      # Generate index keys for angular and exponent breakdown
      for orb in parser.orbs:
        if len(orb.prim_expon)>0:
          l_ = orb.l
          exp_ = orb.prim_expon[0]
          expdict[ (orb.atom_index, l_, exp_) ] = np.zeros(self.npts)

      # Generate density*vabs in AO basis and partition by atom and L
      # for nstep=20 and nprintstep=5, [5,10,15,20]
      prevtime = 0
      for time in range(self.nprintstep, self.nsteps+self.nprintstep, self.nprintstep):
        writecsv = False
        if time in [self.nprintstep, self.nsteps]: writecsv = True
        vdens = self.rc.make_vdens(time,e_idx,d+1, writecsv, prevtime)
        #print("in __init__")
        #print(vdens[:2,:5])
        #print(np.sum(vdens))
        #rate[time] = np.trace(vdens)
        t_idx_max = (self.nsteps//self.nprintstep)-1
        t_idx = (time//self.nprintstep)-1
        rate[t_idx] = np.sum(vdens) # sum all elements
        if (t_idx == t_idx_max ) or (t_idx == 0):
          #print(f"rate for step={time} d={d+1} : {rate[t_idx]} au | {(1/au2fs)*rate[t_idx]:.10f} fs")
          matrix2csv(vdens, nao, f"vdens-d{d+1}-{time}.csv")
        # Atoms/orbitals are 1-indexed!
        #print("Breakdown of rate by atom and angular momentum")

        for orb in parser.orbs:
          idx = orb.orb_index-1
          #orbrate[time][idx] = vdens[idx][idx]
          orbrate[t_idx][idx] = np.sum(vdens[idx])
          if len(orb.prim_expon)>0:
            key = (orb.atom_index, orb.l, orb.prim_expon[0])
            #expdict[ (orb.atom_index, orb.l, orb.prim_expon[0]) ][time] += vdens[idx][idx]
            tmpsum = np.sum(vdens[idx])
            expdict[ key ][t_idx] += tmpsum
            if abs(expdict_max[key]) < abs(tmpsum):
              expdict_max[key] = np.sum(vdens[idx])
            expdict_max[key] = max( abs(expdict_max[key]), abs(np.sum(vdens[idx])) )
        if (t_idx == t_idx_max):
          with open(f"expdict_d{d+1}-{time}.pickle", "wb") as f:
            #pickle.dump(expdict_max, f, fix_imports=True, buffer_callback=None )
            pickle.dump(expdict_max, f, fix_imports=True )
          tempdict = defaultdict(float)
          for (atomidx, l_, exp_), values in expdict.items():
            tempdict[(l_, exp_)] += values[-1] # Only final timestep, sum over atomidx
          d_data.expdict_final = dict(tempdict)
          


        for iatom in range(0,self.natoms):
          for l_ in range(0, lmax):
            for orb in parser.orbs:
              if orb.l == l_ and orb.atom_index == iatom+1:
                #print(orb.orb_index)
                #print(len(a.orbs))
                idx = orb.orb_index-1 # orb_index is 1-indexed
                #lsum[time][iatom][l_] += vdens[idx][idx]
                #atomsum[time][iatom] += vdens[idx][idx]
                lsum[t_idx][iatom][l_] += np.sum(vdens[idx])
                atomsum[t_idx][iatom]  += np.sum(vdens[idx])
          totsum[t_idx] += atomsum[t_idx][iatom]
        if np.abs(totsum[t_idx]-rate[t_idx])>1e-3:
          print(f"totsum doesnt match rate: totsum={totsum[t_idx]} , rate={rate[t_idx]}")
          pass
        prevtime = time # dirty hack for missing files
      self.dataset.append(d_data)
      self.rate_by_time_plots(d_data)
      #self.individual_orb_plot(d_data, parser)
      self.AtomLExp_plot(d_data, parser, expdict)
      #self.RateCheck_CSV(d_data, parser)
      self.AtomLExp_AvgCSV(d_data, parser, expdict)
      self.RmaxRaRbPlot(expdict, ADDLABELS=True, DIRSTRING=str(d_data.direction+1))
      self.RmaxRaRbPlot(expdict, ADDLABELS=False, DIRSTRING=str(d_data.direction+1))
    self.polar_cmap()
    self.polar_angular_plot()
    #combine_csv_files(self.ndir)
    self.AtomLExp_DirAvgCSV()
    print(f"min(expdict_max) : {min(expdict_max)}")
    with open(f"expdict_max.pickle", "wb") as f:
      #pickle.dump(expdict_max, f, fix_imports=True, buffer_callback=None )
      pickle.dump(expdict_max, f, fix_imports=True )
    self.RmaxRaRbPlot(expdict_max)
    self.RmaxRaRbPlot(expdict_max, ADDLABELS=False)


  def RmaxRaRbPlot(self, expdict_max, ADDLABELS=True, DIRSTRING=None):
    #from adjustText import adjust_text
    parser = self.parser
    Ra = self.parser.Ra
    Rb = self.parser.Rb

    exponent_list = []
    atype_list = []
    for key, val in expdict_max.items():
      if ((key[2] < 0.1) and key[2] not in exponent_list):
        exponent_list.append(key[2])
      atype = parser.atomtypes[key[0]]
      if atype not in atype_list: atype_list.append(atype)
    exponent_list.sort(reverse=True) # sort our exponents descending

    atoms_plotted = []
    #print("Ra:")
    #print(Ra)
    for i in range(1,len(Ra)+1):
      if parser.atomtypes[i] in atoms_plotted:
        continue
      else:
        atoms_plotted.append(parser.atomtypes[i])
      #print(f"plotting {parser.atomtypes[i]}!")
      fig = plt.figure()
      ax = fig.add_subplot(1,1,1)

      # Plot the absorbing potential
      ax2 = ax.twinx()
      X_abs = np.linspace(0,50,500)
      Y_abs = np.zeros(len(X_abs))
      for idx_, x_ in enumerate(X_abs):
        if x_ < Ra[i-1]:
          Y_abs[idx_] = 0.0
        elif x_ > Rb[i-1]:
          Y_abs[idx_] = 1.0
        else:
          Y_abs[idx_] = np.sin( (np.pi/2.)*(X_abs[idx_]-Ra[i-1])/(Rb[i-1]-Ra[i-1])  )**2
      ax2.plot(X_abs, Y_abs, color="black", linestyle="--", label=r"$V_{abs}(r)$")
      ax2.set_ylim(0, 1.1)
      ax2.get_yaxis().set_visible(False)

      # Create data and plot for each angular momentum
      texts = []
      X_tot = []
      Y_tot = []
      label_tot = []
      for l_ in range(0,self.lmax):
        l_str = parser.inv_l_map[l_] 
        X = [] # distance
        Y = [] # rate contribution
        exponent = [] # exponent for label
        for exp_ in exponent_list:
          if ( i, l_, exp_ ) not in expdict_max.keys():
            #print(f"missing: {(i,l_,exp_)}")
            continue # Skip
          #print(f"hit: {(i,l_,exp_)}")
          Rmax = np.sqrt((l_+1)/exp_)
          X.append(Rmax)
          X_tot.append(Rmax)
          if DIRSTRING is None:
            Y.append( (1/au2fs)*expdict_max[ (i,l_,exp_) ] )
            Y_tot.append( (1/au2fs)*expdict_max[ (i,l_,exp_) ] )
          else: # Using directional expdict which has an extra time index
            Y.append( abs((1/au2fs)*expdict_max[ (i,l_,exp_) ][-1]) )
            Y_tot.append( abs((1/au2fs)*expdict_max[ (i,l_,exp_) ][-1]) )
          exponent.append( exp_ )
          idxlist = []
          for orb in parser.orbs:
            if ( (i==orb.atom_index) and (l_==orb.l) and (exp_==orb.prim_expon[0]) ):
              idxlist.append( orb.orb_index )
          label_tot.append( str(f"{exp_:.4f}") )
        #ax.plot(X,Y, marker='o', label=l_str, markersize=14 )
        ax.plot(X,Y, marker='o', label=l_str )
      if ADDLABELS:
        for j, (x, y, label_) in enumerate(zip(X_tot, Y_tot, label_tot)):
          ax.annotate(label_, (x,y), fontsize=8, xytext=(0,0),
                      textcoords='offset points', ha='center', va='bottom')
      # Go to the next 5-bohr mark after the last basis function, or Rb.
      xlim_max = max( math.ceil((max(X_tot)+1)/5)*5  ,  math.ceil((Rb[i-1]+1)/5)*5 )
      ax.set_xlim(0.0,xlim_max)

      lines_1, labels_1 = ax.get_legend_handles_labels()
      lines_2, labels_2 = ax2.get_legend_handles_labels()
      ax.legend(lines_1 + lines_2, labels_1 + labels_2)

      #adjust_text(texts, x=X_tot, y=Y_tot, horizontalalignment='center', autoalign='xy', arrowprops=dict(arrowstyle='->', color='red'))
        
      # translate directory naming convention
      names = (os.getcwd()).split('_')[-2:]
      basisname = ""
      paramname = ""
      for name in names:
        if "basis" in name:
          basisname = name
        else:
          paramname = name

      #ax.title.set_text(f"Rmax for {parser.atomtypes[i]}, {basisname}, {paramname}")
      if DIRSTRING is None:
        ax.title.set_text(f"Rmax for {parser.atomtypes[i]}, {self.plotname}")
      else:
        ax.title.set_text(f"Rmax for {parser.atomtypes[i]}, d{DIRSTRING}, {self.plotname}")
      ax.set_ylabel("Max Rate Contribution ($fs^{-1}$)")
      ax.set_xlabel("Rmax Distance (bohr)")
      ax.set_ylim( 0.0, np.max(Y_tot)*1.05 )

      # Vertical lines for Ra and Rb
      ax.axvline(x=Ra[i-1], color='r', linestyle='--')
      #ypos = ax.get_ylim()[0] - 0.015 * (ax.get_ylim()[1] - ax.get_ylim()[0])
      #ax.text(Ra[i-1], -0.03, "Ra", va='top', ha='center', transform=ax.get_xaxis_transform())
      ax.axvline(x=Rb[i-1], color='r', linestyle='--')
      #ax.text(Rb[i-1], -0.03, "Rb", va='top', ha='center', transform=ax.get_xaxis_transform())

      # This is causing a very weird crash...
      #plt.tight_layout()

      if ADDLABELS:
        if DIRSTRING is None:
          plt.savefig(f"Rmax{parser.atomtypes[i]}.png", dpi=400)
        else:
          plt.savefig(f"Rmax{parser.atomtypes[i]}-d{DIRSTRING}.png", dpi=400)
      else:
        if DIRSTRING is None:
          plt.savefig(f"Rmax{parser.atomtypes[i]}_nolabels.png", dpi=400)
        else:
          plt.savefig(f"Rmax{parser.atomtypes[i]}-d{DIRSTRING}_nolabels.png", dpi=400)
      plt.close()
      
    
  # Compare rate calculated from AO density * Vabs_AO
  #   with the rate in RESULTS.dat files
  # Note: the RESULTS value is multiplied by norm^2
  def RateCheck_CSV(self,d_data, parser):
    f = open(f"RESULTS-e1-d{d_data.direction+1}.dat", 'r')
    g = open(f"ratecheck-d{d_data.direction+1}.csv", 'w')
    line = f.readline() ; f.readline() ; f.readline()
    line = f.readline() #first data line, step 1
    rate = d_data.rate
    g.write("step,rate python,rate RESULTS\n")
    for i in range(self.nsteps//self.nprintstep):
      RESULTSrate = float( line.split()[10] )
      g.write(f"{i},{(1/au2fs)*rate[i]:.10f},{RESULTSrate:.10f}\n")
      line = f.readline()
    f.close()
    g.close()
    return 0


  def AtomLExp_plot(self,d_data, parser, expdict):
    orbs = parser.orbs
    npts = self.npts
    norbs = len(orbs)
    #maxrates = np.max(d_data.orbrate, axis=0)
    maxrates = np.zeros( len(expdict) )
    keys = []
    data = []
    for key, val in expdict.items():
      keys.append(key)
      data.append(val)
    data = np.array(data)
    maxrates = np.max(data, axis=1) # Max across time axis
    #print("maxrates:")
    #print(maxrates)
    orb_idx_ratesorted = np.argsort(maxrates)[::-1]
    MAX_LINES = 10 # max number of lines to plot

    X = np.array(list(range(1,npts+1)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)


    for i in range(0,MAX_LINES):
      idx = orb_idx_ratesorted[i]
      key = keys[idx]
      l_str = parser.inv_l_map[key[1]]
      atype = parser.atomtypes[key[0] ] 
      ax.plot( X, (1/au2fs)*data[idx],
        label=f"{atype}{key[0]}, {l_str}, {key[2]}, {(1/au2fs)*maxrates[idx]:0.3f}")
        #label=f"{atype}{key[0]}, {l_str}, {key[2]}, {(1/au2fs)*maxrates[idx]:0.3f}")
    ax.legend(title=f"Atom, L, Exp, MaxRate", loc='upper right')
    ax.set_ylabel("Rate $fs^{-1}$")
    ax.set_xlabel("Timestep")
    ax.title.set_text(f"Top 10 orbitals, Dir: {d_data.direction+1}\n{self.plotname}")
    
    plt.tight_layout()
    plt.savefig(f"AtomLExp{d_data.direction+1}.png", dpi=200)
    plt.close()


    # Write CSV
    f = open(f"AtomLExp_rate{d_data.direction+1}.csv", 'w')
    f.write(f"Direction,{d_data.direction+1},,,\n")
    f.write("AtomIdx,AtomType,OrbIdxList,L,Exponent,MaxRate(fs^-1)\n")
    for i in range(0,len(maxrates)):
      idx = orb_idx_ratesorted[i]
      # key is ( atomidx, l, exponent )
      key = keys[idx]
      atype = parser.atomtypes[ key[0] ]
      l_str = parser.inv_l_map[key[1]]
      idxlist = []
      for orb in parser.orbs:
        #print(key)
        #print(orb.prim_expon)
        if ( (key[0]==orb.atom_index) and (key[1]==orb.l) and (key[2]==orb.prim_expon[0]) ):
          idxlist.append( orb.orb_index )
      idxliststr = str(idxlist).replace(" ","").replace("[","").replace("]","").replace(","," ")
      
      outstr = f"{key[0]},{atype},{idxliststr},{l_str},{key[2]},{(1/au2fs)*maxrates[idx]}\n"
      #print(outstr)
      f.write(outstr)
    f.close()

    # Write CSV
    f = open(f"AtomLExp_FinalRate{d_data.direction+1}.csv", 'w')
    f.write(f"Direction,{d_data.direction+1},,FinalNorm2,{d_data.norm2[-1]},\n")
    f.write("AtomIdx,AtomType,OrbIdxList,L,Exponent,FinalRate(fs^-1),Norm2 Adj. FinalRate(fs^-1)\n")
    for i in range(0,len(maxrates)):
      idx = orb_idx_ratesorted[i]
      # key is ( atomidx, l, exponent )
      key = keys[idx]
      atype = parser.atomtypes[ key[0] ]
      l_str = parser.inv_l_map[key[1]]
      idxlist = []
      for orb in parser.orbs:
        #print(key)
        #print(orb.prim_expon)
        if ( (key[0]==orb.atom_index) and (key[1]==orb.l) and (key[2]==orb.prim_expon[0]) ):
          idxlist.append( orb.orb_index )
      idxliststr = str(idxlist).replace(" ","").replace("[","").replace("]","").replace(","," ")
      finalrate = expdict[key][-1]
      finalrate_norm2 = expdict[key][-1]/d_data.norm2[-1]
      outstr = f"{key[0]},{atype},{idxliststr},{l_str},{key[2]},{(1/au2fs)*finalrate:.10E},{(1/au2fs)*finalrate_norm2:.10E}\n"
      #print(outstr)
      f.write(outstr)
    f.close()




  # CSV Summary Average Table
  def AtomLExp_AvgCSV(self, d_data, parser, expdict): 
    f = open(f"AtomLExp_AvgRate{d_data.direction+1}.csv", 'w')
    exponent_list = []
    atype_list = []
    for key, val in expdict.items():
      if ((key[2] < 0.1) and key[2] not in exponent_list):
        exponent_list.append(key[2])
      atype = parser.atomtypes[key[0]]
      if atype not in atype_list: atype_list.append(atype)
    exponent_list.sort(reverse=True) # sort our exponents descending

      
    #header = f"orbital,exponent{','+f for f in atype_list}"
    header = "orbital,exponent"
    for i in range(0,len(atype_list)):
      header = header + f",{atype_list[i]}(fs-1)"
    f.write(header)
    for l_ in range(0,self.lmax):
      l_str = parser.inv_l_map[l_]
      for exp_ in exponent_list:
        entries = np.zeros(len(atype_list))
        nowrite = False
        for iatype, atype in enumerate(atype_list):
          #print( (iatype, atype) )
          denom = 0 # for averaging, count atoms of each type
          for iatom in range(1,self.natoms+1): # four for loops!
            if parser.atomtypes[iatom] == atype:
              if ( iatom, l_, exp_ ) in expdict.keys():
                #print(iatype)
                #print( expdict[ (iatom,l_,exp_) ] )
                entries[iatype] += max(expdict[ (iatom,l_,exp_) ]) # max across timesteps
                denom += 1
              else: # WARNING: if different atoms have different basis, this needs to be changed
                nowrite = True # Some l_/exp_ pairs do not exist.
              #  print("warning: no "+str( (iatom,l_,exp_) ))
          if denom > 1:
            entries[iatype] = entries[iatype]/float(denom) # Average over same typed atoms
        #outstr = "\n"+f"{l_str},{exp_}{','+str(entries[i]) for i in range(0,len(atype_list))}"
        outstr = "\n"+f"{l_str},{exp_}"
        for i in range(0,len(atype_list)):
          outstr = outstr + f",{(1/au2fs)*entries[i]}"
        if nowrite: # Skip basis elements that do not exist
          pass
        else:
          f.write(outstr)
    f.close()
        
  # Combines the CSV files 
  def AtomLExp_DirAvgCSV(self):
    csv_files = [f"AtomLExp_AvgRate{d}.csv" for d in range(1,self.ndir+1)]

    
    # Read the first file to get dimensions
    with open(csv_files[0], 'r') as f:
      reader = csv.reader(f)
      headers = next(reader)
      first_file_data = list(reader)
      num_rows = len(first_file_data)
      num_cols = len(headers)

    # Initialize numpy arrays
    data_sum = np.zeros((num_rows, num_cols - 2))  # Exclude first two columns
    row_count = np.zeros(num_rows)

    # Process all CSV files
    for fname in csv_files:
      with open(fname, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for i, row in enumerate(reader):
          data_sum[i] += np.array([float(val) for val in row[2:]])
          row_count[i] += 1

    # Calculate averages
    averages = data_sum / row_count[:, np.newaxis]

    # Prepare output data
    output_data = [headers]
    with open(csv_files[-1], 'r') as f:
      reader = csv.reader(f)
      next(reader)  # Skip header
      for i, row in enumerate(reader):
        new_row = row[:2] + [f"{avg:.8e}" for avg in averages[i]]
        output_data.append(new_row)

    # Write output to new CSV file
    with open("RateBasis_Summary.csv", 'w', newline='') as f:
      
      writer = csv.writer(f)
      writer.writerows(output_data)



  def individual_orb_plot(self,d_data, parser):
    orbs = parser.orbs
    npts = self.npts
    norbs = len(orbs)
    maxrates = np.max(d_data.orbrate, axis=0)
    orb_idx_ratesorted = np.argsort(maxrates)[::-1]
    MAX_LINES = 10 # max number of lines to plot

    X = np.array(list(range(1,self.npts+1)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(0,MAX_LINES):
      idx = orb_idx_ratesorted[i]
      orb = orbs[idx]
      l_str = parser.inv_l_map[orb.l]
      ax.plot( X, (1/au2fs)*d_data.orbrate[:,idx],
        label=f"{orb.atom_type}, {l_str}, {orb.prim_expon}, {(1/au2fs)*maxrates[idx]:0.3f}")
    ax.legend(title=f"Atom, L, Exp, MaxRate")
    ax.title.set_text(f"Top 10 orbitals, Dir: {d_data.direction+1}\n{self.plotname}")
    
    plt.tight_layout()
    plt.savefig(f"orbsort{d_data.direction+1}.png", dpi=200)
    plt.close()
      
      
      

  def rate_by_time_plots(self,d_data):
    npts = self.npts
    lmax = self.lmax
    parser = self.parser
    natoms = self.natoms
    maxrate = 1e-5 # forgot where to get this from
    for iatom in range(0,natoms):
      for l_ in range(0,lmax):
        maxrate = max( max(d_data.lsum[:,iatom,l_]), maxrate)
    maxrate = maxrate

    #######################################
    # Plot by Atom and Angular momentum
    #######################################
    X = (self.dt*au2fs*self.nprintstep)*np.array(list(range(1,self.npts+1)))
    #print(f"X[0], X[-1]: {(X[0], X[-1])}")
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for iatom in range(0,natoms):
      for l_ in range(0, lmax):
        if max(d_data.lsum[:,iatom,l_]) < 0.10*maxrate:
          continue # Skip low contributions to avoid clutter
        l_str = parser.inv_l_map[l_]
        #print(parser.atomtypes, flush=True)
        atomstr = str(iatom+1)+parser.atomtypes[iatom+1]
        Y = fs2au * (d_data.lsum[:, iatom, l_] - d_data.lsum[0, iatom, l_])
        ax.plot( X , Y, label=f"{atomstr}: {l_str}" )

    # Plot total rate
    ax.plot(X, fs2au*d_data.rate, color='black', label="Total Rate")

    ax.legend()
    ax.title.set_text(f"Rate by Atom:Angular Momentum. Dir: {d_data.direction+1}\n{self.plotname}")
    ax.set_xlabel("Time (fs)")
    
    plt.tight_layout()
    plt.savefig(f"atomangular{d_data.direction+1}.png", dpi=300)
    plt.close()

    ###########################################
    # Plot by angular momentum, atoms summed
    ###########################################
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    l_noatom = np.sum(d_data.lsum, axis=1) # Now has shape l_noatom[nsteps][lmax]
    for l_ in range(0, lmax):
      #if max(l_noatom[:,l_]) < 0.1*maxrate:
      #  continue # Skip low contributions to avoid clutter
      l_str = parser.inv_l_map[l_]
      Y = (1/au2fs)* (l_noatom[:, l_])
      ax.plot( X , Y, label=f"{l_str}" )

    # Plot total rate
    ax.plot(X, (1/au2fs)*d_data.rate, color='black', label="Total Rate")

    ax.legend()
    ax.title.set_text(f"Rate by Angular Momentum. Dir:{d_data.direction+1}\n{self.plotname}")
    
    plt.tight_layout()
    plt.savefig(f"angular{d_data.direction+1}.png", dpi=300)
    plt.close()

    # CSV angular momentum, no atoms
    f = open(f"angular{d_data.direction+1}.csv",'w')
    f.write("step,S,P,D,F,G,H\n")
    for i, time in enumerate(X):
      writestr = f"{time}"
      for l_ in range(0, lmax):
        #Y = (1/au2fs)*l_noatom[time-1, l_]/d_data.norm2[time-1]
        Y = (1/au2fs)* (l_noatom[i, l_])
        writestr += f",{Y:.12E}"
      f.write(writestr+'\n')
    f.close()


    #######################
    # Figure 1            #
    #######################
    plt.clf()
    fig, ax1 = plt.subplots()

    Y_Efield = self.Emax * np.array( self.Efield[:-1:self.nprintstep] )
    Y_rate = fs2au * (d_data.rate - d_data.rate[0]) 

    for l_ in range(0, lmax):
      l_str = parser.inv_l_map[l_]
      Y = fs2au* (l_noatom[:, l_] - l_noatom[0, l_]) # Subtract time=0 tunneling value
      ax1.plot(X, Y, label=f"{l_str}")

    # Plot total rate
    ax1.plot(X, Y_rate, color='black', label="Total Rate")
    #print(f"max rate : {max(d_data.rate)}")

    ax2 = ax1.twinx() # second y-axis sharing x
    #print(f"size X, Efield: {len(X), len(Y_Efield)}")
    ax2.plot(X, Y_Efield, color='black', label=r"$\overrightarrow{E}$", linestyle='dashdot', linewidth=2)
    #print(f"max Efield: {max(Y_Efield)}")

    ax3 = ax1.twinx()
    #ax3.spines['right'].set_position(('outward', 60))
    #ax3.spines['right'].set_visible(False)
    #ax3.spines['left'].set_visible(False)
    ax3.get_yaxis().set_visible(False)
    
    exp_field_model = lambda Efield_, scale_, K_: K_*np.exp(Efield_ / scale_)
    # Optimize scaling factor. p0 is initial guess.
    #popt, _ = curve_fit(exp_field_model, Y_Efield, Y_rate, p0=[0.0035]) 
    scale = 0.0035
    K = d_data.rate[-1]/np.exp(Y_Efield[-1]/scale)
    try:
      popt, _ = curve_fit(exp_field_model, Y_Efield, Y_rate, p0=[0.0035, d_data.rate[-1]]) 
      scale, K = popt
    except:
      pass
       
    #scale = popt[0]
    exp_field = (d_data.rate[-1]/(np.exp(Y_Efield[-1]/scale)))*np.exp(Y_Efield / scale)
    #print(f"Figure1 Alpha: {(d_data.rate[-1]/(np.exp(Y_Efield[-1]/scale)))}")
    #print(f"max exp_field: {max(exp_field)}, scale: {scale}, scaling(K): {K}")
    
    #ax3.plot(X, exp_field, color='red', linestyle=':', label=f"{K:0.2f}"+"exp($\overrightarrow{E}$/"+f"{scale:0.4f})", linewidth=2)
    ax3.plot(X, exp_field, color='red', linestyle=':', label="exp($\overrightarrow{E}$/"+f"{scale:0.4f})", linewidth=2)


    # Add legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    lines_3, labels_3 = ax3.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2 + lines_3, labels_1 + labels_2 + labels_3)

    ax1.set_ylabel("Rate Contribution ($fs^{-1}$)")
    ax1.set_xlabel("Time (fs)")
    ax1.title.set_text(f"{self.plotname}, d{d_data.direction+1}")

    plt.tight_layout()
    plt.savefig(f"Fig1-d{d_data.direction+1}.png", dpi=300)
    plt.close()

    # End Figure 1.



    ########################################
    # Plot by atom, no angular momentum
    ########################################
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for iatom in range(0,natoms):
      #if max(atomsum[:,iatom]) < 0.1*maxrate:
      #  continue # Skip low contributions to avoid clutter
      #print(parser.atomtypes, flush=True)
      atomstr = str(iatom+1)+parser.atomtypes[iatom+1]
      ax.plot( X , fs2au*(d_data.atomsum[:, iatom] - d_data.atomsum[0, iatom] ), label=f"{atomstr}" )

    # Plot total rate
    ax.plot(X, fs2au*d_data.rate, color='black', label="Total Rate")

    ax.legend()
    ax.title.set_text(f"{self.plotname} Dir:{d_data.direction+1}")
    
    plt.tight_layout()
    plt.savefig(f"atom{d_data.direction+1}.png", dpi=300)
    plt.close()

    return 0 


  def polar_cmap(self):

    lmax = self.lmax
    ndir = len(self.polar_coords)
    X_theta = []
    ratemin = 99.99 # cmap windowing
    ratemax = -99.99
    for x in self.polar_coords:
      #X_theta.append( x[0] * (np.pi/180.)  )
      X_theta.append( x[0] )
    Rate_L = [ [] for _ in range(lmax+1) ]
    Rmax_L = [ [] for _ in range(lmax+1) ]
    X_L = [ [] for _ in range(lmax+1) ]
    for d in range(0,ndir):
      ddat = self.dataset[d]
      #print("items ddat")
      #print(ddat.expdict.items())
      tempdict = ddat.expdict_final.copy()

      # Generate Rmax and data for this direction
      for (l_, exp_), value in tempdict.items():
        if value > ratemax: ratemax = value # cmap windowing
        if value < ratemin: ratemin = value 
        Rate_L[l_].append( value )
        Rmax_L[l_].append( np.sqrt((l_+1)/exp_) )
        X_L[l_].append(X_theta[d])

    # Make plot
    colors = ['blue','white','red']
    #cmap = ListedColormap(colors)
    #bounds = [ratemin, 0, ratemax]
    #norm = BoundaryNorm(bounds, cmap.N)

    bounds = [ratemin, 0, ratemax]
    #cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    #norm = TwoSlopeNorm(vmin=ratemin, vcenter=0, vmax=ratemax)
    

    print(f"ratemin, ratemax : {(ratemin, ratemax)}")
    plt.clf()
    fig, axes = plt.subplots(1,lmax, figsize=(20,7), sharey=True)
    for l in range(0,lmax):
      l_label = self.parser.inv_l_map[l]
      sc = axes[l].scatter(X_L[l], Rmax_L[l], c=Rate_L[l], cmap="Greys", vmin=ratemin, vmax=ratemax)
      #sc = axes[l].scatter(X_L[l], Rmax_L[l], c=Rate_L[l], cmap=cmap, norm=norm, edgecolor='black')
      axes[l].set_title(l_label)

    fig.subplots_adjust(right=0.85, wspace=0.05)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.7])
    fig.supxlabel(r"Direction ($^\circ$)", fontsize=18)
    #cbar = fig.colorbar(sc, ax=axes, orientation='vertical')
    cbar = fig.colorbar(sc, ax=axes, cax=cbar_ax, orientation='vertical')
    plt.suptitle(f"{self.plotname}", fontsize=18)
    #plt.tight_layout(rect=[0,0,0.95,1])
    #plt.tight_layout()
    plt.savefig("polar_cmap.png",dpi=400)
 

    
  

  def polar_angular_plot(self):
    # For now we just ignore phi i think
    #print( (len(self.dataset), len(self.polar_coords)))


    lmax = 0
    for orb in self.parser.orbs:
      if orb.l > lmax: lmax = orb.l
    print(f"lmax: {lmax}, {self.parser.inv_l_map[lmax]}")

    ndir = len(self.polar_coords)
    X_theta = []
    for x in self.polar_coords:
      X_theta.append( x[0] * (np.pi/180.)  )
    Y_L = [ [] for _ in range(lmax+1) ]
    Y_L_norm = [ [] for _ in range(lmax+1) ]
    Y_tot = []
    Y_tot_norm = []
    
    nva99 = []

    #ymin = -0.15
    #ymax = 0.12
    ymin = 0.0
    ymax = 0.001

    time = (self.npts-1)
    for d in range(0,ndir):
      ddat = self.dataset[d]
      nva99.append(get_word_from_last_line(f"RESULTS-e1-d{d+1}.dat",3))
      Y_tot.append( ddat.rate[time] )
      #Y_tot_norm.append( ddat.rate[time]/ddat.norm2[time] )
      #print( (d, ddat.rate[time], ddat.norm2[time], ddat.rate[time]/ddat.norm2[time] ) )
      for l_ in range(0,lmax+1):
        rate_l = 0.0
        #rate_l_norm = 0.0
        for iatom in range(0,self.natoms):
          rate_l += ddat.lsum[time][iatom][l_]
          #rate_l_norm += ddat.lsum[time][iatom][l_]/ddat.norm2[time]
        Y_L[l_].append(rate_l)
        #Y_L_norm[l_].append(rate_l_norm)
        if (1./au2fs)*rate_l>ymax: ymax = (1./au2fs)*rate_l

    Y_tot_ = (1/au2fs)*np.array(Y_tot)
    #print(f"Y_tot_norm before au2fs")
    #print(Y_tot_norm)
    #Y_tot_norm = (1/au2fs)*np.array(Y_tot_norm)
    #print(f"Y_tot_norm after au2fs")
    #print(Y_tot_norm)

    ymin, ymax = 0.0, max(Y_tot_)
    #ymin, ymax = 0.0, max(Y_tot_norm)

    # Do flat plot.
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( X_theta, Y_tot_, label=f"Total Rate")
    ax.set_ylabel("Norm Adj. Rate")
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"Direction ($^\circ$)")
    ax.set_title(f"Angular momentum contribution to rate by incident light angle\n"+self.plotname)
    plt.savefig("polarflat_totrate.png", dpi=600)

    # Do flat angular momenum plot
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for l_ in range(0,lmax+1):
      l_str = self.parser.inv_l_map[l_]
      Y = Y_L[l_]
      #Y = Y_L_norm[l_]
      Y = (1/au2fs)*np.array(Y)
      ax.plot(X_theta, Y, label=f"{l_str}")
    ax.plot( X_theta, Y_tot_, label=f"Total Rate")
    #ax.plot( X_theta, Y_tot_norm, label=f"Total Rate")
    ax.set_ylabel("Norm Adj. Rate")
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"Direction ($^\circ$)")
    ax.legend()
    ax.set_title(f"Total rate by incident light angle\n"+self.plotname)
    plt.savefig("polarflat_totrate.png", dpi=600)
    


    plt.clf()
    fig, ax = plt.subplots( subplot_kw={'projection':'polar'})
    for l_ in range(0,lmax+1):
      l_str = self.parser.inv_l_map[l_]
      #print(X_theta)
      #print(Y_L[l_])
      #ax.plot(X_theta, Y_L[l_], label=f"{l_str}") 
      # 180 degree symmetry
      X_theta_neg = list(-1.0*np.array(X_theta))
      X = X_theta + (X_theta_neg[1:-1][::-1])
      X = np.array(X)
      Y = Y_L[l_] + (Y_L[l_][1:-1][::-1])
      #Y = list(Y_L_norm[l_]) + list(Y_L_norm[l_][1:-1][::-1])
      Y = (1/au2fs)*np.array(Y)
      ax.plot(X,Y, label=f"{l_str}", zorder=2) 

    # Plot total rate
    X_theta_neg = list(-1.0*np.array(X_theta))
    X = X_theta + (X_theta_neg[1:-1][::-1])
    X = np.array(X)

    Y_ = Y_tot + (Y_tot[1:-1][::-1])
    #Y_ = list(Y_tot_norm) + list(Y_tot_norm[1:-1][::-1])
    Y_ = (1/au2fs)*np.array(Y_)

    # Manually set the limits
    #ax.set_rlim(ymin, ymax)
    ax.set_rlim(ymin, ymax)
    ax.set_ylabel("Norm Adj. Rate")

    # I think Magenta is a good color for contrast with the default mpl colors
    ax.plot(X, Y_, color="m", label="Total",zorder=3)

    # Draw a bold circle at radius 0
    #print(f"rmin: {ax.get_rmin()}")
    zero_circle = plt.Circle((0, 0), -1*ax.get_rmin(), transform=ax.transData._b, color='black', linewidth=2.5, fill=False, zorder=3)
    ax.add_artist(zero_circle)

    for label in ax.get_yticklabels():
        label.set_zorder(4)  # Set a high z-order for the tick labels

    # Decrease font size of radial ticks if needed
    ax.tick_params(axis='y', which='major', labelsize=10)
    ax.grid(True)
    ax.legend()
    ax.set_title(f"Incident light angle and Angular rate dependence\n"+self.plotname)
    plt.savefig("polar_angular.png",dpi=600)
    plt.close()


    # Write CSV
    f = open("polar.csv",'w')
    f.write(f"nbasis,{self.parser.nbasis},nbsuse,{self.parser.nbsuse}\n")
    f.write("d,theta,nva99,Total Rate,S,P,D,F\n")
    for i in range(len(X_theta)):
      f.write(f"{i+1},{X_theta[i]},{nva99[i]},{fs2au*Y_tot[i]},{fs2au*Y_L[0][i]},{fs2au*Y_L[1][i]},{fs2au*Y_L[2][i]},{fs2au*Y_L[3][i]}\n")
      #f.write(f"{i+1},{X_theta[i]},{nva99[i]},{fs2au*Y_tot_norm[i]},{fs2au*Y_L_norm[0][i]},{fs2au*Y_L_norm[1][i]},{fs2au*Y_L_norm[2][i]},{fs2au*Y_L_norm[3][i]}\n")
    f.close()
        


def combine_csv_files(ndir):

  all_data = []
  max_rows = 0
  csv_files = []
  for d in range(ndir):
    csv_files.append(f"AtomLExp_rate{d}.csv")
  
  # Read all CSV files
  for file in csv_files:
    with open(file, 'r') as f:
      reader = csv.reader(f)
      data = list(reader)
      all_data.append(data)
      max_rows = max(max_rows, len(data))

  # Combine data
  combined_data = []
  for i in range(max_rows):
    row = []
    for data in all_data:
      if i < len(data):
        row.extend(data[i])
      else:
         row.extend([''] * len(data[0]))  # Add empty cells if this CSV has fewer rows
      row.append('')  # Add blank column between CSVs
    combined_data.append(row)

  # Write combined data to output file
  with open("AtomLExp_rate_combined.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(combined_data)


def __main__():
  gdvlog = sys.argv[1]
  if len(sys.argv)>2:
    rateplotter = RatePlotter(sys.argv[1], sys.argv[2])
  else:
    rateplotter = RatePlotter(sys.argv[1])

if __name__ == "__main__":  
  __main__()





