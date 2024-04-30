#!/usr/bin/env python3


#
# Assume cwd is a completed tdci job, also containing a gaussian output file
# This script partitions the rate by atomic orbital, and can generate plots
#   showing contribution to the rate based on atom/quantum numbers. 
#

import numpy as np
import struct, sys
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

# globals
au2fs = 2.418884326509e-2

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

class orb:
  def __init__(self):
    orb_index = 0
    atom_index = 0
    atom_type = 0
    n = 0
    l = 0
    m = 0

class gdvlog_parser:
  def __init__(self, filename):
    self.filename = filename
    self.orbs = []

  def get_params(self):
    nao = 128
    nmo = 22
    natoms = 3
    return nao, nmo, natoms

  # Return tuple (n, l, m)
  # Examples:
  # "1S" -> (1, 0, 0)
  # "9PX" -> (9, 1, -1)
  # "20D+2" -> (20, 2, 2)
  # "29F 0" -> (29, 3, 0) 
  def decipher_orbstring(self,orb_string):
    l_map = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}

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
  def parse(self):
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
      neworb = orb()

      orb_index = int(line[0])

      if len(line) == 6: # New atom
        atom_index = int(line[1])
        atom_type = line[2]
        orbstring = line[3] # Never will have space like '16D 0'
      elif len(line) == 4: # 154       20D+2        0.12572   1.79229
        orbstring = line[1]
      elif len(line) == 5: # 131       16D 0        0.00000   0.00000
        orbstring = line[1]+line[2]
        
      
      # Save info
      neworb.orb_index = orb_index
      neworb.atom_index = atom_index
      neworb.atom_type = atom_type
      neworb.n, neworb.l, neworb.m = self.decipher_orbstring(orbstring)
      self.orbs.append(neworb) # Append to output data
      # End loop, make sure next line doesn't crash the while loop
      line = f.readline().split()
      if len(line)==0:
        break
      elif len(line[0])==0:
        break
    print("Done Parsing!")

##############################
# End gdvlog_parser class
##############################

class ratechecker:
  def __init__(self, gdvparser):
    self.gdvparser = gdvparser
    self.nao , self.nmo, self.natoms = gdvparser.get_params()
    self.ntri = self.nao*(self.nao+1)//2
    self.CMO = None
    self.Vabs_AO = None
    self.prep_matrices()


  # Read in and transform 
  def prep_matrices(self):

    nao, nmo, ntri = self.nao, self.nmo, self.ntri
    
    vabs_tri = read_bin_array( "matrices/Vabs_AO.bin", ntri )
    #vabs_MO = read_bin_array( "matrices/Vabs_MO.bin", nmo**2)

    CMO = read_bin_array( "matrices/CMO.bin", nmo*nao)
    CMO.resize((nmo,nao)) # This is correct.
    self.CMO = CMO

    # Unpack triangular matrix
    vabs_AO = np.zeros((nao,nao))

    for iao in range(0,nao):
      for jao in range(0,nao):
        ij = (iao+1)*(iao)//2 + (jao+1)
        if (jao>iao) :
          ij = (jao+1)*(jao)//2 + (iao+1)
        vabs_AO[iao][jao] = vabs_tri[ij-1]

    vabs_AO.resize((nao**2))
    self.Vabs_AO = vabs_AO
    return 0

  def make_vdens(self,timestep):
    nao, nmo = self.nao, self.nmo
    CMO = self.CMO
    dens_MO = read_bin_array( "matrices/MO_density-e1-d1."+str(timestep)+"00.bin", nmo**2)
    dens_MO.resize((nmo,nmo))
    mm = np.matmul
    dens_AO = mm(CMO.T, mm(dens_MO, CMO))
    dens_AO.resize((nao**2))
    self.Vabs_AO.resize((nao**2))
    vdens = self.Vabs_AO * dens_AO
    vdens.resize((nao,nao))
    return vdens



def __main__():

  parser = gdvlog_parser("h2o_gaussian.log")

  parser.parse()

  nao = 128
  nmo = 22
  natoms = 3
  nsteps = 160
  lmax = 4

  rc = ratechecker(parser)

  rate = np.zeros(nsteps)
  atomsum = np.zeros((nsteps,natoms))
  lsum = np.zeros((nsteps,natoms,lmax))
  totsum = np.zeros(nsteps) # Sanity check for rate

  # Generate density*vabs in AO basis and partition by atom and L
  for time in range(0, nsteps):
    vdens = rc.make_vdens(time+1)
    rate[time] = np.trace(vdens)
    # Atoms/orbitals are 1-indexed!
    #print("Breakdown of rate by atom and angular momentum")

    for iatom in range(0,natoms):
      for l_ in range(0, lmax):
        for orb in parser.orbs:
          if orb.l == l_ and orb.atom_index == iatom+1:
            #print(orb.orb_index)
            #print(len(a.orbs))
            idx = orb.orb_index-1 # orb_index is 1-indexed
            lsum[time][iatom][l_] += vdens[idx][idx]
            atomsum[time][iatom] += vdens[idx][idx]
      totsum[time] += atomsum[time][iatom]
    
    #for iatom in range(1,natoms+1):
    #  for l_ in range(0, lmax):
    #    print(f"atom {iatom+1}, l={l_}: {lsum[iatom][l_]: 4.2e} ({lsum[iatom][l_]/rate: .2f})")
    #  print(f"atom {iatom+1} : {atomsum[iatom]: 4.2e} ({atomsum[iatom]/rate: .2f})")
    #print(f"Total rate : {totsum: 4.2e} ({totsum/rate: .2f} of Tr(vdens))")

  # Plotting

  # Plot by Atom and Angular momentum
  X = np.array(list(range(1,nsteps+1)))
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  for iatom in range(0,natoms):
    for l_ in range(0, lmax):
      ax.plot( X , lsum[:, iatom, l_], label=f"{iatom+1}: {l_}" )

  ax.legend()
  ax.title.set_text("Rate by Atom:Angular Momentum")
  
  plt.tight_layout()
  plt.savefig("atomangular.png", dpi=200)

  # Plot by angular momentum, no atoms
  plt.clf()
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  l_noatom = np.sum(lsum, axis=1) # Now has shape l_noatom[nsteps][lmax]
  for l_ in range(0, lmax):
    ax.plot( X , l_noatom[:, l_], label=f"{l_}" )

  ax.legend()
  ax.title.set_text("Rate by Angular Momentum")
  
  plt.tight_layout()
  plt.savefig("angular.png", dpi=200)


  # Plot by atom, no angular momentum
  plt.clf()
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  for iatom in range(0,natoms):
    ax.plot( X , atomsum[:, iatom], label=f"{iatom+1}" )

  ax.legend()
  ax.title.set_text("Rate by Atom")
  
  plt.tight_layout()
  plt.savefig("atom.png", dpi=200)

  return 0 
  
__main__()





