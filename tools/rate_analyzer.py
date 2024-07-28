#!/usr/bin/env python3


#
# Assume cwd is a completed tdci job, also containing a gaussian output file
# This script partitions the rate by atomic orbital, and can generate plots
#   showing contribution to the rate based on atom/quantum numbers. 
#

import numpy as np
import struct, sys, re, os, csv
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

# Parses the file 'input' for the theta/phi of each direction
def parse_tdci_input_directions():
  polar_coords = []
  f = open('input', 'r')
  l = f.readline()
  while l.strip() != "&FIELD_directions":
    l = f.readline()
    if l == "":
      print("Can't find &FIELD_directions !!!")
  l = f.readline()
  nemax = int(l.split()[-1])
  l = f.readline()
  while l.strip() != "/":
    if l == "":
      print("Can't find end of directions !!!")
      sys.exit()
    index, theta, phi = parse_tdci_input_line(l)
    polar_coords.append( (theta,phi) )
    l = f.readline()
  #print("Polar coords:")
  #print( (nemax, len(polar_coords)) )
  #print(polar_coords)
  return polar_coords
  
  
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
    self.filename = filename
    self.orbs = []
    # What should I do with 'SP' orbitals? I think just treat them as P?
    self.l_map = {'S': 0, 'P': 1, 'SP': 1, 'D': 2, 'F': 3, 'G': 4, 'H': 5, 'I': 6}
    self.inv_l_map = ['S', 'P', 'D', 'F', 'G', 'H', 'I'] # Only used for labels
    # the 'SP' entry breaks this
    #self.inv_l_map = {value: key for key, value in self.l_map.items()}
    self.atomtypes = {} # map iatom -> atom type (H, He, etc.)

  def get_params(self):
    nao = 0
    nmo = 0
    natoms = 0

    # Parse tdci output for array sizes
    f = open("OUTPUT" , 'r')
    line = f.readline()
    while line != "":
      ls = line.split()
      if len(ls) > 2:
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

    return nao, nmo, natoms

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
        self.atomtypes[atom_index] = atom_type
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

  # This should be executed after parse(), and it assumes that the list self.orbs
  #   is properly in order ( orbs[i] == orbs[i].orb_index ), as it should be
  def parse_prims(self):
    print("Entered parse_prims")
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

      if ( len(ls) == 2 ): # New prim
        for i in range(bfi_min, bfi_max+1):
          # Make sure nothing is fishy
          # Remember orbs is 0-indexed while orb_index, bfi are 1-indexed.
          assert( i == self.orbs[i-1].orb_index )
          assert( atomtype_ == self.orbs[i-1].atom_type )
          assert( atomi_ == self.orbs[i-1].atom_index )
          assert( l_ == self.orbs[i-1].l )
          # Add the prims
          #   Python requires scientific notation have an E instead of D
          self.orbs[i-1].prim_expon.append(float(ls[0].replace('D','E'))) 
          self.orbs[i-1].prim_coeff.append(float(ls[1].replace('D','E')))

      line = f.readline()

  
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
    self.nao , self.nmo, self.natoms = gdvparser.get_params()
    self.ntri = self.nao*(self.nao+1)//2
    self.CMO = None
    self.Vabs_AO = None
    self.prep_matrices()


  # Read in and transform 
  def prep_matrices(self):

    nao, nmo, ntri = self.nao, self.nmo, self.ntri
    
    vabs_tri = read_bin_array( "matrices/Vabs_AO.bin", ntri )
    vabs_MO = read_bin_array( "matrices/Vabs_MO.bin", nmo**2)

    CMO = read_bin_array( "matrices/CMO.bin", nmo*nao)
    CMO.resize((nmo,nao)) # This is correct.
    print( "CMO size:"+ str((nmo,nao)) )
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

    # Verify Vabs_MO
    mm = np.matmul
    vabs_MO.resize((nmo,nmo))
    test_vabs_AO = mm( CMO.T, mm(vabs_MO, CMO))
    print("test vabs_MO")
    print(np.shape(test_vabs_AO))
    test_vabs_AO.resize((nao**2))
    #for i in range(0,len(test_vabs_AO)):
    for i in range(0,5):
      if test_vabs_AO[i] != vabs_AO[i]:
        print(f"Vabs_AO test {i} bad: {vabs_AO[i]} vs {test_vabs_AO[i]}")

    return 0

  def make_vdens(self,timestep,e,d):
    nao, nmo = self.nao, self.nmo
    CMO = self.CMO
    dens_MO = read_bin_array( f"matrices/MO_density-e{e}-d{d}.{timestep}00.bin", nmo**2)
    dens_MO.resize((nmo,nmo))
    mm = np.matmul
    dens_AO = mm(CMO.T, mm(dens_MO, CMO))
    dens_AO.resize((nao**2))
    self.Vabs_AO.resize((nao**2))
    vdens = self.Vabs_AO * dens_AO
    vdens.resize((nao,nao))
    return vdens


class RateTrajectoryData:
  def __init__(self, nsteps, natoms, lmax, norbs):
    self.rate = np.zeros(nsteps)
    self.atomsum = np.zeros((nsteps,natoms))
    self.lsum = np.zeros((nsteps,natoms,lmax))
    self.totsum = np.zeros(nsteps) # Sanity check for rate
    self.orbrate = np.zeros((nsteps,norbs))
    self.maxrate = 0
    self.direction = 0

class RatePlotter:
  def  __init__(self, gdvlog ):
    self.parser = gdvlog_parser(gdvlog)
    parser = self.parser
    parser.parse()
    parser.parse_prims()

    self.polar_coords = parse_tdci_input_directions()
    self.dataset = [] # RateTrajectoryData indexed by direction index.
    self.ndir = len(self.polar_coords)
    print("ndir: "+str(self.ndir))

    self.nao, self.nmo, self.natoms = self.parser.get_params()
    nao, nmo, natoms = self.nao, self.nmo, self.natoms
    self.nsteps = 160
    self.lmax = 6
    nsteps = self.nsteps
    lmax = self.lmax
   
    e_idx = 1 # Field number
    d_idx = 1 # Direction

    self.rc = ratechecker(self.parser)

    for d in range(0,self.ndir):
      d_data = RateTrajectoryData(self.nsteps, self.natoms, self.lmax, len(parser.orbs))
      rate = d_data.rate
      atomsum = d_data.atomsum
      lsum = d_data.lsum
      totsum = d_data.totsum
      orbrate = d_data.orbrate
      d_data.direction = d
      expdict = {}


      # Generate index keys for angular and exponent breakdown
      for orb in parser.orbs:
        if len(orb.prim_expon)>0:
          l_ = orb.l
          exp_ = orb.prim_expon[0]
          expdict[ (orb.atom_index, l_, exp_) ] = np.zeros(self.nsteps)

      # Generate density*vabs in AO basis and partition by atom and L
      for time in range(0, nsteps):
        vdens = self.rc.make_vdens(time+1,e_idx,d+1)
        rate[time] = np.trace(vdens)
        if time == nsteps-1:
          print(f"Final rate for d={d+1} : {rate[time]} | {rate[time]/au2fs:.3f}")
        # Atoms/orbitals are 1-indexed!
        #print("Breakdown of rate by atom and angular momentum")

        for orb in parser.orbs:
          idx = orb.orb_index-1
          orbrate[time][idx] = vdens[idx][idx]
          if len(orb.prim_expon)>0:
            expdict[ (orb.atom_index, orb.l, orb.prim_expon[0]) ][time] += vdens[idx][idx]
          

        for iatom in range(0,self.natoms):
          for l_ in range(0, lmax):
            for orb in parser.orbs:
              if orb.l == l_ and orb.atom_index == iatom+1:
                #print(orb.orb_index)
                #print(len(a.orbs))
                idx = orb.orb_index-1 # orb_index is 1-indexed
                lsum[time][iatom][l_] += vdens[idx][idx]
                atomsum[time][iatom] += vdens[idx][idx]
          totsum[time] += atomsum[time][iatom]
        if np.abs(totsum[time]-rate[time])>1e-3:
          print(f"totsum doesnt match rate: totsum={totsum[time]} , rate={rate[time]}")
      #self.dataset.append(d_data)
      #self.rate_by_time_plots(d_data)
      #self.individual_orb_plot(d_data, parser)
      #self.AtomLExp_plot(d_data, parser, expdict)
      self.AtomLExp_AvgCSV(d_data, parser, expdict)
    #self.polar_angular_plot()
    #combine_csv_files(self.ndir)
    self.AtomLExp_DirAvgCSV()
    

  def AtomLExp_plot(self,d_data, parser, expdict):
    orbs = parser.orbs
    nsteps = self.nsteps
    norbs = len(orbs)
    #maxrates = np.max(d_data.orbrate, axis=0)
    maxrates = np.zeros( len(expdict) )
    keys = []
    data = []
    for key, val in expdict.items():
      keys.append(key)
      data.append(val)
    data = np.array(data)
    maxrates = np.max(data, axis=1)
    #print("maxrates:")
    #print(maxrates)
    orb_idx_ratesorted = np.argsort(maxrates)[::-1]
    MAX_LINES = 10 # max number of lines to plot

    X = np.array(list(range(1,nsteps+1)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)


    for i in range(0,MAX_LINES):
      idx = orb_idx_ratesorted[i]
      key = keys[idx]
      l_str = parser.inv_l_map[key[1]]
      atype = parser.atomtypes[key[0] ] 
      ax.plot( X, data[idx],
        label=f"{atype}{key[0]}, {l_str}, {key[2]}, {maxrates[idx]:0.3f}")
    ax.legend(title=f"Atom, L, Exp, MaxRate")
    ax.title.set_text(f"Top 10 orbitals, Dir: {d_data.direction}")
    
    plt.tight_layout()
    plt.savefig(f"AtomLExp{d_data.direction}.png", dpi=200)
    plt.close()


    # Write CSV
    f = open(f"AtomLExp_rate{d_data.direction}.csv", 'w')
    f.write(f"Direction,{d_data.direction},,,\n")
    f.write("AtomIdx,AtomType,L,Exponent,MaxRate\n")
    for i in range(0,len(maxrates)):
      idx = orb_idx_ratesorted[i]
      # key is ( atomidx, l, exponent )
      key = keys[idx]
      atype = parser.atomtypes[ key[0] ]
      l_str = parser.inv_l_map[key[1]]
      
      outstr = f"{key[0]},{atype},{l_str},{key[2]},{maxrates[idx]}\n"
      #print(outstr)
      f.write(outstr)
    f.close()


  # CSV Summary Average Table
  def AtomLExp_AvgCSV(self, d_data, parser, expdict): 
    f = open(f"AtomLExp_AvgRate{d_data.direction}.csv", 'w')
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
      header = header + f",{atype_list[i]}"
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
          outstr = outstr + f",{entries[i]}"
        if nowrite: # Skip basis elements that do not exist
          pass
        else:
          f.write(outstr)
    f.close()
        
  # Combines the CSV files 
  def AtomLExp_DirAvgCSV(self):
    csv_files = [f"AtomLExp_AvgRate{d}.csv" for d in range(0,self.ndir)]
    
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
    nsteps = self.nsteps
    norbs = len(orbs)
    maxrates = np.max(d_data.orbrate, axis=0)
    orb_idx_ratesorted = np.argsort(maxrates)[::-1]
    MAX_LINES = 10 # max number of lines to plot

    X = np.array(list(range(1,nsteps+1)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(0,MAX_LINES):
      idx = orb_idx_ratesorted[i]
      orb = orbs[idx]
      l_str = parser.inv_l_map[orb.l]
      ax.plot( X, d_data.orbrate[:,idx],
        label=f"{orb.atom_type}, {l_str}, {orb.prim_expon}, {maxrates[idx]:0.3f}")
    ax.legend(title=f"Atom, L, Exp, MaxRate")
    ax.title.set_text(f"Top 10 orbitals, Dir: {d_data.direction}")
    
    plt.tight_layout()
    plt.savefig(f"orbsort{d_data.direction}.png", dpi=200)
    plt.close()
      
      
      

  def rate_by_time_plots(self,d_data):
    nsteps = self.nsteps
    lmax = self.lmax
    parser = self.parser
    natoms = self.natoms
    maxrate = 1e-5 # forgot where to get this from
    for iatom in range(0,natoms):
      for l_ in range(0,lmax):
        maxrate = max( max(d_data.lsum[:,iatom,l_]), maxrate)
    # Plot by Atom and Angular momentum
    X = np.array(list(range(1,nsteps+1)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for iatom in range(0,natoms):
      for l_ in range(0, lmax):
        if max(d_data.lsum[:,iatom,l_]) < 0.10*maxrate:
          continue # Skip low contributions to avoid clutter
        l_str = parser.inv_l_map[l_]
        #print(parser.atomtypes, flush=True)
        atomstr = str(iatom+1)+parser.atomtypes[iatom+1]
        ax.plot( X , d_data.lsum[:, iatom, l_], label=f"{atomstr}: {l_str}" )

    ax.legend()
    ax.title.set_text(f"Rate by Atom:Angular Momentum. Dir: {d_data.direction}")
    
    plt.tight_layout()
    plt.savefig(f"atomangular{d_data.direction}.png", dpi=200)
    plt.close()

    # Plot by angular momentum, no atoms
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    l_noatom = np.sum(d_data.lsum, axis=1) # Now has shape l_noatom[nsteps][lmax]
    for l_ in range(0, lmax):
      #if max(l_noatom[:,l_]) < 0.1*maxrate:
      #  continue # Skip low contributions to avoid clutter
      l_str = parser.inv_l_map[l_]
      ax.plot( X , l_noatom[:, l_], label=f"{l_str}" )

    ax.legend()
    ax.title.set_text(f"Rate by Angular Momentum. Dir:{d_data.direction}")
    
    plt.tight_layout()
    plt.savefig(f"angular{d_data.direction}.png", dpi=200)
    plt.close()


    # Plot by atom, no angular momentum
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for iatom in range(0,natoms):
      #if max(atomsum[:,iatom]) < 0.1*maxrate:
      #  continue # Skip low contributions to avoid clutter
      #print(parser.atomtypes, flush=True)
      atomstr = str(iatom+1)+parser.atomtypes[iatom+1]
      ax.plot( X , d_data.atomsum[:, iatom], label=f"{atomstr}" )

    ax.legend()
    ax.title.set_text(f"Rate by Atom. Dir:{d_data.direction}")
    
    plt.tight_layout()
    plt.savefig(f"atom{d_data.direction}.png", dpi=200)
    plt.close()

    return 0 

  def polar_angular_plot(self):
    # For now we just ignore phi i think
    print( "Length check (dataset, polar_coords)")
    print( (len(self.dataset), len(self.polar_coords)))
    ndir = len(self.polar_coords)
    X_theta = []
    for x in self.polar_coords:
      X_theta.append( x[0] * (np.pi/180.)  )
    Y_L = [ [] for _ in range(self.lmax) ]
    time = (self.nsteps-1)
    for d in range(0,ndir):
      for l_ in range(0,self.lmax):
        rate_l = 0.0
        for iatom in range(0,self.natoms):
          rate_l += self.dataset[d].lsum[time][iatom][l_]
        Y_L[l_].append(rate_l)

    plt.clf()
    fig, ax = plt.subplots( subplot_kw={'projection':'polar'})
    for l_ in range(0,self.lmax):
      l_str = self.parser.inv_l_map[l_]
      print(X_theta)
      print(Y_L[l_])
      #ax.plot(X_theta, Y_L[l_], label=f"{l_str}") 
      # 180 degree symmetry
      X_theta_neg = list(-1.0*np.array(X_theta))
      X = X_theta + (X_theta_neg[1:-1][::-1])
      Y = Y_L[l_] + (Y_L[l_][1:-1][::-1])
      ax.plot(X,Y, label=f"{l_str}") 
    ax.grid(True)
    ax.legend()
    ax.set_title("Incident light angle and Angular rate dependence")
    plt.savefig("polar_angular.png",dpi=600)
    plt.close()
        


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
  rateplotter = RatePlotter(sys.argv[1])

  
__main__()





