
import os, sys, subprocess, re, struct
import numpy as np
import matplotlib.pyplot as plt

# Tools for automating density movies


# Config variables for grid bounds (edit these)
# Angstroms
xmin = ( -12.0, -12.0, -12.0 ) # .cube region minimum corner
xmax = ( 12.0, 12.0, 12.0 )    # .cube region maximum corner
nsteps = 45                    # .cube grid points in each direction (used in both 3D and 1D case)

oneD_maxdist = 16 # In 1D case, this is the length of the sampling line in Angstroms

# TDCI timesteps to create .cube files for
start_timestep_default = 1
end_timestep_default = 160


# Path to folder density2fchk executable
#tools = "/home/andy/tdci-test/TDCI-CAP/tools/" 
tools = os.path.dirname(os.path.abspath(__file__))
print(tools)



from optparse import OptionParser
usage = "usage: python3 movietool.py [options] [gaussian_directory] [tdci_directory] [field direction index]"
parser = OptionParser(usage=usage)
parser.add_option("-v", action="store_true", default=False, dest="verbose", help="Print extra info for debugging")
parser.add_option("--cube_only", action="store_true", default=False, dest="cube_only", help="Create .cube files from density in TDCI directory")
parser.add_option("--render_only", action="store_true", default=False, dest="render_only", help="Run VMD and FFMPEG")
parser.add_option("--vmd_only", action="store_true", default=False, dest="vmd_only", help="Run only VMD")
parser.add_option("--ffmpeg_only", action="store_true", default=False, dest="ffmpeg_only", help="Run only ffmpeg")
parser.add_option("--start_timestep", action="store", type="int", default=start_timestep_default, dest="start_timestep", help=f"Set the timestep to start on. Default {start_timestep_default}")
parser.add_option("--end_timestep", action="store", type="int", default=end_timestep_default, dest="end_timestep", help=f"Set the timestep to end on. Default {end_timestep_default}")
parser.add_option("--isoval1", action="store", type="float", default="0.002", dest="isoval1", help=f"Inner isovalue in VMD render. Default 0.002")
parser.add_option("--isoval2", action="store", type="float", default="0.00095", dest="isoval2", help=f"Outer isovalue in VMD render. Default 0.00095")
parser.add_option("--diff", action="store", type="int", default="0", dest="diff", help=f"0 - Plot Density. 1 - Plot Density minus the HF density. 2 - Plot Density minus the density in the template fchk file. 3 - Plot Density minus the density at timestep 1. NOTE: OPTIONS 2 AND 3 ARE NOT CURRENTLY SUPPORTED")
parser.add_option("--eidx", action="store", type="int", default="1", dest="eidx", help=f"Index of Emax. Ex: --eidx=2 will read MO_density-e2-d[dir].[timestep].bin . Default 1.")
parser.add_option("--1D", action="store_true", default=False, dest="oneD", help="Generate a 1D cubefile for simple plotting. Density is sampled along a 16 Angstrom line from the origin. See theta and phi options.")
parser.add_option("--theta", action="store", type="float", default=None, dest="theta", help=f"In degrees. For 1D cubefile, selects the theta direction of the density slice. Default aligns with the field direction.")
parser.add_option("--phi", action="store", type="float", default=None, dest="phi", help=f"In degrees. For 1D cubefile, selects the phi direction of the density slice. Default aligns with the field direction.")
parser.add_option("--orblist", action="store", type="string", default="", dest="orblist", help="Comma separated indices, no spaces, e.g. 0,3,5 . For 1D plots only. Only plot the density of orbitals in the orblist.")
(opts, args) = parser.parse_args()



# Determine number of CPUs we can use for cube file generation.
p = subprocess.Popen("nproc", shell=True, stdout=subprocess.PIPE)
p.wait()
NCPU = int(p.stdout.read().decode().strip())
del p # remove this symbol from global namespace
if NCPU > 3:
  NCPU = NCPU - 1 # Let's leave one free ^_^



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



def sph2cart(theta, phi, units="Radians"):
  if units.lower()[0] == "d":
    theta *= (180.0/np.pi)
    phi *= (180.0/np.pi)
  return (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))



def load_cube_data(filename):
  # First, determine the number of header lines to skip.
  with open(filename, 'r') as f:
    # Skip two comment lines.
    f.readline()
    f.readline()
    header_line = f.readline().split()
    natoms = abs(int(header_line[0]))
    # Total header lines = 2 (comments) + 1 (atom count + origin)
    #                     + 3 (grid vectors) + natoms (atoms)
    header_lines = 2 + 1 + 3 + natoms

  # Load the data into a NumPy array.
  data = np.concatenate(np.loadtxt(filename, skiprows=header_lines))
  return data

def plot_density_z(density, outfile):
  dz = (xmax[2]-xmin[2])/nsteps  # globals
  z = np.arange(len(density)) * dz
  plt.figure(figsize=(6, 4))
  plt.plot(z, density, marker='o', linestyle='-')
  plt.xlabel('z (Angstrom)')
  plt.ylabel('Density')
  plt.title('Density Along the z-Axis.')
  plt.tight_layout()
  plt.savefig(outfile)
  plt.close()

def plot_density_z_series(start_timestep, end_timestep):
  for i in range(start_timestep, end_timestep+1):
    density = load_cube_data(f"out{i}.cube")
    plot_density_z(density, f"densz{i}.png")

def plot_density_1D_allinone(start_timestep, end_timestep,thetaphi=(0,0),cubeprefix="out"):
  dz = (xmax[2]-xmin[2])/nsteps  # globals
  plt.figure(figsize=(6, 4))
  plt.xlabel('z (Angstrom)')
  plt.ylabel('Density')
  plt.title(f'Density Along sample line theta={thetaphi[0]}, phi={thetaphi[1]}. Blue=Early time, Red=Late time')

  for i in range(start_timestep, end_timestep+1):
    density = load_cube_data(f"{cubeprefix}{i}.cube")
    z = np.arange(len(density)) * dz
    color = interpolate_color("#0000FF", "#FF0000", (end_timestep-start_timestep), i-start_timestep)
    plt.plot(z, density, marker='o', linestyle='-', color=color, alpha=0.6)

  plt.legend()
  plt.tight_layout()
  plt.savefig(f"{cubeprefix}_dens.png", dpi=300)
  plt.close()


def interpolate_color(start_hex: str, end_hex: str, Nstep: int, step: int) -> str:
  start_hex = start_hex.lstrip('#')
  end_hex = end_hex.lstrip('#')
  # Convert hex to integer values for each color channel (R, G, B)
  start_rgb = [int(start_hex[i:i+2], 16) for i in (0, 2, 4)]
  end_rgb = [int(end_hex[i:i+2], 16) for i in (0, 2, 4)]
  t = step / Nstep
  interpolated_rgb = [round(start + (end - start) * t) for start, end in zip(start_rgb, end_rgb)]
  return '#' + ''.join(f'{val:02x}' for val in interpolated_rgb)



# assuming there is only one in the directory, return path of chk file or fchk file
def find_chk_file(directory):
  for filename in os.listdir(directory):
    if filename.endswith('.fchk'):
      return os.path.join(directory, filename)
  for filename in os.listdir(directory):
    if filename.endswith('.chk'):
      return os.path.join(directory, filename)
  return None


# Get the molecular coordinates from a .fchk file
def read_fchk_coords(filename):
  coords = []
  atomic_nums = []
  coord_vals = []
  reading_atomic = False
  reading_coords = False
  
  with open(filename, 'r') as f:
    for line in f:
      if 'Atomic numbers' in line:
        reading_atomic = True
        continue
      
      if reading_atomic:
        if line.strip() and not line.strip()[0].isalpha():
          atomic_nums.extend(map(int, line.split()))
        else:
          reading_atomic = False
      
      if 'Current cartesian coordinates' in line:
        reading_coords = True
        continue
      
      if reading_coords:
        if line.strip() and not line.strip()[0].isalpha():
          coord_vals.extend(map(float, line.split()))
        else:
          break

  for i in range(len(atomic_nums)):
    coords.append((atomic_nums[i], coord_vals[i*3], coord_vals[i*3+1], coord_vals[i*3+2]))
  
  return coords

class cubemaker:
  def __init__(self, gaussian_dir, tdci_dir):
    self.gaussian_dir = gaussian_dir
    self.tdci_dir = tdci_dir
    self.Z_map = {'H': 1, 'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Cl':17, 'Br':35, 'Ar':18, 'I':53 } 
    self.inv_Z_map = {value: key for key, value in self.Z_map.items()}
    # i assume this charge map is basis dependent, but it probably doesnt matter for
    #  plotting densities
    self.Chg_map = {'H':1.0, 'C':6.0, 'O':8.0, 'Ar':18.0, 'I':25.0 }
    self.Chg_Zmap = {1:1.0, 6:6.0, 8:8.0, 18:18.0, 35:35.0, 53:25.0 }

  # https://gaussian.com/cubegen/
  # cubegen needs a template cube file to specify the boundaries
  #   and step-size of the grid. only header is needed.

  # x_min, x_max : 3-tuples of the bounds of the cube
  # nstep : number of steps
  # coords: list of 4-tuples, tup[0] : Z, tup[1-3] : x,y,z coords in angstroms
  def generate_template_cube(self, x_min, x_max, nstep, coords ):
    natoms = len(coords)
    f = open("input.cube", 'w')
    f.write(" HCCI opt FDENSITY=SCF\n Electron density from Total SCF Density\n")
    # NAtoms, X-Origin, Y-Origin, Z-Origin, NVal  (NVal is the # of values/point)
    f.write(f"{natoms:>4d}   {x_min[0]:9.4f}   {x_min[1]:9.4f}   {x_min[2]:9.4f}    1\n")
    f.write(f"{nstep:>4d}   {(x_max[0]-x_min[0])/nstep:9.4f}   {0.0:9.4f}   {0.0:9.4f}\n")
    f.write(f"{nstep:>4d}   {0.0:9.4f}   {(x_max[1]-x_min[1])/nstep:9.4f}   {0.0:9.4f}\n")
    f.write(f"{nstep:>4d}   {0.0:9.4f}   {0.0:9.4f}   {(x_max[2]-x_min[2])/nstep:9.4f}\n")
    for atom in coords:
      #z = self.inv_Z_map[atom[0]]
      z = atom[0]
      chg = self.Chg_Zmap[atom[0]]
      #print( (z,chg,coords) )
      f.write(f"{z:>4d}   {chg:9.4f}   {atom[1]:9.4f}   {atom[2]:9.4f}   {atom[3]:9.4f}\n")
    # Normally the grid values follow, but not needed for this template
    f.close()
    return 0

  # For making a 1-d strip in the z-direction
  def generate_template_cube_1D_polar(self, coords, theta=0.0, phi=0.0, maxdist=oneD_maxdist, nstep=nsteps):
    natoms = len(coords)
    nx, ny, nz = sph2cart(theta,phi)
    dx, dy, dz = nx*maxdist/nstep, ny*maxdist/nstep, nz*maxdist/nstep
    print(f"1DCube: {nstep} steps in direction theta={theta}, phi={phi}")
    print(f"    dx,dy,dx=({dx},{dy},{dz}) ")
    f = open("input.cube", 'w')
    f.write(" HCCI opt FDENSITY=SCF\n Electron density from Total SCF Density\n")
    # NAtoms, X-Origin, Y-Origin, Z-Origin, NVal  (NVal is the # of values/point)
    f.write(f"{natoms:>4d}   {0.0:9.4f}   {0.0:9.4f}   {0.0:9.4f}    1\n")
    f.write(f"{1:>4d}   {0.0:9.4f}   {0.0:9.4f}   {0.0:9.4f}\n")
    f.write(f"{1:>4d}   {0.0:9.4f}   {0.0:9.4f}   {0.0:9.4f}\n")
    f.write(f"{nstep:>4d}   {dx:9.4f}   {dy:9.4f}   {dz:9.4f}\n")
    for atom in coords:
      #z = self.inv_Z_map[atom[0]]
      z = atom[0]
      chg = self.Chg_Zmap[atom[0]]
      #print( (z,chg,coords) )
      f.write(f"{z:>4d}   {chg:9.4f}   {atom[1]:9.4f}   {atom[2]:9.4f}   {atom[3]:9.4f}\n")
    # Normally the grid values follow, but not needed for this template
    f.close()
    return 0


# Call the density2fchk fortran program to insert the density from tdci into the template fchk file
# 'direction' and 'timestep' are ignored if 'filename' is present.
def density2fchk(tdci_dir, nocc, norb, direction=1, timestep=1, diff=0, filename=None):
  #!: density2fchk.f90 Usage:
  #!:   density2fchk infile MOdensityfile timestep diff iscale
  #!:
  #!: Command line input:
  #!:    infile -  full name of the template formatted checkpoint file (with ".fchk")
  #!:    MOdensityfile - full name of the ion coefficient file (e.g. MO_density-e1-d1.dat)
  #!:    timestep - step number to be used for the constructiion of the NTOs and hole density
  #!:    diff = 0 - calculate the density
  #!:         = 1 - calculate the density minus the density of the field-free neutral
  #!:         = 2 - calculate the density minus the density in the template fchk file
  #!:         = 3 - calculate the density minus the density at time step 1
  #!:    iscale = 0 - do not scale NTOs or hole density
  #!:           = 1 - scale by the ratio of the rate to the maximum rate
  #!:
  #!: Output
  #!:    infile-1.fchk - formatted checkpoint file with total density
  #!:    total SCF density replaced bt the density or density difference
  # See density2fchk.f90 for more details

  # iscale=0 : dont scale
  # iscale=1 : read from maxrate.dat and scale by 1./maxrate (disabled)
  iscale = 0

  if filename is None:
    print(timestep)
    densfile = f"{tdci_dir}/matrices/MO_density-e1-d{direction}.{timestep}00.bin"
  else:
    densfile = filename

  if diff==3:
    densfile_t0 = f"{tdci_dir}/matrices/MO_density-e1-d{direction}.100.bin"
  else:
    densfile_t0 = ""

  # density2fchk template_fchk density_binary, timestep, diff, iscale, nocc, norb
  cmd = f"{tools}/density2fchk temp.fchk {densfile} {timestep} {diff} {iscale} {nocc} {norb} {densfile_t0}"
  print(f"Running: {cmd}") 
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
  p.wait()
  out, err = p.communicate()
  if verbose:
    print(out)
    print(err)
  return 0

# Run cubegen using previously generated input.cube as a template for the region to sample 
def cubegen(outname):
  cpus = NCPU
  ngrid = "-1"
  print(f"Generating cube : {outname}.cube")
  p = subprocess.Popen(f"cubegen {cpus} FDENSITY=SCF temp-1.fchk {outname}.cube {ngrid} h input.cube >> output", shell=True)
  p.wait()


# Call VMD to visualize the cubefile and render to a bitmap bmp file
# You will need to modify camera settings in the template tcl file
# We plot two nested isovalues. They are the same color.
def vmd_cube2bmp(start_timestep, end_timestep, isoval1=0.00095, isoval2=0.002):
  start_timestepp1 = int(start_timestep)+1
  sh_cmds = []
  sh_cmds.append("rm -rf *.bmp") # Clean old output
  sh_cmds.append(f"cp {tools}/VMD_template.tcl VMD_run.tcl") # Copy template
  # Run replacements on template
  sh_cmds.append(f"sed -i 's%TSTARTKEY%{start_timestep}%g' VMD_run.tcl")
  sh_cmds.append(f"sed -i 's%TSTARTP1KEY%{start_timestepp1}%g' VMD_run.tcl")
  sh_cmds.append(f"sed -i 's%TENDKEY%{end_timestep}%g' VMD_run.tcl")
  sh_cmds.append(f"sed -i 's%ISOKEY1%{isoval1:0.8f}%g' VMD_run.tcl")
  sh_cmds.append(f"sed -i 's%ISOKEY2%{isoval2:0.8f}%g' VMD_run.tcl")
  sh_cmds.append(f"vmd -dispdev text -e VMD_run.tcl")
  p = subprocess.Popen(";".join(sh_cmds),shell=True)
  p.wait()

# See if any vmd processes are still running
def check_vmd_orphans():
  p = subprocess.Popen("ps -ef | grep 'vmd' | grep -v grep", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
  out, err = p.communicate()
  if out.strip():
    for line in out.strip().split("\n"):
      if not line.strip(): continue
      print(f"Orphan vmd process: {line.split()[1]}")

# https://ffmpeg.org/ffmpeg.html
def ffmpeg_bmp2mp4(start_timestep):
  sh_cmds = []
  # if the input -framerate is lower than the output -r then ffmpeg will duplicate frames to reach your desired output frame rate
  sh_cmds.append(f'ffmpeg -y -framerate 5 -start_number {start_timestep} -i "render%d.bmp" -c:v libx264 -r 30 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" animation.mp4')
  p = subprocess.Popen(";".join(sh_cmds),shell=True)
  p.wait()


# Parses the file 'input' for the theta/phi of each direction
def parse_tdci_input(tdci_dir):
  polar_coords = []
  f = open(f'{tdci_dir}/input', 'r')
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




# Grab orbital space sizes from tdci OUTPUT file
def parse_tdci_output(tdcidir):
  f = open(tdcidir+"OUTPUT", 'r')
  nocc = 0
  nva = 0
  norb = 0  
  for line in f:
    ls = line.split()
    if len(ls) == 3:
      if ls[0] == "noa":
        nocc = int(ls[2])
      if ls[0] == "nva":
        nva = int(ls[2])
  norb = nocc+nva
  return nocc, norb

# set up necessary templates and generate cube 
#   at each timestep from tdci density bin files
def generate_cube_series(gaussian_dir, tdci_dir, direction, start_timestep, end_timestep, diff=0):
  nocc, norb = parse_tdci_output(tdci_dir)
  # make fchk and get coords from it
  chkpath = find_chk_file(gaussian_dir)  
  if chkpath.endswith('.fchk'):
    print(f"Found fchk! Copying: {chkpath} to be used.")
    p = subprocess.Popen(f"cp {chkpath} temp.fchk", shell=True)
    p.wait()
  else:
    print(f"Found chk! Copying: {chkpath} to be formatted and then used.")
    p = subprocess.Popen(f"cp {chkpath} temp.chk", shell=True)
    p.wait()
    p = subprocess.Popen(f"formchk temp.chk", shell=True)
    p.wait()

  # create input.cube template
  coords = read_fchk_coords("temp.fchk")
  cm = cubemaker( gaussian_dir, tdci_dir )
  cm.generate_template_cube( xmin, xmax, nsteps, coords )
  
  for i in range(start_timestep, end_timestep+1):
    density2fchk(tdci_dir, nocc, norb, direction=direction, timestep=i,  diff=diff)
    cubegen(f"out{i}")

# set up necessary templates and generate cube
#   at each timestep from tdci density bin files
def generate_cube_series_1D(gaussian_dir, tdci_dir, direction, start_timestep, end_timestep, thetaphi, diff=0):
  nocc, norb = parse_tdci_output(tdci_dir)
  # make fchk and get coords from it
  chkpath = find_chk_file(gaussian_dir)
  if chkpath.endswith('.fchk'):
    print(f"Found fchk! Copying: {chkpath} to be used.")
    p = subprocess.Popen(f"cp {chkpath} temp.fchk", shell=True)
    p.wait()
  else:
    print(f"Found chk! Copying: {chkpath} to be formatted and then used.")
    p = subprocess.Popen(f"cp {chkpath} temp.chk", shell=True)
    p.wait()
    p = subprocess.Popen(f"formchk temp.chk", shell=True)
    p.wait()

  # create input.cube template
  coords = read_fchk_coords("temp.fchk")
  cm = cubemaker( gaussian_dir, tdci_dir )
  # Global oneD_maxdist and nsteps
  cm.generate_template_cube_1D( coords, thetaphi[0], thetaphi[1], maxdist=oneD_maxdist, nsteps=nsteps )

  for i in range(start_timestep, end_timestep+1):
    density2fchk(tdci_dir, direction, i, nocc, norb)
    cubegen(f"out{i}")

# orbsubset should be a list of 0-indexed orbital indices to include in the density matrix
# everything else will be zero'd out to give a density plot of the selected orbitals
def generate_cube_series_1D_orbsubset(gaussian_dir, tdci_dir, direction, start_timestep, end_timestep, thetaphi, orbsubset=[], diff=0):
  nocc, norb = parse_tdci_output(tdci_dir)
  # make fchk and get coords from it
  chkpath = find_chk_file(gaussian_dir)
  if chkpath.endswith('.fchk'):
    print(f"Found fchk! Copying: {chkpath} to be used.")
    p = subprocess.Popen(f"cp {chkpath} temp.fchk", shell=True)
    p.wait()
  else:
    print(f"Found chk! Copying: {chkpath} to be formatted and then used.")
    p = subprocess.Popen(f"cp {chkpath} temp.chk", shell=True)
    p.wait()
    p = subprocess.Popen(f"formchk temp.chk", shell=True)
    p.wait()

  # create input.cube template
  coords = read_fchk_coords("temp.fchk")
  cm = cubemaker( gaussian_dir, tdci_dir )
  # Global oneD_maxdist and nsteps
  cm.generate_template_cube_1D_polar( coords, thetaphi[0], thetaphi[1], maxdist=oneD_maxdist, nstep=nsteps )

  for i in range(start_timestep, end_timestep+1):
    # Read density file, copy selected orbitals
    densfile = f"{tdci_dir}/matrices/MO_density-e1-d{direction}.{i}00.bin"
    dens = read_bin_array(densfile, norb**2)
    dens.resize((norb,norb))
    if len(orbsubset) > 0:
      tempdens = np.zeros((norb,norb))
      for orb in orbsubset:
        tempdens[orb, orbsubset] = dens[orb, orbsubset]
        tempdens[orbsubset, orb] = dens[orbsubset, orb]
      tempdens.resize((norb**2))
      write_bin_array(tempdens, "tempdens.bin")
      density2fchk(tdci_dir, nocc, norb, filename="tempdens.bin")
    else: # Don't modify density matrix
      density2fchk(tdci_dir, direction, i, nocc, norb)
    cubegen(f"out{i}")



def __main__():
  if not any([opts.cube_only, opts.render_only, opts.vmd_only, opts.ffmpeg_only]):
    do_all = True
  else:
    do_all = False

  gaussian_dir = args[0]
  tdci_dir = args[1]
  dir_idx = int(args[2])

  if opts.oneD: # One dimensional cubefile sampled along a line
    if opts.diff != 0:
      print("WARNING: NONZERO DIFF PARAMETER DETECTED IN 1-D INPUT -- THIS MAY NOT BE SUPPORTED")
    # Get sample direction from cli options or tdci input.
    thetaphi = [0.0, 0.0]
    if ((opts.theta is None) and (opts.phi is None)):
      polar_coords, nsteps, nprintstep, dt, Emax = parse_tdci_input(tdci_dir)
      thetaphi = polar_coords[dir_idx]
    else:
      if opts.theta != None:
        thetaphi[0] = float(opts.theta)
      if opts.phi != None:
        thetaphi[1] = float(opts.phi)

    if opts.cube_only or do_all:
      if opts.orblist == "":
        generate_cube_series_1D(gaussian_dir, tdci_dir, dir_idx, opts.start_timestep, opts.end_timestep, thetaphi, diff=opts.diff)
      else:
        orblist = [int(x) for x in opts.orblist.split(",") if x]
        generate_cube_series_1D_orbsubset(gaussian_dir, tdci_dir, dir_idx, opts.start_timestep, opts.end_timestep, thetaphi, orblist, diff=opts.diff)

    if opts.render_only or do_all:
      #plot_density_z_series(opts.start_timestep, opts.end_timestep)
      if opts.orblist == "":
        plot_density_1D_allinone(opts.start_timestep, opts.end_timestep, thetaphi)
      else:
        orblist = [int(x) for x in opts.orblist.split(",") if x]
        for orb in orblist:
          plot_density_1D_allinone(opts.start_timestep, opts.end_timestep, thetaphi,cubeprefix="out")
      
  else: # 3D cubefile
    if opts.cube_only or do_all:
      generate_cube_series(gaussian_dir, tdci_dir, dir_idx, opts.start_timestep, opts.end_timestep, diff=opts.diff)
    if opts.render_only or opts.vmd_only or do_all:
      vmd_cube2bmp(opts.start_timestep, opts.end_timestep, opts.isoval1, opts.isoval2)

    if opts.render_only or opts.ffmpeg_only or do_all:
      ffmpeg_bmp2mp4(opts.start_timestep)
    check_vmd_orphans()


__main__()








