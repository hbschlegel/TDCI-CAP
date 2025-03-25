
import os, sys, subprocess

# Tools for automating density movies


# Config variables for grid bounds (edit these)
# Angstroms
xmin = ( -16.0, -16.0, -16.0 ) # .cube region minimum corner
xmax = ( 16.0, 16.0, 48.0 )    # .cube region maximum corner
nsteps = 55                    # .cube grid points in each direction

# TDCI timesteps to create .cube files for
start_timestep_default = 1
end_timestep_default = 160


# Path to folder with executables
#tools = "/home/andy/tdci-test/TDCI-CAP/tools/" 
tools = os.path.dirname(os.path.abspath(__file__))
print(tools)



from optparse import OptionParser
usage = "usage: python3 movietool.py [options] [gaussian_directory] [tdci_directory] [field direction index]"
parser = OptionParser(usage=usage)
parser.add_option("--cube_only", action="store_true", default=False, dest="cube_only", help="Create .cube files from density in TDCI directory")
parser.add_option("--render_only", action="store_true", default=False, dest="render_only", help="Run VMD and FFMPEG")
parser.add_option("--vmd_only", action="store_true", default=False, dest="vmd_only", help="Run only VMD")
parser.add_option("--ffmpeg_only", action="store_true", default=False, dest="ffmpeg_only", help="Run only ffmpeg")
parser.add_option("--start_timestep", action="store", type="int", default=start_timestep_default, dest="start_timestep", help=f"Set the timestep to start on. Default {start_timestep_default}")
parser.add_option("--end_timestep", action="store", type="int", default=end_timestep_default, dest="end_timestep", help=f"Set the timestep to end on. Default {end_timestep_default}")
parser.add_option("--isoval1", action="store", type="float", default="0.002", dest="isoval1", help=f"Inner isovalue in VMD render. Default 0.002")
parser.add_option("--isoval2", action="store", type="float", default="0.00095", dest="isoval2", help=f"Outer isovalue in VMD render. Default 0.00095")
(opts, args) = parser.parse_args()



# Determine number of CPUs we can use for cube file generation.
p = subprocess.Popen("nproc", shell=True, stdout=subprocess.PIPE)
p.wait()
NCPU = int(p.stdout.read().decode().strip())
del p # remove this symbol from global namespace
if NCPU > 3:
  NCPU = NCPU - 1 # Let's leave one free ^_^




# assuming there is only one in the directory, return path of chk file
def find_chk_file(directory):
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
    self.Z_map = {'H': 1, 'He':2, 'C':6, 'N':7, 'O':8, 'F':9, 'Cl':17, 'Br':35, 'I':53 } 
    self.inv_Z_map = {value: key for key, value in self.Z_map.items()}
    # i assume this charge map is basis dependent, but it probably doesnt matter for
    #  plotting densities
    self.Chg_map = {'H':1.0, 'C':6.0, 'O':8.0, 'I':25.0 }
    self.Chg_Zmap = {1:1.0, 6:6.0, 8:8.0, 35:35.0, 53:25.0 }

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
      print( (z,chg,coords) )
      f.write(f"{z:>4d}   {chg:9.4f}   {atom[1]:9.4f}   {atom[2]:9.4f}   {atom[3]:9.4f}\n")
    # Normally the grid values follow, but not needed for this template
    f.close()
    return 0


# Call the density2fchk fortran program to insert the density from tdci into the template fchk file
def density2fchk(tdci_dir, direction, timestep, nstep, nocc, norb):
  #!: density2fchk.f90 Usage:
  #!:   density2fchk infile MOdensityfile nstep diff iscale
  #!:
  #!: Command line input:
  #!:    infile -  full name of the template formatted checkpoint file (with ".fchk")
  #!:    MOdensityfile - full name of the ion coefficient file (e.g. MO_density-e1-d1.dat)
  #!:    nstep - step number to be used for the constructiion of the NTOs and hole density
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

  diff = 0
  iscale = 0
  print(timestep)
  # density2fchk template_fchk density_binary, nstep, diff, iscale, nocc, norb
  densfile = f"{tdci_dir}/matrices/MO_density-e1-d{direction}.{timestep}00.bin"
  # iscale=0 : dont scale
  # iscale=1 : read from maxrate.dat and scale by 1./maxrate (disabled)
  p = subprocess.Popen(f"{tools}/density2fchk temp.fchk {densfile} {nstep} {diff} {iscale} {nocc} {norb} ", shell=True)
  p.wait()
  return 0

# Run cubegen using previously generated input.cube as a template for the region to sample 
def cubegen(timestep):
  cpus = NCPU
  ngrid = "-1"
  print(f"Generating cube : out{timestep}.cube")
  p = subprocess.Popen(f"cubegen {cpus} FDENSITY=SCF temp-1.fchk out{timestep}.cube {ngrid} h input.cube >> output", shell=True)
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
def generate_cube_series(gaussian_dir, tdci_dir, direction, start_timestep, end_timestep):
  nocc, norb = parse_tdci_output(tdci_dir)
  # make fchk and get coords from it
  chkpath = find_chk_file(gaussian_dir)  
  p = subprocess.Popen(f"cp {chkpath} temp.chk", shell=True)
  p.wait()
  p = subprocess.Popen(f"formchk temp.chk", shell=True)
  p.wait()

  # create input.cube template
  coords = read_fchk_coords("temp.fchk")
  cm = cubemaker( gaussian_dir, tdci_dir )
  cm.generate_template_cube( xmin, xmax, nsteps, coords )
  
  for i in range(start_timestep, end_timestep+1):
    density2fchk(tdci_dir, direction, i, nsteps, nocc, norb) 
    cubegen(i)
                                         
                                         
def __main__():
  if not any([opts.cube_only, opts.render_only, opts.vmd_only, opts.ffmpeg_only]):
    do_all = True
  else:
    do_all = False

  if opts.cube_only or do_all:
    #               gaussian_dir, tdci_dir, direction
    generate_cube_series(args[0], args[1], args[2], opts.start_timestep, opts.end_timestep)

  if opts.render_only or opts.vmd_only or do_all:
    vmd_cube2bmp(opts.start_timestep, opts.end_timestep, opts.isoval1, opts.isoval2)

  if opts.render_only or opts.ffmpeg_only or do_all:
    ffmpeg_bmp2mp4(opts.start_timestep)

  check_vmd_orphans()
    

__main__()








