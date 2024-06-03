import trimesh, os, glob
import numpy as np
import pandas as pd

def create_directory_structure(name, CLASS, varied_parameters, launch_directories):

    """
    Creates directories that contain computer experiments + saves them
    in a datafile so that it is easier to launch them using bash.
    """
    
    def include_attribute(obj, written, varied_parameters, name):
        index = 0
        for v in varied_parameters:
            if(hasattr(obj, v)) and written[index] == False:
                temp_name = v.replace("_", "")
                name += "_"+temp_name+"_{}".format(getattr(obj, v))
                written[index] = True
            index+=1

        return written, name
    counter = 0

    written = np.full(len(varied_parameters), False)
    written, name = include_attribute(CLASS, written, varied_parameters, name)

    new_path = name
    isExist = os.path.exists(new_path)
    if not isExist:
       os.makedirs(new_path)

    isExist = os.path.exists(new_path+"/checkpoints")
    if not isExist:
       os.makedirs(new_path+"/checkpoints")

    launch_directories.write(new_path)
    launch_directories.write('\n')
    return new_path

def prepare_bash_launch():

    """
    Considers how many simulations have been launched already and tells you
    what directory to pass to slurm script in order to launch new simulations.
    """
    
    bash_path = "bashdirs"
    isExist = os.path.exists(bash_path)
    if not isExist:
       os.makedirs(bash_path)

    # check how many directory files there are
    counter_directories = 1
    for files in glob.glob(bash_path+'/directories*'):
        counter_directories += 1

    print("LAUNCH ---> bashdirs/directories_"+str(counter_directories))
    launch_directories = open(bash_path+'/directories_'+str(counter_directories)+'.dat', 'w')
    return launch_directories

def check_num_runs(outkeyword):
    
    # check how many directory files there are
    counter_directories = 1
    for files in glob.glob(outkeyword+"/run*"):
        counter_directories += 1
    return counter_directories

def append_lines_to_list(text):
    """
    Appends each line from a block of text to a list.

    :param text: A string containing multiple lines of text.
    :return: A list of lines.
    """
    lines = []  # Initialize an empty list to store the lines
    for line in text.split('\n'):  # Split the text into lines
        lines.append(line)  # Append each line to the list
    return lines

class SimLauncher_Vesicle:

    """
    This class generates a Python script to launch a TriLMP vesicle.

    The default values of the class come from TriLMP simulations run
    in the past. Additional details can be found on Github.

    This class as well as the Python script it gives rise to
    have been written with a focus on readability.
    """

    def __init__(self, 
                 
                 # use trimesh to generate triangulation
                 use_trimesh = True,

                 # use in-house triangulation
                 mesh_coordinates_file=None,
                 mesh_faces_file=None,

                 # mesh properties
                 resolution = 4,
                 sigma = 1.0,
                 
                 # membrane properties
                 kappa_b = 20,
                 kappa_a = 2.5e5,
                 kappa_v = 0.0,
                 kappa_c = 0.0,
                 kappa_t = 1e4,
                 kappa_r = 1e3,
                 
                 # MD properties
                 step_size     = 0.001,
                 traj_steps    = 50,
                 langevin_damp = 1.0,
                 temperature   = 1.0,
                 langevin_seed = 123,
                 total_sim_time = 1000,
                 discret_snapshots = 10,
                 print_program_iterations = 10,
                 eq_rounds = 10,

                 # MC/TRIMEM bond flipping properties
                 flip_ratio=0.1,
                 switch_mode='random',

                 # simulation box
                 xlo = -30,
                 xhi = 30,
                 ylo = -30,
                 yhi = 30,
                 zlo = -30,
                 zhi = 30,
                 
                 # additional fixes
                 fix_COM = "#"
                 ):

        print("WARNING: Default membrane group is 'vertices'.")

        # basic lengthscale
        self.sigma = sigma

        self.mesh_coordinates_file = mesh_coordinates_file
        self.mesh_faces_file = mesh_faces_file
        
        # if you count on trimesh for initialization
        if use_trimesh:
                self.flag_use_trimesh=""
                self.flag_use_inhouse="#"
                # mesh properties
                self.resolution = resolution

                # mesh initialization
                self.mesh = trimesh.creation.icosphere(self.resolution)
        else:
             self.flag_use_inhouse=""
             self.flag_use_trimesh="#"
             mesh_coordinates = pd.read_csv(mesh_coordinates_file, header = None, index_col = False, sep = ' ')
             mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
             mesh_faces = pd.read_csv(mesh_faces_file, header = None, index_col = False, sep = ' ')
             mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
             self.mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

        self.N = len(self.mesh.vertices)

        # rescaling mesh distances
        desired_average_distance = 2**(1.0/6.0) * self.sigma
        current_average_distance = np.mean(self.mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        self.mesh.vertices *=scaling
        
        # extract radius membrane
        self.radius_membrane = np.mean(np.sqrt(self.mesh.vertices[:, 0]**2 + self.mesh.vertices[:, 1]**2 + self.mesh.vertices[:, 2]**2))

        # leftmost vesicle
        self.edge_vesicle = np.max(self.mesh.vertices[:, 0])

        # membrane mechanical properties
        self.kappa_b = kappa_b
        self.kappa_a = kappa_a
        self.kappa_v = kappa_v
        self.kappa_c = kappa_c
        self.kappa_t = kappa_t
        self.kappa_r = kappa_r

        # MD properties
        self.step_size     = step_size
        self.traj_steps    = traj_steps
        self.langevin_damp = langevin_damp
        self.temperature   = temperature
        self.langevin_seed = langevin_seed

        # MC/TRIMEM bond flipping properties
        self.flip_ratio=flip_ratio
        self.switch_mode=switch_mode

        # simulation step structure (given in time units)
        self.total_sim_time = total_sim_time
        self.MD_simulation_steps = int(self.total_sim_time/self.step_size)
        
        # ouput and printing (given in time units)
        self.discret_snapshots = discret_snapshots
        self.print_frequency = int(self.discret_snapshots/self.step_size)
        self.print_program_iterations = print_program_iterations

        # simulation box
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

        # equilibration steps for the membrane (how many MD stages need to pass)
        self.eq_rounds = eq_rounds
        self.equilibration_steps = self.eq_rounds*self.traj_steps

        # additional fixes
        self.fix_COM = fix_COM

    def write_function(self, filename):
        
        text_to_write = f"""# ---------------------------------------------------------------------#
# TriLMP vesicle                                                       #
# Author: Maitane Muñoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at) #
#                                                                      #
# This code launches TriLMP (clean version) for a fluid membrane.      #
# ---------------------------------------------------------------------#

import trimesh
import numpy as np
import pandas as pd
from trimem.mc.trilmp import TriLmp

# mesh initialization
{self.flag_use_trimesh}mesh = trimesh.creation.icosphere({self.resolution})
{self.flag_use_inhouse}mesh_coordinates = pd.read_csv({self.mesh_coordinates_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
{self.flag_use_inhouse}mesh_faces = pd.read_csv({self.mesh_faces_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
{self.flag_use_inhouse}mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

N = len(mesh.vertices)

# rescaling mesh distances
desired_average_distance = 2**(1.0/6.0) * {self.sigma}
current_average_distance = np.mean(mesh.edges_unique_length)
scaling = desired_average_distance/current_average_distance
mesh.vertices *=scaling

print(f"MESH VERTICES : ", len(mesh.vertices))
print(f"MESH FACES    : ", len(mesh.faces))
print(f"MESH EDGES    : ", len(mesh.edges))

# initialization of the trilmp object
trilmp=TriLmp(initialize=True,                            # use mesh to initialize mesh reference
            debug_mode=False,                             # DEBUGGING: print everything
            num_particle_types=1,                         # PART. SPECIES: total particle species in system 
            mass_particle_type=[1.0],                     # PART. SPECIES: mass of species in system
            group_particle_type=['vertices'],             # PART. SPECIES: group names for species in system

            mesh_points=mesh.vertices,                  # input mesh vertices 
            mesh_faces=mesh.faces,                      # input of the mesh faces

            kappa_b={self.kappa_b},                            # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a={self.kappa_a},                            # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v={self.kappa_v},                            # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c={self.kappa_c},                            # MEMBRANE MECHANICS: constraint on area difference change (kB T)
            kappa_t={self.kappa_t},                            # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r={self.kappa_r},                            # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            
            step_size={self.step_size},                        # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps={self.traj_steps},                      # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio={self.flip_ratio},                      # MC PART SIMULATION: fraction of edges to flip
            initial_temperature={self.temperature},            # MD PART SIMULATION: temperature of the system
            pure_MD=True,                                   # MD PART SIMULATION: accept every MD trajectory
            switch_mode="random",                    # MD/MC PART SIMULATION: 'random' or 'alternating' flip-or-move
            box=({self.xlo}, {self.xhi}, {self.ylo}, {self.yhi}, {self.zlo}, {self.zhi}),         # MD PART SIMULATION: simulation box properties, periodic

            equilibration_rounds={self.equilibration_steps},   # MD PART SIMULATION: HOW LONG DO WE LET THE MEMBRANE EQUILIBRATE
            
            info={self.print_frequency},                              # OUTPUT: frequency output in shell
            thin={self.print_frequency},                              # OUTPUT: frequency trajectory output
            performance_increment={self.print_program_iterations},    # OUTPUT: output performace stats to prefix_performance.dat file - PRINTED MD+MC FREQUENCY
            energy_increment={self.print_program_iterations},         # OUTPUT: output energies to energies.dat file - PRINTED MD FREQUENCY
            checkpoint_every=100*{self.print_frequency},              # OUTPUT: interval of checkpoints (alternating pickles) - PRINTED MD+MC FREQUENCY
            )

# -------------------------#
#  LAMMPS MODIFICATIONS    #
# -------------------------#

# .................................................
#                 PAIR STYLES
# .................................................

# cleanup pair style in case
trilmp.lmp.command("pair_style none")

# pair interactions
trilmp.lmp.command(f"pair_style hybrid/overlay table linear 2000 cosine/squared 1.5")

# compulsory lines
trilmp.lmp.command("pair_modify pair table special lj/coul 0.0 0.0 0.0 tail no")
trilmp.lmp.command("pair_coeff 1 1 table trimem_srp.table trimem_srp")

# set all interactions to zero just in case for added potentials
trilmp.lmp.command("pair_coeff * * cosine/squared 0.0 0.0")

# .................................................
#         COMPUTES, FIXES, ETC
# .................................................

# dump particle trajectories (vertex coordinates)
trilmp.lmp.command(f"dump XYZ all custom {self.print_frequency} trajectory.gz id type x y z")

# compute potential energy
trilmp.lmp.command("compute PeMembrane vertices pe/atom pair")
trilmp.lmp.command("compute pe vertices reduce sum c_PeMembrane")

# compute position CM vesicle
trilmp.lmp.command("compute MembraneCOM vertices com")

# compute shape of the vesicle
trilmp.lmp.command("compute RadiusGMem vertices gyration")
trilmp.lmp.command("compute MemShape vertices gyration/shape RadiusGMem")

# compute temperature of the vesicle
trilmp.lmp.command("compute TempComputeMem vertices temp")

# print out all the computations
trilmp.lmp.command(f"fix  aveMEM all ave/time {self.print_frequency} 1 {self.print_frequency} c_TempComputeMem c_pe c_MembraneCOM[1] c_MembraneCOM[2] c_MembraneCOM[3] c_MemShape[1] c_MemShape[2] c_MemShape[3] c_MemShape[4] c_MemShape[5] c_MemShape[6] file 'membrane_CM.dat'")

# .................................................#
#       PRE-EQUILIBRATION INTEGRATION              #
# .................................................#

# include the integrators (pre-equilibration)
trilmp.lmp.command("fix NVEMEM vertices nve")
trilmp.lmp.command(f"fix LGVMEM vertices langevin {self.temperature} {self.temperature} {self.langevin_damp} {self.langevin_seed} zero yes")

# fix the CM of the vesicle - by default at the center
{self.fix_COM}trilmp.lmp.command(f"fix COMFIX vertices recenter 0 0 0")

# !!!!!!!!!!!!!!!!!!!!!!!!!#
# -------------------------#
#    POST-EQUILIBRATION    #
# -------------------------#
# !!!!!!!!!!!!!!!!!!!!!!!!!#

postequilibration_commands = []

# cleanup of the fixes
postequilibration_commands.append("unfix NVEMEM")
postequilibration_commands.append("unfix LGVMEM")
{self.fix_COM}postequilibration_commands.append("unfix COMFIX")

# .................................................#
#  INTEGRATION OF EQS OF MOTION (All beads)        #
# .................................................#

postequilibration_commands.append("fix NVEMEM all nve")
postequilibration_commands.append(f"fix LGVMEM all langevin {self.temperature} {self.temperature} {self.langevin_damp} {33+self.langevin_seed} zero yes")
{self.fix_COM}postequilibration_commands.append(f"fix COMFIX vertices recenter 0 0 0")

# -------------------------#
#         RUN              #
# -------------------------#

# RUN THE SIMULATION
trilmp.run({self.MD_simulation_steps}, integrators_defined=True, fix_symbionts_near=False, 
        postequilibration_lammps_commands = postequilibration_commands)

print("End of the simulation.")

        """

        file = open(filename, "w")
        file.write(text_to_write)
        file.close()
        print("File has been produced")

class SimLauncher_BasicGradient:

    """
    This class generates a Python script to launch TriLMP together
    with other LAMMPS functionalities. The main usage of the class
    is to place the vesicle in a metabolite gradient and see
    how it deforms.

    The default values of the class come from TriLMP simulations run
    in the past. Additional details can be found on Github.

    This class as well as the Python script it gives rise to
    have been written with a focus on readability.
    """

    def __init__(self, 
                 
                 # use trimesh to generate triangulation
                 use_trimesh = True,

                 # use in-house triangulation
                 mesh_coordinates_file=None,
                 mesh_faces_file=None,

                 # mesh properties
                 resolution = 4,
                 sigma = 1.0,
                 
                 # membrane properties
                 kappa_b = 20,
                 kappa_a = 2.5e5,
                 kappa_v = 0.0,
                 kappa_c = 0.0,
                 kappa_t = 1e4,
                 kappa_r = 1e3,
                 
                 # MD properties
                 step_size     = 0.001,
                 traj_steps    = 50,
                 langevin_damp = 1.0,
                 temperature   = 1.0,
                 langevin_seed = 123,
                 total_sim_time = 1000,
                 discret_snapshots = 10,
                 print_program_iterations = 10,
                 eq_rounds = 10,

                 # MC/TRIMEM bond flipping properties
                 flip_ratio=0.1,
                 switch_mode='random',

                 # simulation box
                 xlo = -30,
                 xhi = 30,
                 ylo = -30,
                 yhi = 30,
                 zlo = -30,
                 zhi = 30,
                 
                 # bead properties
                 mass_metabolites   = 1.0,
                 sigma_metabolites  = 1.0,

                 # additional fixes
                 fix_COM = "#",
                 channel_height = 0.2
                 ):

        print("WARNING: Default membrane group is 'vertices' and metabolites 'metabolites'.")

        # basic lengthscale
        self.sigma = sigma

        self.mesh_coordinates_file = mesh_coordinates_file
        self.mesh_faces_file = mesh_faces_file
        
        # if you count on trimesh for initialization
        if use_trimesh:
                self.flag_use_trimesh=""
                self.flag_use_inhouse="#"
                # mesh properties
                self.resolution = resolution

                # mesh initialization
                self.mesh = trimesh.creation.icosphere(self.resolution)
        else:
             self.flag_use_inhouse=""
             self.flag_use_trimesh="#"
             mesh_coordinates = pd.read_csv(mesh_coordinates_file, header = None, index_col = False, sep = ' ')
             mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
             mesh_faces = pd.read_csv(mesh_faces_file, header = None, index_col = False, sep = ' ')
             mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
             self.mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

        self.N = len(self.mesh.vertices)
        # channel height, this is only here for filewriting
        self.channel_height = channel_height
        # rescaling mesh distances
        desired_average_distance = 2**(1.0/6.0) * self.sigma
        current_average_distance = np.mean(self.mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        self.mesh.vertices *=scaling
        
        # extract radius membrane
        self.radius_membrane = np.mean(np.sqrt(self.mesh.vertices[:, 0]**2 + self.mesh.vertices[:, 1]**2 + self.mesh.vertices[:, 2]**2))

        # leftmost vesicle
        self.edge_vesicle = np.max(self.mesh.vertices[:, 0])

        # membrane mechanical properties
        self.kappa_b = kappa_b
        self.kappa_a = kappa_a
        self.kappa_v = kappa_v
        self.kappa_c = kappa_c
        self.kappa_t = kappa_t
        self.kappa_r = kappa_r

        # MD properties
        self.step_size     = step_size
        self.traj_steps    = traj_steps
        self.langevin_damp = langevin_damp
        self.temperature   = temperature
        self.langevin_seed = langevin_seed

        # MC/TRIMEM bond flipping properties
        self.flip_ratio=flip_ratio
        self.switch_mode=switch_mode

        # simulation step structure (given in time units)
        self.total_sim_time = total_sim_time
        self.MD_simulation_steps = int(self.total_sim_time/self.step_size)
        
        # ouput and printing (given in time units)
        self.discret_snapshots = discret_snapshots
        self.print_frequency = int(self.discret_snapshots/self.step_size)
        self.print_program_iterations = print_program_iterations

        # simulation box
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

        # equilibration steps for the membrane (how many MD stages need to pass)
        self.eq_rounds = eq_rounds
        self.equilibration_steps = self.eq_rounds*self.traj_steps

        # bead properties
        self.mass_metabolites   = mass_metabolites
        self.sigma_metabolites  = sigma_metabolites
        self.sigma_tilde = 0.5*(self.sigma+self.sigma_metabolites)

        # additional fixes
        self.fix_COM = fix_COM

    def initialize_commands_gradient(self, params):
        
        def add_blocks_commands(params, params_keyword):
                
            lines = []
            for cmd in params[params_keyword]:
                lines.append(f"postequilibration_commands.append('{cmd}')")
            return lines
        
        self.interaction_membrane_metabolites = params['interaction_membrane_metabolites']
        self.rc_membrane_metabolites = params['rc_membrane_metabolites']
        self.rc_tilde_membrane_metabolites = self.rc_membrane_metabolites*self.sigma_tilde

        # add blocks of text
        self.metabolite_command_block  = "\n".join(add_blocks_commands(params, 'gradient_commands'))
        self.additional_external_fixes = "\n".join(add_blocks_commands(params, 'additional_fixes'))

    def write_function(self, filename):
        
        text_to_write = f"""# ---------------------------------------------------------------------#
# TriLMP vesicle in metabolite gradients with various geometries       #
# Author: Maitane Muñoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at) #
#                                                                      #
# This code launches TriLMP (clean version) for a fluid membrane.      #
# Additionally, we consider gradient sink and sources in the system.   #                                                                      #
# ---------------------------------------------------------------------#

import trimesh
import numpy as np
import pandas as pd
from trimem.mc.trilmp import TriLmp

# mesh initialization
{self.flag_use_trimesh}mesh = trimesh.creation.icosphere({self.resolution})
{self.flag_use_inhouse}mesh_coordinates = pd.read_csv({self.mesh_coordinates_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
{self.flag_use_inhouse}mesh_faces = pd.read_csv({self.mesh_faces_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
{self.flag_use_inhouse}mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

N = len(mesh.vertices)

# rescaling mesh distances
desired_average_distance = 2**(1.0/6.0) * {self.sigma}
current_average_distance = np.mean(mesh.edges_unique_length)
scaling = desired_average_distance/current_average_distance
mesh.vertices *=scaling

print(f"MESH VERTICES : ", len(mesh.vertices))
print(f"MESH FACES    : ", len(mesh.faces))
print(f"MESH EDGES    : ", len(mesh.edges))

# initialization of the trilmp object
trilmp=TriLmp(initialize=True,                            # use mesh to initialize mesh reference
            debug_mode=False,                             # DEBUGGING: print everything
            num_particle_types=2,                         # PART. SPECIES: total particle species in system 
            mass_particle_type=[1.0, {self.mass_metabolites}],          # PART. SPECIES: mass of species in system
            group_particle_type=['vertices', 'metabolites'],     # PART. SPECIES: group names for species in system

            mesh_points=mesh.vertices,                  # input mesh vertices 
            mesh_faces=mesh.faces,                      # input of the mesh faces

            kappa_b={self.kappa_b},                            # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a={self.kappa_a},                            # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v={self.kappa_v},                            # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c={self.kappa_c},                            # MEMBRANE MECHANICS: constraint on area difference change (kB T)
            kappa_t={self.kappa_t},                            # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r={self.kappa_r},                            # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            
            step_size={self.step_size},                        # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps={self.traj_steps},                      # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio={self.flip_ratio},                      # MC PART SIMULATION: fraction of edges to flip
            initial_temperature={self.temperature},            # MD PART SIMULATION: temperature of the system
            pure_MD=True,                                   # MD PART SIMULATION: accept every MD trajectory
            switch_mode="random",                    # MD/MC PART SIMULATION: 'random' or 'alternating' flip-or-move
            box=({self.xlo}, {self.xhi}, {self.ylo}, {self.yhi}, {self.zlo}, {self.zhi}),         # MD PART SIMULATION: simulation box properties, periodic

            equilibration_rounds={self.equilibration_steps},   # MD PART SIMULATION: HOW LONG DO WE LET THE MEMBRANE EQUILIBRATE
            
            info={self.print_frequency},                              # OUTPUT: frequency output in shell
            thin={self.print_frequency},                              # OUTPUT: frequency trajectory output
            performance_increment={self.print_program_iterations},    # OUTPUT: output performace stats to prefix_performance.dat file - PRINTED MD+MC FREQUENCY
            energy_increment={self.print_program_iterations},         # OUTPUT: output energies to energies.dat file - PRINTED MD FREQUENCY
            checkpoint_every=100*{self.print_frequency},              # OUTPUT: interval of checkpoints (alternating pickles) - PRINTED MD+MC FREQUENCY
            )

# -------------------------#
#  LAMMPS MODIFICATIONS    #
# -------------------------#

# set masses (do two methods)
trilmp.lmp.command("mass 1 1")
trilmp.lmp.command("mass 2 1")

# .................................................
#                 PAIR STYLES
# .................................................

# cleanup pair style in case
trilmp.lmp.command("pair_style none")

# pair interactions
trilmp.lmp.command(f"pair_style hybrid/overlay table linear 2000 cosine/squared 1.5")

# compulsory lines
trilmp.lmp.command("pair_modify pair table special lj/coul 0.0 0.0 0.0 tail no")
trilmp.lmp.command("pair_coeff 1 1 table trimem_srp.table trimem_srp")

# set all interactions to zero just in case for added potentials
trilmp.lmp.command("pair_coeff * * cosine/squared 0.0 0.0")

# pair coefficients for cosine square interaction
trilmp.lmp.command(f"pair_coeff 1 2 cosine/squared {self.interaction_membrane_metabolites} 0 {self.rc_tilde_membrane_metabolites}")

# .................................................
#         COMPUTES, FIXES, ETC
# .................................................

# dump particle trajectories (vertex coordinates)
trilmp.lmp.command(f"dump XYZ all custom {self.print_frequency} trajectory.gz id type x y z")

# compute potential energy
trilmp.lmp.command("compute PeMembrane vertices pe/atom pair")
trilmp.lmp.command("compute pe vertices reduce sum c_PeMembrane")

# compute position CM vesicle
trilmp.lmp.command("compute MembraneCOM vertices com")

# compute shape of the vesicle
trilmp.lmp.command("compute RadiusGMem vertices gyration")
trilmp.lmp.command("compute MemShape vertices gyration/shape RadiusGMem")

# compute temperature of the vesicle
trilmp.lmp.command("compute TempComputeMem vertices temp")

# print out all the computations
trilmp.lmp.command(f"fix  aveMEM all ave/time {self.print_frequency} 1 {self.print_frequency} c_TempComputeMem c_pe c_MembraneCOM[1] c_MembraneCOM[2] c_MembraneCOM[3] c_MemShape[1] c_MemShape[2] c_MemShape[3] c_MemShape[4] c_MemShape[5] c_MemShape[6] file 'membrane_CM.dat'")

# .................................................#
#       PRE-EQUILIBRATION INTEGRATION              #
# .................................................#

# include the integrators (pre-equilibration)
trilmp.lmp.command("fix NVEMEM vertices nve")
trilmp.lmp.command(f"fix LGVMEM vertices langevin {self.temperature} {self.temperature} {self.langevin_damp} {self.langevin_seed} zero yes")

# fix the CM of the vesicle - by default at the center
{self.fix_COM}trilmp.lmp.command(f"fix COMFIX vertices recenter 0 0 0")

# !!!!!!!!!!!!!!!!!!!!!!!!!#
# -------------------------#
#    POST-EQUILIBRATION    #
# -------------------------#
# !!!!!!!!!!!!!!!!!!!!!!!!!#

postequilibration_commands = []

# cleanup of the fixes
postequilibration_commands.append("unfix NVEMEM")
postequilibration_commands.append("unfix LGVMEM")
{self.fix_COM}postequilibration_commands.append("unfix COMFIX")

# .................................................#
#           INSERTION OF METABOLITES               #
#        (geometry depends on problem)             #
# .................................................#

{self.metabolite_command_block}

# .................................................#
#           ADDITIONAL EXTERNAL FIXES              #
# .................................................#

{self.additional_external_fixes}

# .................................................#
#  INTEGRATION OF EQS OF MOTION (All beads)        #
# .................................................#

postequilibration_commands.append("fix NVEMEM all nve")
postequilibration_commands.append(f"fix LGVMEM all langevin {self.temperature} {self.temperature} {self.langevin_damp} {33+self.langevin_seed} zero yes")
{self.fix_COM}postequilibration_commands.append(f"fix COMFIX vertices recenter 0 0 0")

# -------------------------#
#         RUN              #
# -------------------------#

# RUN THE SIMULATION
trilmp.run({self.MD_simulation_steps}, integrators_defined=True, fix_symbionts_near=False, 
        postequilibration_lammps_commands = postequilibration_commands)

print("End of the simulation.")

        """

        file = open(filename, "w")
        file.write(text_to_write)
        file.close()
        print("File has been produced")

class SimLauncher_HeterogeneousMembrane:

    """
    This class generates a Python script to launch TriLMP together
    with other LAMMPS functionalities. The main usage of the class
    is to place the vesicle in a metabolite gradient and see
    how it deforms, given that the membrane is heterogeneous;
    that is, composed of two different particle types.

    Within the membrane both types are identical. Differences
    between them come when it comes to interacting with the
    metabolites.
    """

    def __init__(self, 
                 
                 # main feature: heterogeneous membrane
                 heterogeneous_membrane = True,
                 heterogeneous_membrane_id = None,
                 fraction_heterogeneous = 0,
                 randomly_distribute_patches = False,

                 # use trimesh to generate triangulation
                 use_trimesh = True,

                 # use in-house triangulation
                 mesh_coordinates_file=None,
                 mesh_faces_file=None,

                 # mesh properties
                 resolution = 4,
                 sigma = 1.0,
                 
                 # membrane properties
                 kappa_b = 20,
                 kappa_a = 2.5e5,
                 kappa_v = 0.0,
                 kappa_c = 0.0,
                 kappa_t = 1e4,
                 kappa_r = 1e3,
                 
                 # MD properties
                 step_size     = 0.001,
                 traj_steps    = 50,
                 langevin_damp = 1.0,
                 temperature   = 1.0,
                 langevin_seed = 123,
                 total_sim_time = 1000,
                 discret_snapshots = 10,
                 print_program_iterations = 10,
                 eq_rounds = 10,

                 # MC/TRIMEM bond flipping properties
                 flip_ratio=0.1,
                 switch_mode='random',

                 # simulation box
                 xlo = -30,
                 xhi = 30,
                 ylo = -30,
                 yhi = 30,
                 zlo = -30,
                 zhi = 30,
                 
                 # bead properties
                 mass_metabolites   = 1.0,
                 sigma_metabolites  = 1.0,

                 # additional fixes
                 fix_COM = "#",
                 channel_height = 0.2,
                 probability = 0.
                 ):

        print("WARNING: Heterogeneous membrane simulations. Default membrane group (1 + 2) is 'vertices' and metabolites (3) 'metabolites'.")
        print("Particle type 1 is 'inert_vertices' and particle type 2 is 'active_vertices'.")
        
        # main feature: heterogeneous membrane
        self.heterogeneous_membrane = heterogeneous_membrane
        self.heterogeneous_membrane_id = heterogeneous_membrane_id

        # basic lengthscale
        self.sigma = sigma

        self.mesh_coordinates_file = mesh_coordinates_file
        self.mesh_faces_file = mesh_faces_file
        
        # if you count on trimesh for initialization
        if use_trimesh:
                self.flag_use_trimesh=""
                self.flag_use_inhouse="#"
                # mesh properties
                self.resolution = resolution

                # mesh initialization
                self.mesh = trimesh.creation.icosphere(self.resolution)
        else:
             self.flag_use_inhouse=""
             self.flag_use_trimesh="#"
             mesh_coordinates = pd.read_csv(mesh_coordinates_file, header = None, index_col = False, sep = ' ')
             mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
             mesh_faces = pd.read_csv(mesh_faces_file, header = None, index_col = False, sep = ' ')
             mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
             self.mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

        self.N = len(self.mesh.vertices)
        self.fraction_heterogeneous = fraction_heterogeneous
        self.N_patch = int(self.N*self.fraction_heterogeneous)

        # channel height, this is only here for filewriting
        self.channel_height = channel_height
        # transition probability, also only for filewriting
        self.probability = probability

        if randomly_distribute_patches:
            self.heterogeneous_membrane_id = np.random.randint(0, self.N, self.N_patch)
            self.heterogeneous_membrane_id = self.heterogeneous_membrane_id.tolist()
        elif not randomly_distribute_patches:
            index_max = np.where(self.mesh.vertices[:, 0] == np.max(self.mesh.vertices[:, 0]))[0]
            distances_to_point = np.sqrt((self.mesh.vertices[:, 0]-self.mesh.vertices[index_max, 0])**2 + (self.mesh.vertices[:, 1]-self.mesh.vertices[index_max, 1])**2 + (self.mesh.vertices[:, 2]-self.mesh.vertices[index_max, 2])**2)
            sorted_distances = np.argsort(distances_to_point)
            self.heterogeneous_membrane_id = sorted_distances[:int(self.fraction_heterogeneous*self.N)+1].tolist()

        # rescaling mesh distances
        desired_average_distance = 2**(1.0/6.0) * self.sigma
        current_average_distance = np.mean(self.mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        self.mesh.vertices *=scaling
        
        # extract radius membrane
        self.radius_membrane = np.mean(np.sqrt(self.mesh.vertices[:, 0]**2 + self.mesh.vertices[:, 1]**2 + self.mesh.vertices[:, 2]**2))

        # leftmost vesicle
        self.edge_vesicle = np.max(self.mesh.vertices[:, 0])

        # membrane mechanical properties
        self.kappa_b = kappa_b
        self.kappa_a = kappa_a
        self.kappa_v = kappa_v
        self.kappa_c = kappa_c
        self.kappa_t = kappa_t
        self.kappa_r = kappa_r

        # MD properties
        self.step_size     = step_size
        self.traj_steps    = traj_steps
        self.langevin_damp = langevin_damp
        self.temperature   = temperature
        self.langevin_seed = langevin_seed

        # MC/TRIMEM bond flipping properties
        self.flip_ratio=flip_ratio
        self.switch_mode=switch_mode

        # simulation step structure (given in time units)
        self.total_sim_time = total_sim_time
        self.MD_simulation_steps = int(self.total_sim_time/self.step_size)
        
        # ouput and printing (given in time units)
        self.discret_snapshots = discret_snapshots
        self.print_frequency = int(self.discret_snapshots/self.step_size)
        self.print_program_iterations = print_program_iterations

        # simulation box
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

        # equilibration steps for the membrane (how many MD stages need to pass)
        self.eq_rounds = eq_rounds
        self.equilibration_steps = self.eq_rounds*self.traj_steps

        # bead properties
        self.mass_metabolites   = mass_metabolites
        self.sigma_metabolites  = sigma_metabolites
        self.sigma_tilde = 0.5*(self.sigma+self.sigma_metabolites)

        # additional fixes
        self.fix_COM = fix_COM

    def initialize_commands_gradient(self, params):
        
        def add_blocks_commands(params, params_keyword):
                
            lines = []
            for cmd in params[params_keyword]:
                lines.append(f"postequilibration_commands.append('{cmd}')")
            return lines
        
        self.interaction_membrane_metabolites = params['interaction_membrane_metabolites']
        self.rc_membrane_metabolites = params['rc_membrane_metabolites']
        self.rc_tilde_membrane_metabolites = self.rc_membrane_metabolites*self.sigma_tilde

        # add blocks of text
        self.metabolite_command_block  = "\n".join(add_blocks_commands(params, 'gradient_commands'))
        self.additional_external_fixes = "\n".join(add_blocks_commands(params, 'additional_fixes'))

    def write_function(self, filename):
        
        text_to_write = f"""# ---------------------------------------------------------------------#
# TriLMP vesicle in metabolite gradients with various geometries       #
# Author: Maitane Muñoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at) #
#                                                                      #
# This code launches TriLMP (clean version) for a fluid membrane.      #
# Additionally, we consider gradient sink and sources in the system.   #                                                                      #
# ---------------------------------------------------------------------#

import trimesh
import numpy as np
import pandas as pd
from trimem.mc.trilmp import TriLmp

# mesh initialization
{self.flag_use_trimesh}mesh = trimesh.creation.icosphere({self.resolution})
{self.flag_use_inhouse}mesh_coordinates = pd.read_csv({self.mesh_coordinates_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
{self.flag_use_inhouse}mesh_faces = pd.read_csv({self.mesh_faces_file}, header = None, index_col = False, sep = ' ')
{self.flag_use_inhouse}mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
{self.flag_use_inhouse}mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array) 

N = len(mesh.vertices)

# rescaling mesh distances
desired_average_distance = 2**(1.0/6.0) * {self.sigma}
current_average_distance = np.mean(mesh.edges_unique_length)
scaling = desired_average_distance/current_average_distance
mesh.vertices *=scaling

print(f"MESH VERTICES : ", len(mesh.vertices))
print(f"MESH FACES    : ", len(mesh.faces))
print(f"MESH EDGES    : ", len(mesh.edges))

# initialization of the trilmp object
trilmp=TriLmp(initialize=True,                            # use mesh to initialize mesh reference
            debug_mode=False,                             # DEBUGGING: print everything
            heterogeneous_membrane = {self.heterogeneous_membrane},
            heterogeneous_membrane_id = {self.heterogeneous_membrane_id},
            num_particle_types=3,                         # PART. SPECIES: total particle species in system 
            mass_particle_type=[1.0, 1.0, {self.mass_metabolites}],          # PART. SPECIES: mass of species in system
            group_particle_type=['inert_vertices', 'active_vertices','metabolites'],     # PART. SPECIES: group names for species in system

            mesh_points=mesh.vertices,                  # input mesh vertices 
            mesh_faces=mesh.faces,                      # input of the mesh faces

            kappa_b={self.kappa_b},                            # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a={self.kappa_a},                            # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v={self.kappa_v},                            # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c={self.kappa_c},                            # MEMBRANE MECHANICS: constraint on area difference change (kB T)
            kappa_t={self.kappa_t},                            # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r={self.kappa_r},                            # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            
            step_size={self.step_size},                        # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps={self.traj_steps},                      # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio={self.flip_ratio},                      # MC PART SIMULATION: fraction of edges to flip
            initial_temperature={self.temperature},            # MD PART SIMULATION: temperature of the system
            pure_MD=True,                                   # MD PART SIMULATION: accept every MD trajectory
            switch_mode="random",                    # MD/MC PART SIMULATION: 'random' or 'alternating' flip-or-move
            box=({self.xlo}, {self.xhi}, {self.ylo}, {self.yhi}, {self.zlo}, {self.zhi}),         # MD PART SIMULATION: simulation box properties, periodic

            equilibration_rounds={self.equilibration_steps},   # MD PART SIMULATION: HOW LONG DO WE LET THE MEMBRANE EQUILIBRATE
            
            info={self.print_frequency},                              # OUTPUT: frequency output in shell
            thin={self.print_frequency},                              # OUTPUT: frequency trajectory output
            performance_increment={self.print_program_iterations},    # OUTPUT: output performace stats to prefix_performance.dat file - PRINTED MD+MC FREQUENCY
            energy_increment={self.print_program_iterations},         # OUTPUT: output energies to energies.dat file - PRINTED MD FREQUENCY
            checkpoint_every=100*{self.print_frequency},              # OUTPUT: interval of checkpoints (alternating pickles) - PRINTED MD+MC FREQUENCY
            )

# -------------------------#
#  LAMMPS MODIFICATIONS    #
# -------------------------#

trilmp.lmp.command("group vertices type 1 2")

# .................................................
#                 PAIR STYLES
# .................................................

# cleanup pair style in case
trilmp.lmp.command("pair_style none")

# pair interactions
trilmp.lmp.command(f"pair_style hybrid/overlay table linear 2000 cosine/squared 1.5")

# compulsory lines
trilmp.lmp.command("pair_modify pair table special lj/coul 0.0 0.0 0.0 tail no")
trilmp.lmp.command("pair_coeff 1 1 table trimem_srp.table trimem_srp11")
trilmp.lmp.command("pair_coeff 2 2 table trimem_srp.table trimem_srp22")
trilmp.lmp.command("pair_coeff 1 2 table trimem_srp.table trimem_srp12")

# set all interactions to zero just in case for added potentials
trilmp.lmp.command("pair_coeff * * cosine/squared 0.0 0.0")

# pair coefficients for cosine square interaction
trilmp.lmp.command(f"pair_coeff 2 3 cosine/squared {self.interaction_membrane_metabolites} 0 {self.rc_tilde_membrane_metabolites}")

# .................................................
#         COMPUTES, FIXES, ETC
# .................................................

# dump particle trajectories (vertex coordinates)
trilmp.lmp.command(f"dump XYZ all custom {self.print_frequency} trajectory.gz id type x y z")

# compute potential energy
trilmp.lmp.command("compute PeMembrane vertices pe/atom pair")
trilmp.lmp.command("compute pe vertices reduce sum c_PeMembrane")

# compute position CM vesicle
trilmp.lmp.command("compute MembraneCOM vertices com")

# compute shape of the vesicle
trilmp.lmp.command("compute RadiusGMem vertices gyration")
trilmp.lmp.command("compute MemShape vertices gyration/shape RadiusGMem")

# compute temperature of the vesicle
trilmp.lmp.command("compute TempComputeMem vertices temp")

# print out all the computations
trilmp.lmp.command(f"fix  aveMEM all ave/time {self.print_frequency} 1 {self.print_frequency} c_TempComputeMem c_pe c_MembraneCOM[1] c_MembraneCOM[2] c_MembraneCOM[3] c_MemShape[1] c_MemShape[2] c_MemShape[3] c_MemShape[4] c_MemShape[5] c_MemShape[6] file 'membrane_CM.dat'")

# .................................................#
#       PRE-EQUILIBRATION INTEGRATION              #
# .................................................#

# include the integrators (pre-equilibration)
trilmp.lmp.command("fix NVEMEM vertices nve")
trilmp.lmp.command(f"fix LGVMEM vertices langevin {self.temperature} {self.temperature} {self.langevin_damp} {self.langevin_seed} zero yes")

# fix the CM of the vesicle - by default at the center
{self.fix_COM}trilmp.lmp.command(f"fix COMFIX vertices recenter 0 0 0")

# !!!!!!!!!!!!!!!!!!!!!!!!!#
# -------------------------#
#    POST-EQUILIBRATION    #
# -------------------------#
# !!!!!!!!!!!!!!!!!!!!!!!!!#

postequilibration_commands = []

# cleanup of the fixes
postequilibration_commands.append("unfix NVEMEM")
postequilibration_commands.append("unfix LGVMEM")
{self.fix_COM}postequilibration_commands.append("unfix COMFIX")

# .................................................#
#           INSERTION OF METABOLITES               #
#        (geometry depends on problem)             #
# .................................................#

{self.metabolite_command_block}

# .................................................#
#           ADDITIONAL EXTERNAL FIXES              #
# .................................................#

{self.additional_external_fixes}

# .................................................#
#  INTEGRATION OF EQS OF MOTION (All beads)        #
# .................................................#

postequilibration_commands.append("fix NVEMEM all nve")
postequilibration_commands.append(f"fix LGVMEM all langevin {self.temperature} {self.temperature} {self.langevin_damp} {33+self.langevin_seed} zero yes")
{self.fix_COM}postequilibration_commands.append(f"fix COMFIX vertices recenter 0 0 0")

# -------------------------#
#         RUN              #
# -------------------------#

# RUN THE SIMULATION
trilmp.run({self.MD_simulation_steps}, integrators_defined=True, fix_symbionts_near=False, 
        postequilibration_lammps_commands = postequilibration_commands)

print("End of the simulation.")

        """

        file = open(filename, "w")
        file.write(text_to_write)
        file.close()
        print("File has been produced")

