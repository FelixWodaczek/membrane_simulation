import numpy as np
import os, sys

def write_gcmc_experiment_regions_and_parameters(f, parameters):

    """
    This function appends to the python file the required elements to
    conduct GCMC simulations in different geometries

    To see in detail what each experiment means --> GradientExperimentTable Keynote

    EXP 1:
        Chemostated regions: 2
            - One with metabolites
            - One without metabolites
        Geometry:
            - Region w/o metabolites constrains region with metabolites

    EXP 2:
        Chemostated regions: 2
            - One with metabolites
            - One without metabolites
        Geometry:
            - Each chemostated region is at a extreme of the simulation box

    EXP 3:
        Chemostrated regions: 1
            - One with metabolites
        Geometry:
            - Single chemostated region at a distance from the membrane

    """

    dictionary_parameters_experiment = parameters['dictionary_parameters_experiment']
    if parameters['experiment_number']==1 or parameters['experiment_number']==2:

        f.writelines("# region where particles are added\n")
        f.writelines(f"N_gcmc_1 ={dictionary_parameters_experiment['N_gcmc_1']}\n")
        f.writelines(f"X_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"seed_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"mu_gcmc_1={dictionary_parameters_experiment['mu_gcmc_1']}\n")
        f.writelines(f"max_gcmc_1={dictionary_parameters_experiment['max_gcmc_1']}\n")
        f.writelines("\n")

        f.writelines("# region where particles are removed\n")
        f.writelines(f"N_gcmc_2 ={dictionary_parameters_experiment['N_gcmc_2']}\n")
        f.writelines(f"X_gcmc_2 ={dictionary_parameters_experiment['X_gcmc_2']}\n")
        f.writelines(f"seed_gcmc_2 ={dictionary_parameters_experiment['X_gcmc_2']}\n")
        f.writelines(f"mu_gcmc_2={dictionary_parameters_experiment['mu_gcmc_2']}\n")
        f.writelines(f"max_gcmc_2={dictionary_parameters_experiment['max_gcmc_2']}\n")
        f.writelines("\n")

        f.writelines(f"sigma_metabolites                = {dictionary_parameters_experiment['sigma_metabolites']}\n")
        f.writelines(f"interaction_range_metabolites    = {dictionary_parameters_experiment['interaction_range_metabolites']}\n")
        f.writelines(f"interaction_strength_metabolites = {dictionary_parameters_experiment['interaction_strength_metabolites']}\n")
        f.writelines("\n")

        f.writelines("# define the regions\n")
        f.writelines(f"region_1_geo=({dictionary_parameters_experiment['xlo_1']}, {dictionary_parameters_experiment['xhi_1']}, {dictionary_parameters_experiment['ylo_1']}, {dictionary_parameters_experiment['yhi_1']}, {dictionary_parameters_experiment['zlo_1']}, {dictionary_parameters_experiment['zhi_1']})\n")
        f.writelines(f"region_2_geo=({dictionary_parameters_experiment['xlo_2']}, {dictionary_parameters_experiment['xhi_2']}, {dictionary_parameters_experiment['ylo_2']}, {dictionary_parameters_experiment['yhi_2']}, {dictionary_parameters_experiment['zlo_2']}, {dictionary_parameters_experiment['zhi_2']})\n")
        f.writelines("parameters_region_1=(N_gcmc_1, X_gcmc_1, seed_gcmc_1, 1.0, mu_gcmc_1, max_gcmc_1)\n")
        f.writelines("parameters_region_2=(N_gcmc_2, X_gcmc_2, seed_gcmc_2, 1.0, mu_gcmc_2, max_gcmc_2)\n")

    if parameters['experiment_number']==3:

        f.writelines("# region where particles are added\n")
        f.writelines(f"N_gcmc_1 ={dictionary_parameters_experiment['N_gcmc_1']}\n")
        f.writelines(f"X_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"seed_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"mu_gcmc_1={dictionary_parameters_experiment['mu_gcmc_1']}\n")
        f.writelines(f"max_gcmc_1={dictionary_parameters_experiment['max_gcmc_1']}\n")
        f.writelines("\n")

        f.writelines(f"sigma_metabolites                = {dictionary_parameters_experiment['sigma_metabolites']}\n")
        f.writelines(f"interaction_range_metabolites    = {dictionary_parameters_experiment['interaction_range_metabolites']}\n")
        f.writelines(f"interaction_strength_metabolites = {dictionary_parameters_experiment['interaction_strength_metabolites']}\n")
        f.writelines("\n")

        f.writelines("# define the region\n")
        f.writelines(f"region_1_geo=({dictionary_parameters_experiment['xlo_1']}, {dictionary_parameters_experiment['xhi_1']}, {dictionary_parameters_experiment['ylo_1']}, {dictionary_parameters_experiment['yhi_1']}, {dictionary_parameters_experiment['zlo_1']}, {dictionary_parameters_experiment['zhi_1']})\n")
        f.writelines("parameters_region_1=(N_gcmc_1, X_gcmc_1, seed_gcmc_1, 1.0, mu_gcmc_1, max_gcmc_1)\n")

    if parameters['experiment_number']==4:
        f.writelines("# region where particles are added\n")
        f.writelines(f"N_gcmc_1 ={dictionary_parameters_experiment['N_gcmc_1']}\n")
        f.writelines(f"X_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"seed_gcmc_1 ={dictionary_parameters_experiment['X_gcmc_1']}\n")
        f.writelines(f"mu_gcmc_1={dictionary_parameters_experiment['mu_gcmc_1']}\n")
        f.writelines(f"max_gcmc_1={dictionary_parameters_experiment['max_gcmc_1']}\n")
        f.writelines("\n")

        f.writelines(f"sigma_metabolites                = {dictionary_parameters_experiment['sigma_metabolites']}\n")
        f.writelines(f"interaction_range_metabolites    = {dictionary_parameters_experiment['interaction_range_metabolites']}\n")
        f.writelines(f"interaction_strength_metabolites = {dictionary_parameters_experiment['interaction_strength_metabolites']}\n")
        f.writelines("\n")

        f.writelines("# define the region\n")
        f.writelines(f"region_1_geo=({dictionary_parameters_experiment['xlo_1']}, {dictionary_parameters_experiment['xhi_1']}, {dictionary_parameters_experiment['ylo_1']}, {dictionary_parameters_experiment['yhi_1']}, {dictionary_parameters_experiment['zlo_1']}, {dictionary_parameters_experiment['zhi_1']})\n")
        f.writelines("parameters_region_1=(N_gcmc_1, X_gcmc_1, seed_gcmc_1, 1.0, mu_gcmc_1, max_gcmc_1)\n")

def write_gcmc_experiment_fixes(f, parameters):

    if parameters['experiment_number']==1:
        f.writelines("              fix_gcmc=True,\n")
        f.writelines("              fix_gcmc_num_regions=2,\n")
        f.writelines("              fix_gcmc_region_type=('block','intersection'),\n")
        f.writelines("              fix_gcmc_region_parameters=(region_1_geo, (region_2_geo, region_1_geo)),\n")
        f.writelines("              fix_gcmc_fix_parameters=(parameters_region_1, parameters_region_2),\n")
        f.writelines("              fix_gcmc_interaction_parameters=(sigma_metabolites, interaction_range_metabolites, interaction_strength_metabolites),\n")
        f.writelines(" ")

    if parameters['experiment_number']==2:
        f.writelines("              fix_gcmc=True,\n")
        f.writelines("              fix_gcmc_num_regions=2,\n")
        f.writelines("              fix_gcmc_region_type=('block','block'),\n")
        f.writelines("              fix_gcmc_region_parameters=(region_1_geo, region_2_geo),\n")
        f.writelines("              fix_gcmc_fix_parameters=(parameters_region_1, parameters_region_2),\n")
        f.writelines("              fix_gcmc_interaction_parameters=(sigma_metabolites, interaction_range_metabolites, interaction_strength_metabolites),\n")
        f.writelines(" ")

    if parameters['experiment_number']==3:
        f.writelines("              fix_gcmc=True,\n")
        f.writelines("              fix_gcmc_num_regions=1,\n")
        f.writelines("              fix_gcmc_region_type=('block'),\n")
        f.writelines("              fix_gcmc_region_parameters=(region_1_geo),\n")
        f.writelines("              fix_gcmc_fix_parameters=(parameters_region_1),\n")
        f.writelines("              fix_gcmc_interaction_parameters=(sigma_metabolites, interaction_range_metabolites, interaction_strength_metabolites),\n")
        f.writelines(" ")

    if parameters['experiment_number']==4:
        f.writelines("              fix_gcmc=True,\n")
        f.writelines("              fix_gcmc_num_regions=1,\n")
        f.writelines("              fix_gcmc_region_type=('block'),\n")
        f.writelines("              fix_gcmc_region_parameters=(region_1_geo),\n")
        f.writelines("              fix_gcmc_fix_parameters=(parameters_region_1),\n")
        f.writelines("              fix_gcmc_interaction_parameters=(sigma_metabolites, interaction_range_metabolites, interaction_strength_metabolites),\n")
        f.writelines(" ")

def write_trilmp_file(filename, parameters):

    f = open(filename, 'w')

    f.writelines("import trimesh, os, pickle\n")
    f.writelines("import numpy as np\n")
    f.writelines("import pandas as pd\n")
    f.writelines("from trimem.mc.trilmp import TriLmp\n")
    f.writelines("from trimem.core import TriMesh\n")
    f.writelines("\n")

    f.writelines("# create directory to save data of interest\n")
    f.writelines("os.system('mkdir -p data')\n")
    f.writelines("\n")

    f.writelines("# initialization of the membrane mesh\n")
    f.writelines(f"sigma={parameters['sigma']}\n")
    f.writelines("mesh_coordinates = pd.read_csv('mesh_coordinates_N_5072_.dat', header = None, index_col = False, sep = ' ')\n")
    f.writelines("mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()\n")
    f.writelines("mesh_faces = pd.read_csv('mesh_faces_N_5072_.dat', header = None, index_col = False, sep = ' ')\n")
    f.writelines("mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()\n")
    f.writelines("# turning it into a trimesh object\n")
    f.writelines("mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array)\n")
    f.writelines("# rescaling it so that we start from the right distances\n")
    f.writelines("desired_average_distance = 2**(1.0/6.0) * sigma\n")
    f.writelines("current_average_distance = np.mean(mesh.edges_unique_length)\n")
    f.writelines("scaling = desired_average_distance/current_average_distance\n")
    f.writelines("mesh.vertices *=scaling\n")
    f.writelines("\n")

    f.writelines("# mechanical properties of the membrane\n")
    f.writelines(f"kappa_b = {parameters['kappa_b']}\n")
    f.writelines(f"kappa_a = {parameters['kappa_a']}\n")
    f.writelines(f"kappa_v = {parameters['kappa_v']}\n")
    f.writelines(f"kappa_c = {parameters['kappa_c']}\n")
    f.writelines(f"kappa_t = {parameters['kappa_t']}\n")
    f.writelines(f"kappa_r = {parameters['kappa_r']}\n")
    f.writelines("\n")

    f.writelines("# MD properties\n")
    f.writelines(f"step_size     = {parameters['step_size']}\n")
    f.writelines(f"traj_steps    = {parameters['traj_steps']}\n")
    f.writelines(f"langevin_damp = {parameters['langevin_damp']}\n")
    f.writelines(f"langevin_seed = {parameters['langevin_seed']}\n")
    f.writelines(f"total_sim_time={parameters['total_sim_time']}\n")
    f.writelines("total_number_steps=int(total_sim_time/(step_size*traj_steps))\n")
    f.writelines("\n")

    f.writelines("# MC/TRIMEM bond flipping properties\n")
    f.writelines(f"flip_ratio={parameters['flipratio']}\n")
    f.writelines(f"switch_mode='{parameters['switch_mode']}'\n")
    f.writelines("\n")

    f.writelines("# ouput and printing\n")
    f.writelines(f"discret_snapshots={parameters['discret_snapshots']}\n")
    f.writelines("print_frequency = int(discret_snapshots/(step_size*traj_steps))\n")
    f.writelines("\n")

    f.writelines("# simulation box parameters\n")
    f.writelines(f"xlo={parameters['xlo']}\n")
    f.writelines(f"xhi={parameters['xhi']}\n")
    f.writelines(f"ylo={parameters['xlo']}\n")
    f.writelines(f"yhi={parameters['yhi']}\n")
    f.writelines(f"zlo={parameters['zlo']}\n")
    f.writelines(f"zhi={parameters['zhi']}\n")
    f.writelines("\n")

    # include fix gcmc regions
    write_gcmc_experiment_regions_and_parameters(f, parameters)

    f.writelines("# initialization of the trilmp object -- writing out all values for initialization\n")
    f.writelines("trilmp=TriLmp(initialize=True,                          # use mesh to initialize mesh reference\n")
    f.writelines("              mesh_points=mesh.vertices,                # input mesh vertices\n")
    f.writelines("              mesh_faces=mesh.faces,                    # input of the mesh faces\n")

    f.writelines("              kappa_b=kappa_b,                          # MEMBRANE MECHANICS: bending modulus (kB T)\n")
    f.writelines("              kappa_a=kappa_a,                          # MEMBRANE MECHANICS: constraint on area change from target value (kB T)\n")
    f.writelines("              kappa_v=kappa_v,                          # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)\n")
    f.writelines("              kappa_c=kappa_c,                          # MEMBRANE MECHANICS: constraint on area difference change (understand meaning) (kB T)\n")
    f.writelines("              kappa_t=kappa_t,                          # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)\n")
    f.writelines("              kappa_r=kappa_r,                          # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)\n")

    f.writelines("              step_size=step_size,                      # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation\n")
    f.writelines("              traj_steps=traj_steps,                    # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping\n")
    f.writelines("              flip_ratio=flip_ratio,                    # MC PART SIMULATION: fraction of edges to flip?\n")
    f.writelines("              initial_temperature=1.0,                  # MD PART SIMULATION: temperature of the system\n")
    f.writelines("              langevin_damp=langevin_damp,              # MD PART SIMULATION: damping of the Langevin thermostat (as in LAMMPS)\n")
    f.writelines("              langevin_seed=langevin_seed,              # MD PART SIMULATION: seed for langevin dynamics\n")
    f.writelines("              pure_MD=True,                             # MD PART SIMULATION: accept every MD trajectory?\n")
    f.writelines("              switch_mode=switch_mode,                  # MD/MC PART SIMULATION: 'random' or 'alternating' flip-or-move\n")
    f.writelines("              box=(xlo,xhi,ylo,yhi,zlo, zhi),           # MD PART SIMULATION: simulation box properties, periodic\n")

    f.writelines("              info=print_frequency,                     # OUTPUT: frequency output in shell\n")
    f.writelines("              thin=print_frequency,                     # OUTPUT: frequency trajectory output\n")
    f.writelines("              output_prefix='data/data',                # OUTPUT: prefix for output filenames\n")
    f.writelines("              restart_prefix='data/data',               # OUTPUT: name for checkpoint files\n")
    f.writelines("              checkpoint_every=print_frequency,         # OUTPUT: interval of checkpoints (alternating pickles)\n")
    f.writelines("              output_format='lammps_txt',               # OUTPUT: choose different formats for 'lammps_txt', 'lammps_txt_folder' or 'h5_custom'\n")
    f.writelines("              output_counter=0,                         # OUTPUT: initialize trajectory number in writer class\n")
    f.writelines("              performance_increment=print_frequency,    # OUTPUT: output performace stats to prefix_performance.dat file\n")
    f.writelines("              energy_increment=print_frequency,         # OUTPUT: output energies to energies.dat file\n")

    f.writelines("              check_neigh_every=1,                      # NEIGHBOUR LISTS\n")
    f.writelines("              equilibration_rounds=1,                   # MEMBRANE EQUILIBRATION ROUNDS\n")

    write_gcmc_experiment_fixes(f, parameters)

    f.writelines("              )\n")


    f.writelines("\n")
    f.writelines("# program run\n")
    f.writelines(f"trilmp.run(total_number_steps, equilibration_gcmc={parameters['equilibration_gcmc']})\n")
    f.close()

def create_directory_name(dataset_name, names_parameters_to_vary, values_parameters_to_vary):

    dataset_name += "/data_"

    for i in range(len(names_parameters_to_vary)):
        dataset_name+=names_parameters_to_vary[i]+"_"+str(values_parameters_to_vary[i])+"_"

    # create the internal directory where we vary parameters
    isExist = os.path.exists(dataset_name)
    if not isExist:
        os.makedirs(dataset_name)

    return dataset_name

def create_parameter_log(dirname, parameters):

    # open file
    f = open(dirname+"/parameter_log.dat", "w")

    # write down parameters in log file
    for key, value in parameters.items():
        if key!='dictionary_parameters_experiment':
            f.writelines("{} {}\n".format(key, value))

    # write down parameters characterizing the experiment
    dictionary_parameters_experiment = parameters['dictionary_parameters_experiment']
    for key, value in dictionary_parameters_experiment.items():
        f.writelines("{} {}\n".format(key, value))

def copy_necessary_for_file_to_run(dirname):

    os.system(f'cp mesh_* {dirname}')
