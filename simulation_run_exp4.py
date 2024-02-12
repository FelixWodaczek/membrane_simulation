from pathlib import Path
import json
import os

import numpy as np
import pandas as pd

import trimesh
from trimem.mc.trilmp import TriLmp
from trimem.core import TriMesh

MOLECULE_TEMPLATES_PATH = Path(__file__).resolve().parent.joinpath("MoleculeTemplates_LAMMPSFixBondReact", "ReactionsExperiment4").absolute()
PREDEGRADE_PATH = MOLECULE_TEMPLATES_PATH.joinpath("pre_DegradeWaste.txt")
POSTDEGRADE_PATH = MOLECULE_TEMPLATES_PATH.joinpath("post_DegradeWaste.txt")
MAPDEGRADE_PATH = MOLECULE_TEMPLATES_PATH.joinpath("map_DegradeWaste.txt")
PRETRANSFORM_PATH = MOLECULE_TEMPLATES_PATH.joinpath("pre_TransformMetabolites.txt")
POSTTRANSFORM_PATH = MOLECULE_TEMPLATES_PATH.joinpath("post_TransformMetabolites.txt")
MAPTRANSFORM_PATH = MOLECULE_TEMPLATES_PATH.joinpath("map_TransformMetabolites.txt")

class SimulationManager:
    def __init__(self, target_dir: Path):
        self.target_dir = target_dir.resolve().absolute()
        self.parent_dir = target_dir.parent

        os.chdir(str(target_dir)) # Not very clean, this should be avoided technically
    
        with open(target_dir.joinpath("parameter_dict.json")) as f:
            self.param_dict = json.load(f)
            f.close()

        # create directory to save data of interest
        target_dir.joinpath("data").mkdir(exist_ok=True)

    @staticmethod
    def add_walls_to_lammps(trilmp, ylo, yhi, zlo, zhi):
        trilmp.lmp.commands_string(f"fix ConfinementMet metabolites wall/reflect xhi EDGE ylo {ylo} yhi {yhi} zlo {zlo} zhi {zhi}")

    
    @staticmethod
    def add_chemistry_to_lammps(
        trilmp, prob_transform, N_gcmc_1: int,
        sigma_tilde_membrane_metabolite, rc_tilde_membrane_metabolite,
        prob_degrade: float=1.0, Nevery=100, seedR=2430
    ):
        chemistry_commands = f"""
# FIX BOND/REACT SECTION

# Transform metabolites to waste
molecule mpreTransform {PRETRANSFORM_PATH}
molecule mpostTransform {POSTTRANSFORM_PATH}

# Degrade waste
molecule mpreDegradeWaste {PREDEGRADE_PATH}
molecule mpostDegradeWaste {POSTDEGRADE_PATH}

fix freact all bond/react reset_mol_ids no &
  react DegradeWaste all {Nevery} {sigma_tilde_membrane_metabolite} {rc_tilde_membrane_metabolite} mpreDegradeWaste mpostDegradeWaste {MAPDEGRADE_PATH} prob {prob_degrade} {seedR} &
  react TransformWaste all {Nevery+1} {sigma_tilde_membrane_metabolite} {rc_tilde_membrane_metabolite} mpreTransform mpostTransform {MAPTRANSFORM_PATH} prob {prob_transform} {seedR+19}

fix  aveMET all ave/time {N_gcmc_1} 1 {N_gcmc_1} v_n2 c_TempCompute c_pe file 'metabolite_properties.dat'
fix  aveREC all ave/time {N_gcmc_1} 1 {N_gcmc_1} f_freact file 'reactions.dat' mode vector
"""
        trilmp.lmp.commands_string(chemistry_commands)

    def create_parameter_log(self, dirname):
        # open file
        f = open(dirname+"/parameter_log.dat", "w")

        # write down parameters in log file
        for key, value in self.param_dict.items():
            if key!='dictionary_parameters_experiment':
                f.writelines("{} {}\n".format(key, value))

        # write down parameters characterizing the experiment
        dictionary_parameters_experiment = self.param_dict['dictionary_parameters_experiment']
        for key, value in dictionary_parameters_experiment.items():
            f.writelines("{} {}\n".format(key, value))

    def run(self):
        sigma_membrane = self.param_dict['sigma_membrane']
        # initialization of the membrane mesh
        mesh_coordinates = pd.read_csv(self.parent_dir.joinpath('mesh_coordinates_N_5072_.dat'), header = None, index_col = False, sep = ' ')
        mesh_coordinates_array = mesh_coordinates[[1, 2, 3]].to_numpy()
        mesh_faces = pd.read_csv(self.parent_dir.joinpath('mesh_faces_N_5072_.dat'), header = None, index_col = False, sep = ' ')
        mesh_faces_array = mesh_faces[[0, 1, 2]].to_numpy()
        # turning it into a trimesh object
        mesh = trimesh.Trimesh(vertices=mesh_coordinates_array, faces = mesh_faces_array)
        # rescaling it so that we start from the right distances
        desired_average_distance = 2**(1.0/6.0) * sigma_membrane
        current_average_distance = np.mean(mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        mesh.vertices *= scaling

        # mechanical properties of the membrane
        kappa_b = self.param_dict['kappa_b']
        kappa_a = self.param_dict['kappa_a']
        kappa_v = self.param_dict['kappa_v']
        kappa_c = self.param_dict['kappa_c']
        kappa_t = self.param_dict['kappa_t']
        kappa_r = self.param_dict['kappa_r']

        # MD properties
        step_size = self.param_dict['step_size']
        traj_steps = self.param_dict['traj_steps']
        langevin_damp = self.param_dict['langevin_damp']
        langevin_seed = self.param_dict['langevin_seed']
        total_sim_time = self.param_dict['total_sim_time']
        total_number_steps = int(total_sim_time/(step_size*traj_steps))

        # MC/TRIMEM bond flipping properties
        flip_ratio = self.param_dict['flip_ratio']
        switch_mode = self.param_dict['switch_mode']

        # ouput and printing
        discrete_snapshots=self.param_dict['discrete_snapshots']
        print_frequency = int(discrete_snapshots/(step_size*traj_steps))

        # simulation box parameters
        xlo = self.param_dict['xlo']
        xhi = self.param_dict['xhi']
        ylo = self.param_dict['xlo']
        yhi = self.param_dict['yhi']
        zlo = self.param_dict['zlo']
        zhi = self.param_dict['zhi']

        experiment_param_dict: dict = self.param_dict['dictionary_parameters_experiment']

        # region where particles are added
        N_gcmc_1 = experiment_param_dict["N_gcmc_1"]
        X_gcmc_1 = experiment_param_dict["X_gcmc_1"]
        seed_gcmc_1 = experiment_param_dict["seed_gcmc_1"]
        mu_gcmc_1 = experiment_param_dict["mu_gcmc_1"]
        max_gcmc_1 = experiment_param_dict["max_gcmc_1"]

        # interaction parameters
        sigma_metabolites = experiment_param_dict["sigma_metabolites"]
        interaction_range_metabolites = experiment_param_dict["interaction_range_metabolites"]
        interaction_strength_metabolites = experiment_param_dict["interaction_strength_metabolites"]

        sigma_tilde_membrane_metabolite = 0.5*(sigma_metabolites+sigma_membrane)
        rc_tilde_membrane_metabolite = interaction_range_metabolites*sigma_tilde_membrane_metabolite

        # define the region
        x_membrane_max = np.max(mesh_coordinates_array[:, 0])
        r_mean = np.mean(
            np.linalg.norm(
                np.mean(mesh_coordinates_array, axis=0)-mesh_coordinates_array, axis=1
            )
        )
        height_width = r_mean*experiment_param_dict['geometric_factor']
        region_1_geo = (
            x_membrane_max, self.param_dict['xhi'], # xlo, xhi
            -height_width/2, height_width/2, # ylo, yhi
            -height_width/2, height_width/2, # zlo, zhi
        )
        parameters_region_1 = (N_gcmc_1, X_gcmc_1, seed_gcmc_1, 1.0, mu_gcmc_1, max_gcmc_1)
        
        # chemistry
        Nevery = 100
        seedR = 2430
        chemistry_parameters = [
            # DegradeWaste
            [ 
                str(PREDEGRADE_PATH), # pre path
                str(POSTDEGRADE_PATH), # post path
                Nevery, # NeveryR
                sigma_tilde_membrane_metabolite, # RminR (minimum reaction distance?)
                rc_tilde_membrane_metabolite, # RmaxR (maximum reaction distance?)
                1.0, # reaction probability
                seedR, # seedR
                str(MAPDEGRADE_PATH) # map path? does this go here?
            ], 
            # Transform
            [
                str(PRETRANSFORM_PATH), # pre path
                str(POSTTRANSFORM_PATH), # post path
                Nevery+1, # NeveryR
                sigma_tilde_membrane_metabolite, # RminR (minimum reaction distance?)
                rc_tilde_membrane_metabolite, # RmaxR (maximum reaction distance?)
                experiment_param_dict['prob_transform'], # reaction probability
                seedR+19, # seedR
                str(MAPTRANSFORM_PATH) # map path? does this go here?
            ]
        ]

        # initialization of the trilmp object -- writing out all values for initialization
        trilmp = TriLmp(
            initialize=True,                          # use mesh to initialize mesh reference
            mesh_points=mesh.vertices,                # input mesh vertices
            mesh_faces=mesh.faces,                    # input of the mesh faces
            kappa_b=kappa_b,                          # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a=kappa_a,                          # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v=kappa_v,                          # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c=kappa_c,                          # MEMBRANE MECHANICS: constraint on area difference change (understand meaning) (kB T)
            kappa_t=kappa_t,                          # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r=kappa_r,                          # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            step_size=step_size,                      # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps=traj_steps,                    # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio=flip_ratio,                    # MC PART SIMULATION: fraction of edges to flip?
            initial_temperature=1.0,                  # MD PART SIMULATION: temperature of the system
            langevin_damp=langevin_damp,              # MD PART SIMULATION: damping of the Langevin thermostat (as in LAMMPS)
            langevin_seed=langevin_seed,              # MD PART SIMULATION: seed for langevin dynamics
            pure_MD=True,                             # MD PART SIMULATION: accept every MD trajectory?
            switch_mode=switch_mode,                  # MD/MC PART SIMULATION: 'random' or 'alternating' flip-or-move
            box=(xlo,xhi,ylo,yhi,zlo, zhi),           # MD PART SIMULATION: simulation box properties, periodic
            info=print_frequency,                     # OUTPUT: frequency output in shell
            thin=print_frequency,                     # OUTPUT: frequency trajectory output
            output_prefix='data/data',                # OUTPUT: prefix for output filenames
            restart_prefix='data/data',               # OUTPUT: name for checkpoint files
            checkpoint_every=print_frequency,         # OUTPUT: interval of checkpoints (alternating pickles)
            output_format='lammps_txt',               # OUTPUT: choose different formats for 'lammps_txt', 'lammps_txt_folder' or 'h5_custom'
            output_counter=0,                         # OUTPUT: initialize trajectory number in writer class
            performance_increment=print_frequency,    # OUTPUT: output performace stats to prefix_performance.dat file
            energy_increment=print_frequency,         # OUTPUT: output energies to energies.dat file
            check_neigh_every=1,                      # NEIGHBOUR LISTS
            equilibration_rounds=1,                   # MEMBRANE EQUILIBRATION ROUNDS
            fix_gcmc=True,
            fix_gcmc_num_regions=1,
            fix_gcmc_region_type=('block'),
            fix_gcmc_region_parameters=(region_1_geo),
            fix_gcmc_fix_parameters=(parameters_region_1),
            fix_gcmc_interaction_parameters=(sigma_metabolites, interaction_range_metabolites, interaction_strength_metabolites),
            # FIX BOND/REACT SECTION
            fix_bond_react=True,
            fix_bond_react_reactions=len(chemistry_parameters),
            fix_bond_react_reactions_parameters=chemistry_parameters,
            fix_bond_react_interactions=None, # Not used
            # fix_bond_react_flag=0, # Sets itself anyways
        )

        self.add_walls_to_lammps(trilmp, -height_width/2, height_width/2, -height_width/2, height_width/2)
        if False:
            self.add_chemistry_to_lammps(
                trilmp, prob_transform=experiment_param_dict['prob_transform'],
                N_gcmc_1=experiment_param_dict['N_gcmc_1'], 
                sigma_tilde_membrane_metabolite=sigma_tilde_membrane_metabolite,
                rc_tilde_membrane_metabolite=rc_tilde_membrane_metabolite
            )

        trilmp.run(total_number_steps, equilibration_gcmc=self.param_dict['equilibration_gcmc'])

def main():
    target_dir = Path("experiment4_varying_factors_distance_boxheight").resolve()
    for target_dir in target_dir.glob("data*/"):

        sim_manager = SimulationManager(target_dir=target_dir)
        sim_manager.run()
        return
        print("Simulation finished")

if __name__ == '__main__':
    main()