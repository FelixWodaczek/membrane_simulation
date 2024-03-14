import argparse
from pathlib import Path
import os
import shutil
import json

import numpy as np

import py_src.simulation_manager as sm
import py_src.simulation_utils as sutils

MESH_PATH = Path(__file__).resolve().parent

def read_parameters(target_path: Path):
    with open(target_path.joinpath('parameter_dict.json'), 'r') as f:
        param_dict = json.load(f)
        f.close()

    return (
        sutils.InteractionParameters(**param_dict['interaction_parameters']),
        sutils.MembraneParameters(**param_dict['membrane_parameters']),
        sutils.TrimemParameters(**param_dict['trimem_parameters']),
        sutils.GCMCParameters(**param_dict['gcmc_parameters']),
        sutils.DynamicGroup(**param_dict['chemistry_group']),
        sutils.GCMC(**param_dict['chemistry_gcmc'])
    )

def main():
    interaction_parameters, membrane_parameters, trimem_parameters, gcmc_parameters, chemistry_group, chemistry_gcmc = read_parameters(Path('.'))

    simulation_manager = sm.SimulationManager(
        trimem_params=trimem_parameters, 
        membrane_params=membrane_parameters, 
        mesh=MESH_PATH
    )
    simulation_manager.init_trilmp()

    # Define pair styles
    table_ps = sutils.TablePairStyle(
        style='linear', N=2000,
        coeff_commands = ['1 1 table trimem_srp.table trimem_srp'],
        modify_commands = ['pair table special lj/coul 0.0 0.0 0.0 tail no']
    )

    harmonic_ps = sutils.HarmonicCutPairStyle()
    harmonic_ps.set_all_repulsive_commands(
        2, simulation_manager.membrane_params.sigma_vertex, 
        interaction_parameters.sigma_tilde_membrane_metabolites(membrane_parameters.sigma_vertex)
    )
    lj_ps = sutils.LJCutPairStyle(cutoff=2.5)
    lj_ps.set_membrane_attraction(
        2, interaction_strength=interaction_parameters.interaction_strength, 
        sigma_tilde=interaction_parameters.sigma_tilde_membrane_metabolites(membrane_parameters.sigma_vertex), 
        interaction_range=interaction_parameters.interaction_range_tilde(membrane_parameters.sigma_vertex)
    )
    
    # Add pair styles to self
    simulation_manager.pair_styles = [table_ps, harmonic_ps]

    x_membrane_max = np.max(simulation_manager.mesh.vertices[:, 0])
    r_mean = np.mean(
        np.linalg.norm(
            np.mean(simulation_manager.mesh.vertices, axis=0)-simulation_manager.mesh.vertices, axis=1
        )
    )
    height_width = r_mean*gcmc_parameters.geometric_factor*2

    vtotal_region = (simulation_manager.box.xhi-x_membrane_max)*(height_width)*(height_width)
    maxp = int((gcmc_parameters.vfrac*vtotal_region*3)/(4*np.pi*(interaction_parameters.sigma_metabolites*0.5)**3))

    gcmc_region_1 = sutils.GCMCRegion(
        name = '1',
        N = int(gcmc_parameters.variable_factor*gcmc_parameters.langevin_damp/trimem_parameters.step_size),
        X = 100,
        seed = gcmc_parameters.seed,
        mu = 0,
        box = sutils.Box(
            x_membrane_max, simulation_manager.box.xhi,
            -height_width/2, height_width/2,
            -height_width/2, height_width/2,
        ),
        maxp = maxp
    )

    pre_equilibration_lammps_commands = []

    # Set up the chemostatted region
    pre_equilibration_lammps_commands.append(gcmc_region_1.region_command())
    pre_equilibration_lammps_commands.append(gcmc_region_1.fix_gcmc_command())

    # Add walls before equilibriation
    pre_equilibration_lammps_commands.append(simulation_manager.walls_command(
        gcmc_region_1.box.ylo, gcmc_region_1.box.yhi, gcmc_region_1.box.zlo, gcmc_region_1.box.zhi
    ))

    # Set equilibriation thermostat
    # TODO: only different because of tally yes, remove?
    pre_equilibration_lammps_commands.append(simulation_manager.equilibriation_thermostat_command(
        initial_temperature=trimem_parameters.initial_temperature, langevin_damp=gcmc_parameters.langevin_damp, langevin_seed=gcmc_parameters.seed, fix_gcmc=True
    ))

    # Set interaction parameters
    pre_equilibration_lammps_commands.append(simulation_manager.get_pair_style_commands())

    simulation_manager.trilmp.lmp.commands_string('\n'.join(pre_equilibration_lammps_commands))

    postequilibration_lammps_commands = []
    # Unfix stuff for some reason if in equilibriation
    postequilibration_lammps_commands.append(f"unfix vertexnve")
    postequilibration_lammps_commands.append(f"unfix lvt")

    # Fix thermostat again without tally yes
    postequilibration_lammps_commands.append(simulation_manager.langevin_commands(
        membrane_vertex_mass=simulation_manager.membrane_params.vertex_mass, initial_temperature=trimem_parameters.initial_temperature,
        langevin_damp=gcmc_parameters.langevin_damp, langevin_seed=gcmc_parameters.seed
    ))

    # Activate attraction and reset interactions
    harmonic_ps.set_metabolite_repulsive_commands(n_types=2, sigma_metabolites=interaction_parameters.sigma_metabolites)
    simulation_manager.pair_styles += [lj_ps]
    postequilibration_lammps_commands.append(simulation_manager.get_pair_style_commands())

    postequilibration_lammps_commands.append(chemistry_group.get_group_commands())

    postequilibration_lammps_commands.append(chemistry_gcmc.fix_gcmc_command())
    if False:
        postequilibration_lammps_commands.append('\n'.join([
            # print number of waste atoms
            "compute n_waste waste count/type atom",
            "thermo_style custom step c_th_pe c_th_ke c_n_waste[*]",
            "thermo_modify norm no"
        ]))

    simulation_manager.trilmp.run(trimem_parameters.total_sim_time, fix_symbionts_near=False, integrators_defined=True, postequilibration_lammps_commands=postequilibration_lammps_commands)

if __name__ == '__main__':
    # get target directory from command line using argparse and chdir there
    parser = argparse.ArgumentParser(description="Log File Analysis")
    parser.add_argument(
        "-t", "--target_directory", 
        type=str,
        help="Path to the target directory"
    )
    args = parser.parse_args()
    target_dir = Path(args.target_directory).resolve()

    if not target_dir.is_dir():
        raise ValueError(f'Target directory {target_dir} does not exist')

    if os.environ.get("SLURM_ARRAY_TASK_ID") is not None:
        n_dir = int(os.environ.get("SLURM_ARRAY_TASK_ID")) - 1
        with open(target_dir.joinpath("directories_in_directory.dat"), 'r') as f:
            subdirs = f.read().splitlines()
            f.close()
        
        target_dir = target_dir.joinpath(subdirs[n_dir])
    
    os.chdir(target_dir)

    main()