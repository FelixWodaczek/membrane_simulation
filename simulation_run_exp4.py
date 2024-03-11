import argparse
from pathlib import Path
import os
import shutil

import numpy as np

import py_src.simulation_manager as sm
import py_src.simulation_utils as sutils

def main():
    simulation_manager = sm.SimulationManager(resolution=2)
    simulation_manager.init_trilmp()

    template_path = Path(__file__).resolve().parent.joinpath('reaction_templates')

    # Define pair styles
    table_ps = sutils.TablePairStyle(
        style='linear', N=2000,
        coeff_commands = ['1 1 table trimem_srp.table trimem_srp'],
        modify_commands = ['pair table special lj/coul 0.0 0.0 0.0 tail no']
    )
    harmonic_ps = sutils.HarmonicCutPairStyle()
    harmonic_ps.set_all_repulsive_commands(2, simulation_manager.membrane_params.sigma_vertex, simulation_manager.sigma_tilde_membrane_metabolites)
    lj_ps = sutils.LJCutPairStyle(cutoff=2.5)
    lj_ps.set_membrane_attraction(2, interaction_strength=simulation_manager.interaction_strength, sigma_tilde=simulation_manager.sigma_tilde_membrane_metabolites, interaction_range=simulation_manager.interaction_range_tilde)
    
    # Add pair styles to self
    simulation_manager.pair_styles = [table_ps, harmonic_ps]

    # add a gcmc region
    variable_factor = 10

    vfrac = 0.05
    geometric_factor = 1.0

    x_membrane_max = np.max(simulation_manager.mesh.vertices[:, 0])
    r_mean = np.mean(
        np.linalg.norm(
            np.mean(simulation_manager.mesh.vertices, axis=0)-simulation_manager.mesh.vertices, axis=1
        )
    )
    height_width = r_mean*geometric_factor

    vtotal_region = (simulation_manager.box.xhi-x_membrane_max)*(height_width)*(height_width)
    maxp = int((vfrac*vtotal_region*3)/(4*np.pi*(simulation_manager.sigma_metabolites*0.5)**3))

    gcmc_region_1 = sutils.GCMCRegion(
        name = '1',
        N = int(variable_factor*simulation_manager.langevin_damp/simulation_manager.step_size),
        X = 100,
        seed = simulation_manager.langevin_seed,
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
        initial_temperature=simulation_manager.initial_temperature, langevin_damp=simulation_manager.langevin_damp, langevin_seed=simulation_manager.langevin_seed, fix_gcmc=True
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
        membrane_vertex_mass=simulation_manager.membrane_params.vertex_mass, initial_temperature=simulation_manager.initial_temperature,
        langevin_damp=simulation_manager.langevin_damp, langevin_seed=simulation_manager.langevin_seed
    ))

    # Activate attraction and reset interactions
    harmonic_ps.set_metabolite_repulsive_commands(n_types=2, sigma_metabolites=simulation_manager.sigma_metabolites)
    simulation_manager.pair_styles += [lj_ps]
    postequilibration_lammps_commands.append(simulation_manager.get_pair_style_commands())

    simulation_manager.trilmp.run(simulation_manager.total_sim_time, fix_symbionts_near=False, integrators_defined=True, postequilibration_lammps_commands=postequilibration_lammps_commands)

if __name__ == '__main__':
    # get target directory from command line using argparse and chdir there
    parser = argparse.ArgumentParser(description='Run simulation')
    parser.add_argument('-t', '--target-dir', default=None, type=str, help='Target directory')
    args = parser.parse_args()

    if args.target_dir is not None:
        target_dir = Path(args.target_dir).resolve()
        if not target_dir.is_dir():
            raise ValueError(f'Target directory {target_dir} does not exist')
        
        shutil.rmtree(target_dir)
        os.mkdir(target_dir)
        os.mkdir(target_dir.joinpath('data'))

        os.chdir(target_dir)
    main()