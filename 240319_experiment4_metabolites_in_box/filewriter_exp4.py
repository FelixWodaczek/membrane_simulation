import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
from itertools import product

import py_src.simulation_utils as sutils

def write_files(target_dir_name: Path):
    # general names
    dataset_name = target_dir_name
    
    # metadata collection for sbatch
    fdir = open(dataset_name.joinpath("directories_in_directory.dat"), "w")

    param_dict = {}

    # mechanical properties membrane
    membrane_parameters = sutils.MembraneParameters(
        sigma_vertex = 1.0,
        vertex_mass = 1.0,
        kappa_b = 20.0,
        kappa_a = 2.5e5,
        kappa_v = 2.5e5,
        kappa_c = 0.0,
        kappa_t = 1.0e4,
        kappa_r = 1.0e4
    )
    param_dict["membrane_parameters"] = asdict(membrane_parameters)

    # trimem parameters
    trimem_parameters = sutils.TrimemParameters(
        traj_steps = 50,
        total_sim_time = 100000, # in time units
        step_size = 0.001,
        discrete_snapshots = 10,   # in time units
        flip_ratio = 0.1,
        initial_temperature=1.0, # MD PART SIMULATION: temperature of the system
        pure_MD=True # MD PART SIMULATION: accept every MD trajectory?
    )
    param_dict["trimem_parameters"] = asdict(trimem_parameters)

    # Parameters for the GCMC
    equilibration_gcmc=0 # in sim step

    rcs_membrane_metabolite = [1.5, 2.5]
    interaction_strengths = [10.0] # [1.0, 5.0, 10.0, 20.0]
    geometric_factors = [0.6, 0.4, 0.2] # [1.0, 0.5, 0.2]
    prob_transforms = [1.0, 0.1, 0.01, 0.0]
    variable_factors = [10] #[1, 5, 10]
    vfracs = [0.05]

    for rc_mm, interaction_strength, geometric_factor, prob_transform, variable_factor, vfrac in product(
        rcs_membrane_metabolite, interaction_strengths,
        geometric_factors, prob_transforms, variable_factors, vfracs
    ):  
        gcmc_parameters = sutils.GCMCParameters(
            langevin_damp = 1.,
            X = 100,
            seed = 123,
            mu = 0,
            geometric_factor = geometric_factor,
            vfrac = vfrac,
            variable_factor = variable_factor
        )
        param_dict['gcmc_parameters'] = asdict(gcmc_parameters)

        interaction_parameters = sutils.InteractionParameters(
            interaction_range = rc_mm,
            interaction_strength = interaction_strength,
            sigma_metabolites = 1.0
        )
        param_dict['interaction_parameters'] = asdict(interaction_parameters)

        chemistry_group = sutils.DynamicGroup(
            interaction_range=interaction_parameters.interaction_range_tilde(membrane_parameters.sigma_vertex),
            name='waste', probability=prob_transform, target_type=2, check_group_every=100,
            seed=gcmc_parameters.seed
        )
        chemistry_gcmc = sutils.GCMC(
            name='waste', target_group='waste', N=100, X=100, seed=gcmc_parameters.seed, mu=-100, maxp=1
        )
        param_dict['chemistry_group'] = asdict(chemistry_group)
        param_dict['chemistry_gcmc'] = asdict(chemistry_gcmc)

        dirname = "data_trajsteps_"+str(trimem_parameters.traj_steps)
        dirname +="_flipratio_"+str(trimem_parameters.flip_ratio)
        dirname +="_stepsize_"+str(trimem_parameters.step_size)
        dirname +="_seed_"+str(gcmc_parameters.seed)
        dirname +="_kb_"+str(membrane_parameters.kappa_b)+"_kv_"+str(membrane_parameters.kappa_v)
        dirname +="_ka_"+str(membrane_parameters.kappa_a)+"_kc_"+str(membrane_parameters.kappa_c)
        dirname +="_kt_"+str(membrane_parameters.kappa_t)
        dirname +="_kr_"+str(membrane_parameters.kappa_r)
        dirname +="_lgvdamp_"+str(gcmc_parameters.langevin_damp)
        dirname +="_exp_4"
        dirname +="_eqgcmc_"+str(equilibration_gcmc)
        dirname +="_vfrac_"+str(vfrac)
        dirname +="_intrnge_"+str(rc_mm)
        dirname +="_intstrgth_"+str(interaction_strength)
        dirname +="_geomfact_"+str(geometric_factor)
        dirname +="_probtrans_"+str(prob_transform)
        dirname +="_factors_"+str(variable_factor)

        dirname = dataset_name.joinpath(dirname)
        fdir.writelines(dirname.name+"\n")

        # create the internal directory where we vary parameters
        if not dirname.exists():
            dirname.mkdir()
        if not dirname.joinpath('data').exists():
            dirname.joinpath('data').mkdir()
        # shutil.copy("trilmp_srp_pot.py", dirname.joinpath("trilmp_srp_pot.py"))

        with open(dirname.joinpath("parameter_dict.json"), "w") as f:
            json.dump(param_dict, f, indent=4)
            f.close()

def main():
    parser = argparse.ArgumentParser(description="Write experiment files")
    parser.add_argument(
        "-t", "--target_directory", 
        type=str,
        help="Name of directory to be created"
    )
    parser.add_argument(
        '-o', '--overwrite',
        action='store_true',
        help='Overwrite directory tree if it already exists'
    )
    args = parser.parse_args()

    if args.target_directory is None:
        print("No target name given. Supply a directory name using '-t <directory_name>'")
        return

    overwrite = bool(args.overwrite)
    target_dir = Path(args.target_directory).resolve()

    # create header directory
    if target_dir.exists() and target_dir.is_dir():
        if not overwrite:
            print(f"Target directory {target_dir} already exists. To overwrite it supply '-o'.")
            return
        shutil.rmtree(target_dir)

    target_dir.mkdir()
    target_dir.joinpath('logs').mkdir()
    write_files(target_dir_name=target_dir)

if __name__ == "__main__":
    main()