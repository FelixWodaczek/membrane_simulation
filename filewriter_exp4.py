import os
import json
import shutil
from pathlib import Path
from itertools import product

import filewriter_utils

def main():
    # general names
    dataset_name = Path("240216_experiment4_varying_intrange_mu_nevery").resolve()

    # create header directory
    isExist = dataset_name.exists()
    if not isExist:
        dataset_name.mkdir()
    
    shutil.copy("mesh_coordinates_N_5072_.dat", dataset_name.joinpath("mesh_coordinates_N_5072_.dat"))
    shutil.copy("mesh_faces_N_5072_.dat", dataset_name.joinpath("mesh_faces_N_5072_.dat"))
    # metadata collection for sbatch
    fdir = open(dataset_name.joinpath("directories_in_directory.dat"), "w")

    # mechanical properties membrane
    sigma_membrane=1.0
    kappa_b=20.0
    kappa_a=2.5e5
    kappa_v=2.5e5
    kappa_c=0.0
    kappa_t=1.0e4
    kappa_r=1.0e4

    # trimem parameters
    traj_steps=50
    flip_ratio=0.1

    # MD parameters
    xlo, xhi = -50, 50
    ylo, yhi = -50, 50
    zlo, zhi = -50, 50
    step_size=0.001
    langevin_damp=1.0
    langevin_seed=123
    switch_mode="random"
    total_sim_time=100000  # in time units
    discrete_snapshots=10   # in time units

    # Parameters for the GCMC
    equilibration_gcmc=0 # in sim step

    rcs_membrane_metabolite = [1.5, 2.5]            # 2 handled
    interaction_strengths = [10.0] # [1.0, 5.0, 10.0, 20.0]  # 4 handled
    geometric_factors = [0.2] # [1.0, 0.5, 0.2]             # 3 handled
    prob_transforms = [1.0] # [1.0, 0.1, 0.0]               # 3 handled
    variable_factors = [10] #[1, 5, 10]                   # 3 handled
    vfracs = [0.01] # [0.05, 0.01]
    mu_gcmc_1s = [0, -5, -10]
    Neverys = [100, 50, 10]

    for rc_mm, interaction_strength, geometric_factor, prob_transform, variable_factor, vfrac, mu_gcmc_1, Nevery in product(
        rcs_membrane_metabolite, interaction_strengths,
        geometric_factors, prob_transforms, variable_factors, 
        vfracs, mu_gcmc_1s, Neverys
    ):
        N_gcmc_1=int(variable_factor*langevin_damp/step_size)
        X_gcmc_1=100
        seed_gcmc_1 = langevin_seed
        # mu_gcmc_1=0
        # xlo_1 = -10
        # xhi_1 = 10
        # ylo_1 = -10
        # yhi_1 = 10
        # zlo_1 = -10
        # zhi_1 = 10

        sigma_metabolites = 1.0

        dirname = "data_trajsteps_"+str(traj_steps)
        dirname +="_flipratio_"+str(flip_ratio)
        dirname +="_stepsize_"+str(step_size)
        dirname +="_seed_"+str(langevin_seed)
        dirname +="_kb_"+str(kappa_b)+"_kv_"+str(kappa_v)
        dirname +="_ka_"+str(kappa_a)+"_kc_"+str(kappa_c)
        dirname +="_kt_"+str(kappa_t)
        dirname +="_kr_"+str(kappa_r)
        dirname +="_lgvdamp_"+str(langevin_damp)
        dirname +="_exp_4"
        dirname +="_eqgcmc_"+str(equilibration_gcmc)
        dirname +="_vfrac_"+str(vfrac)
        dirname +="_mu_"+str(mu_gcmc_1)
        dirname +="_nevery_"+str(Nevery)
        dirname +="_intrnge_"+str(rc_mm)
        dirname +="_intstrgth_"+str(interaction_strength)
        dirname +="_geomfact_"+str(geometric_factor)
        dirname +="_probtrans_"+str(prob_transform)
        dirname +="_factors_"+str(variable_factor)

        dirname = dataset_name.joinpath(dirname)
        fdir.writelines(dirname.name+"\n")

        # create the internal directory where we vary parameters
        isExist = dirname.exists()
        if not isExist:
            dirname.mkdir()
        shutil.copy("trilmp_srp_pot.py", dirname.joinpath("trilmp_srp_pot.py"))

        # experiment dictionaries
        experiment_dictionary={
            'N_gcmc_1': N_gcmc_1,
            'X_gcmc_1': X_gcmc_1,
            'seed_gcmc_1': seed_gcmc_1,
            'mu_gcmc_1': mu_gcmc_1,
            'Nevery': Nevery,
            'vfrac': vfrac,
            # 'xlo_1': xlo_1,
            # 'xhi_1': xhi_1,
            'prob_transform': prob_transform,
            'geometric_factor': geometric_factor,
            'sigma_metabolites': sigma_metabolites,
            'interaction_range_metabolites': rc_mm,
            'interaction_strength_metabolites': interaction_strength
        }
        dictionary_parameters_experiment=experiment_dictionary


        parameters = {
            'traj_steps': traj_steps,
            'flip_ratio': flip_ratio,
            'sigma_membrane': sigma_membrane,
            'kappa_b': kappa_b,
            'kappa_a': kappa_a,
            'kappa_c': kappa_c,
            'kappa_r': kappa_r,
            'kappa_t': kappa_t,
            'kappa_v': kappa_v,
            'step_size': step_size,
            'langevin_damp': langevin_damp,
            'langevin_seed': langevin_seed,
            'switch_mode': switch_mode,
            'total_sim_time': total_sim_time,
            'discrete_snapshots': discrete_snapshots,
            'xlo': xlo,
            'xhi': xhi,
            'ylo': ylo,
            'yhi': yhi,
            'zlo': zlo,
            'zhi': zhi,
            'equilibration_gcmc': equilibration_gcmc,
            'experiment_number': 4,
            'dictionary_parameters_experiment': dictionary_parameters_experiment
        }

        with open(dirname.joinpath("parameter_dict.json"), "w") as f:
            json.dump(parameters, f, indent=4)
            f.close()

if __name__ == "__main__":
    main()