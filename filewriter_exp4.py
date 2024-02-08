import os
from itertools import product

import filewriter_utils

def main():
    # general names
    run_fname = "trilmp_gcmc.py"
    dataset_name = "experiment4_varying_factors_distance_boxheight"

    # create header directory
    isExist = os.path.exists(dataset_name)
    if not isExist:
        os.makedirs(dataset_name)

    # metadata collection for sbatch
    fdir = open(dataset_name+"/directories_in_directory.dat", "w")

    # mechanical properties membrane
    sigma=1.0
    kappa_b=20.0
    kappa_a=2.5e5
    kappa_v=2.5e5
    kappa_c=0.0
    kappa_t=1.0e4
    kappa_r=1.0e4

    # trimem parameters
    traj_steps=50
    flipratio=0.1

    # MD parameters
    xlo, xhi = -50, 50
    ylo, yhi = -50, 50
    zlo, zhi = -50, 50
    step_size=0.001
    langevin_damp=1.0
    langevin_seed=123
    switch_mode="random"
    total_sim_time=100000  # in time units
    discret_snapshots=10   # in time units

    # Parameters for the GCMC
    equilibration_gcmc=500000 # in sim step

    rcs_membrane_metabolite = [1.5, 2.5]            # 2
    interaction_strengths = [1.0, 5.0, 10.0, 20.0]  # 4
    geometric_factors = [1.0, 0.5, 0.2]             # 3
    prob_transforms = [1.0, 0.1, 0.0]               # 3
    variable_factors = [1, 5, 10]                   # 3

    for rc_mm, interaction_strength, geometric_factor, prob_transform, variable_factor in product(
        rcs_membrane_metabolite, interaction_strengths,
        geometric_factors, prob_transforms, variable_factors
    ):
        
        N_gcmc_1=int(variable_factor*langevin_damp/step_size)
        X_gcmc_1=100
        seed_gcmc_1 = langevin_seed
        mu_gcmc_1=0
        max_gcmc_1=1000000
        xlo_1 = -10 # TODO: this the correct place?
        xhi_1 = 10
        ylo_1 = -10
        yhi_1 = 10
        zlo_1 = -10
        zhi_1 = 10

        sigma_metabolites = 1.0
        interaction_range_metabolites = 2.5
        interaction_strength_metabolites = 10

        dirname = "data_trajsteps_"+str(traj_steps)
        dirname +="_flipratio_"+str(flipratio)
        dirname +="_stepsize_"+str(step_size)
        dirname +="_seed_"+str(langevin_seed)
        dirname +="_kb_"+str(kappa_b)+"_kv_"+str(kappa_v)
        dirname +="_ka_"+str(kappa_a)+"_kc_"+str(kappa_c)
        dirname +="_kt_"+str(kappa_t)
        dirname +="_kr_"+str(kappa_r)
        dirname +="_lgvdamp_"+str(langevin_damp)
        dirname +="_exp_4"
        dirname +="_eqgcmc_"+str(equilibration_gcmc)
        dirname +="_rcmm_"+str(rc_mm)
        dirname +="_intstrgth_"+str(interaction_strength)
        dirname +="_geomfact_"+str(geometric_factor)
        dirname +="_probtrans_"+str(prob_transform)
        dirname +="_factors_"+str(variable_factor)

        fdir.writelines(dirname+"\n")
        dirname = os.path.join(dataset_name, dirname)
        # create the internal directory where we vary parameters
        isExist = os.path.exists(dirname)
        if not isExist:
            os.makedirs(dirname)

        # experiment dictionaries
        experiment_dictionary={
            'N_gcmc_1':N_gcmc_1,
            'X_gcmc_1':X_gcmc_1,
            'seed_gcmc_1':seed_gcmc_1,
            'mu_gcmc_1':mu_gcmc_1,
            'max_gcmc_1':max_gcmc_1,
            'xlo_1':xlo_1,
            'xhi_1':xhi_1,
            'ylo_1':ylo_1,
            'yhi_1':yhi_1,
            'zlo_1':zlo_1,
            'zhi_1':zhi_1,
            'sigma_metabolites':sigma_metabolites,
            'interaction_range_metabolites':interaction_range_metabolites,
            'interaction_strength_metabolites':interaction_strength_metabolites
        }
        dictionary_parameters_experiment=experiment_dictionary


        parameters = {
            'traj_steps':traj_steps,
            'flipratio':flipratio,
            'sigma':sigma,
            'kappa_b':kappa_b,
            'kappa_a':kappa_a,
            'kappa_c':kappa_c,
            'kappa_r':kappa_r,
            'kappa_t':kappa_t,
            'kappa_v':kappa_v,
            'step_size':step_size,
            'langevin_damp':langevin_damp,
            'langevin_seed':langevin_seed,
            'switch_mode':switch_mode,
            'total_sim_time':total_sim_time,
            'discret_snapshots':discret_snapshots,
            'xlo':xlo,
            'xhi':xhi,
            'ylo':ylo,
            'yhi':yhi,
            'zlo':zlo,
            'zhi':zhi,
            'equilibration_gcmc':equilibration_gcmc,
            'experiment_number': 4,
            'dictionary_parameters_experiment':dictionary_parameters_experiment
        }

        # create a parameter log in the directory
        filewriter_utils.create_parameter_log(dirname, parameters)

        # write the actual file
        filewriter_utils.write_trilmp_file(dirname+"/"+run_fname, parameters)

        # include mesh if needed
        filewriter_utils.copy_necessary_for_file_to_run(dirname)

        return

if __name__ == "__main__":
    main()