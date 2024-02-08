import os, filewriter_utils

if __name__=="__main__":

    # general simulation variables
    filename = "trilmp_gcmc.py"
    dataset_name = "experiment2_varying_factors_distance"

    # create header directory
    isExist = os.path.exists(dataset_name)
    if not isExist:
        os.makedirs(dataset_name)

    # fdir
    fdir = open(dataset_name+"/directories_in_directory.dat", "w")

    variable_distance = [50, 100]
    variable_factors = [1, 5, 10]

    for vf in variable_factors:

        for vd in variable_distance:

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
            xlo, xhi = -vd, vd
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

            # ---------------------------
            #         EXPERIMENTS
            # ---------------------------

            experiment_number = 2

            # ..... experiment 1
            if experiment_number==1:
                factor_introduce=vf

                N_gcmc_1=int(factor_introduce*langevin_damp/step_size)
                X_gcmc_1=100
                seed_gcmc_1 = 123
                mu_gcmc_1=0
                max_gcmc_1=1000000
                xlo_1 = -10
                xhi_1 = 10
                ylo_1 = -10
                yhi_1 = 10
                zlo_1 = -10
                zhi_1 = 10

                N_gcmc_2=int(factor_introduce*langevin_damp/step_size)
                X_gcmc_2=100
                seed_gcmc_2 = 123
                mu_gcmc_2=-100
                max_gcmc_2=0
                xlo_2 = xlo
                xhi_2 = xhi
                ylo_2 = ylo
                yhi_2 = yhi
                zlo_2 = zlo
                zhi_2 = zhi

                sigma_metabolites = 1.0
                interaction_range_metabolites = 2.5
                interaction_strength_metabolites = 10

            # ..... experiment 2
            if experiment_number==2:
                factor_introduce=vf

                N_gcmc_1=int(factor_introduce*langevin_damp/step_size)
                X_gcmc_1=100
                seed_gcmc_1 = 123
                mu_gcmc_1=0
                max_gcmc_1=1000000
                xlo_1 = vd - 10
                xhi_1 = vd
                ylo_1 = ylo
                yhi_1 = yhi
                zlo_1 = zlo
                zhi_1 = zhi

                N_gcmc_2=int(factor_introduce*langevin_damp/step_size)
                X_gcmc_2=100
                seed_gcmc_2 = 123
                mu_gcmc_2=-100
                max_gcmc_2=0
                xlo_2 = -vd
                xhi_2 = -vd + 10
                ylo_2 = ylo
                yhi_2 = yhi
                zlo_2 = zlo
                zhi_2 = zhi

                sigma_metabolites = 1.0
                interaction_range_metabolites = 2.5
                interaction_strength_metabolites = 10

            # ..... experiment 3
            if experiment_number==3:
                factor_introduce=vf

                N_gcmc_1=int(factor_introduce*langevin_damp/step_size)
                X_gcmc_1=100
                seed_gcmc_1 = 123
                mu_gcmc_1=0
                max_gcmc_1=1000000
                xlo_1 = -10
                xhi_1 = 10
                ylo_1 = -10
                yhi_1 = 10
                zlo_1 = -10
                zhi_1 = 10

                sigma_metabolites = 1.0
                interaction_range_metabolites = 2.5
                interaction_strength_metabolites = 10

            # -------------------------------
            #    WRITE THE DIRECTORY NAME
            # -------------------------------

            dirname = "data_trajsteps_"+str(traj_steps)
            dirname +="_flipratio_"+str(flipratio)
            dirname +="_stepsize_"+str(step_size)
            dirname +="_seed_"+str(langevin_seed)
            dirname +="_kb_"+str(kappa_b)+"_kv_"+str(kappa_v)
            dirname +="_ka_"+str(kappa_a)+"_kc_"+str(kappa_c)
            dirname +="_kt_"+str(kappa_t)
            dirname +="_kr_"+str(kappa_r)
            dirname +="_lgvdamp_"+str(langevin_damp)
            dirname +="_exp_"+str(experiment_number)
            dirname +="_eqgcmc_"+str(equilibration_gcmc)
            dirname +="_factors_"+str(vf)
            dirname +="_distance_"+str(vd)

            fdir.writelines(dirname+"\n")
            dirname = dataset_name +dirname
            # create the internal directory where we vary parameters
            isExist = os.path.exists(dirname)
            if not isExist:
                os.makedirs(dirname)

            # -----------------------------------------------------------
            # INPUT DATA
            # -----------------------------------------------------------

            dictionary_parameters_experiment = None

            if experiment_number == 1 or experiment_number == 2:
                # experiment dictionaries
                experiment_dictionary={'N_gcmc_1':N_gcmc_1,
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
                                        'N_gcmc_2':N_gcmc_2,
                                        'X_gcmc_2':X_gcmc_2,
                                        'seed_gcmc_2':seed_gcmc_2,
                                        'mu_gcmc_2':mu_gcmc_2,
                                        'max_gcmc_2':max_gcmc_2,
                                        'xlo_2':xlo_2,
                                        'xhi_2':xhi_2,
                                        'ylo_2':ylo_2,
                                        'yhi_2':yhi_2,
                                        'zlo_2':zlo_2,
                                        'zhi_2':zhi_2,
                                        'sigma_metabolites':sigma_metabolites,
                                        'interaction_range_metabolites':interaction_range_metabolites,
                                        'interaction_strength_metabolites':interaction_strength_metabolites}

                dictionary_parameters_experiment=experiment_dictionary

            if experiment_number == 3:
                # experiment dictionaries
                experiment_dictionary={'N_gcmc_1':N_gcmc_1,
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
                                        'interaction_strength_metabolites':interaction_strength_metabolites}

                dictionary_parameters_experiment=experiment_dictionary

            parameters = {'traj_steps':traj_steps,
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
                        'experiment_number':experiment_number,
                        'dictionary_parameters_experiment':dictionary_parameters_experiment}

            # -----------------------------------------------------------
            # PROGRAM
            # -----------------------------------------------------------

            # create a parameter log in the directory
            filewriter_utils.create_parameter_log(dirname, parameters)

            # write the actual file
            filewriter_utils.write_trilmp_file(dirname+"/"+filename, parameters)

            # include mesh if needed
            filewriter_utils.copy_necessary_for_file_to_run(dirname)

    fdir.close()
