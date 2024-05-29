# COMMENTED VERSION OF THE CODE FOR FELIX W.

from ClassesGradient import *

if __name__=="__main__":

    """
    Basic gradient simulations
    - Source and sink pair
    - Metabolites confined by walls

    Heterogeneous membrane simulations
    - Source and sink pair
    - Metabolites confined by walls
    - Only a portion of the membrane can interact
    """

    # prepare filenames for slurm launch
    launch_directories = prepare_bash_launch()

    # directory names etc
    main_path = "control"
    index_run = check_num_runs(main_path)
    total_simulations = 0
    gradients = False 

    if gradients:

        # ------------------------ insert loop here of what you want to iterate
        varied_parameters = ['langevin_seed', 'fraction_heterogeneous']

        # to simulate heterogeneous membrane (might not apply for first tests)
        for fh in [0.05, 0.1, 0.2]:
            
            # launch different seeds
            for ss in [123, 456, 789]:
                
                # count number of simulations launched
                total_simulations+=1

                # mesh size (keep it at 5 for 10k beads)
                resolution = 5
                
                # MD variables
                langevin_seed = ss
                total_sim_time = 5000 # units of time
                discret_snapshots = 20 # units of time
                print_program_iterations = 100 # time steps
                eq_rounds = 10 # must be integer number

                # (you can change this) simulation box
                xlo = -40
                xhi = 70
                ylo = -40
                yhi = 40
                zlo = -40
                zhi = 40

                # (you can change this too) membrane properties (see other default parameters in ClassesGradient class)
                kappa_a = 2.5e5
                kappa_v = 2.5e5
                
                # decide whether to fix ('') or not (#) the COM
                fix_COM = ''

                # do not touch this value - the program takes care of it
                metabolite_type = 2

                # heterogeneous membrane, fraction of beads of other kind
                heterogeneous_membrane = False
                fraction_heterogeneous = fh

                if not heterogeneous_membrane:
                    # initialize sim class
                    Simulation = SimLauncher_BasicGradient(langevin_seed=langevin_seed, 
                                                total_sim_time=total_sim_time, 
                                                discret_snapshots=discret_snapshots, 
                                                print_program_iterations=print_program_iterations,
                                                eq_rounds=eq_rounds,
                                                xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi,
                                                kappa_a=kappa_a,
                                                kappa_v=kappa_v,
                                                resolution=resolution,
                                                fix_COM=fix_COM)

                if heterogeneous_membrane:
                    metabolite_type+=1
                    # initialize sim class
                    Simulation = SimLauncher_HeterogeneousMembrane(langevin_seed=langevin_seed, 
                                                total_sim_time=total_sim_time, 
                                                discret_snapshots=discret_snapshots, 
                                                print_program_iterations=print_program_iterations,
                                                eq_rounds=eq_rounds,
                                                xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi,
                                                kappa_a=kappa_a,
                                                kappa_v=kappa_v,
                                                resolution=resolution,
                                                fix_COM=fix_COM,
                                                fraction_heterogeneous=fraction_heterogeneous)
                

                # DEFINE GEOMETRY OF THE SYSTEM
                geometry_factor = 0.25 # <--- to change, size of the metabolite channel
                dimensions_boxes = Simulation.radius_membrane*geometry_factor
                xlo_source = xhi -12 # <--- to gauge: position of the source
                xhi_source = xhi -2  # <------------   ""
                xlo_sink   = -40     # <--- to change: sink should probably be a sphere (to adapt it to the vesicle)
                xhi_sink   = 0       # <---- atm it is also a rectangle (although we don't expect this to be too important)

                # DEFINE CHEMOSTAT PROPERTIES
                N_gcmc = int(10*Simulation.langevin_damp/Simulation.step_size) # do not change
                X_gcmc = 100 # can be changed (the bigger, the more metabolites you should have - also, costlier sims)
                seed_gcmc = langevin_seed*2
                mugcmc = 0 # can be changed (the bigger, the more metabolitres)

                # INPUT ALL THE PARAMETERS
                params = {}
                params['interaction_membrane_metabolites'] = 1.0 # keep it at 1
                params['rc_membrane_metabolites'] = 1.5 # <----- this should go from 1 to 10

                use_blocks = False   # rectangular chanel
                use_cylinders = True # cylindrical chanel

                # ideally here below the shape of the sink would be changed to a sphere

                if use_blocks:
                    gradient_commands = f"""region SOURCE block {xlo_source} {xhi_source} {-dimensions_boxes} {dimensions_boxes} {-dimensions_boxes} {dimensions_boxes} side in
                    region SINK block {xlo_sink} {xhi_sink} {-dimensions_boxes} {dimensions_boxes} {-dimensions_boxes} {dimensions_boxes} side in
                    fix GCMCSOURCE metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {seed_gcmc} {Simulation.temperature} {mugcmc} 0 region SOURCE
                    fix GCMCSINK metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {11*seed_gcmc} {Simulation.temperature} -100 0 region SINK
                    region CONFINEMENT block EDGE EDGE {-dimensions_boxes-2} {dimensions_boxes+2} {-dimensions_boxes-2} {dimensions_boxes+2} side in
                    fix CONFWALL metabolites wall/region CONFINEMENT harmonic 5000 1 1
                    region CHANNEL block {Simulation.edge_vesicle} EDGE {-dimensions_boxes-2} {dimensions_boxes+2} {-dimensions_boxes-2} {dimensions_boxes+2} side in"""
                    
                elif use_cylinders:
                    gradient_commands = f"""region SOURCE cylinder x 0 0 {dimensions_boxes} {xlo_source} {xhi_source} side in
                    region SINK cylinder x 0 0 {dimensions_boxes} {xlo_sink} {xhi_sink} side in
                    fix GCMCSOURCE metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {seed_gcmc} {Simulation.temperature} {mugcmc} 0 region SOURCE
                    fix GCMCSINK metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {11*seed_gcmc} {Simulation.temperature} -100 0 region SINK
                    region CONFINEMENT cylinder x 0 0 {dimensions_boxes+2} EDGE EDGE side in
                    fix CONFWALL metabolites wall/region CONFINEMENT harmonic 5000 1 1
                    region CHANNEL cylinder x 0 0 {dimensions_boxes+2} {Simulation.edge_vesicle} EDGE side in"""
                    
                # introduce gradient commands
                params['gradient_commands'] = append_lines_to_list(gradient_commands)

                # different LAMMPS calculations that make life easy
                # you can add anything you see fitting
                additional_fixes = f"""variable meminchannel equal count(vertices,CHANNEL)
                compute TempCompute metabolites temp
                compute_modify TempCompute dynamic/dof yes
                variable    typemet atom "type=={metabolite_type}"
                group       typemet dynamic metabolites var typemet
                variable    nmet equal count(typemet)
                compute BinMet metabolites chunk/atom bin/1d x lower 10
                compute BinMetX metabolites property/chunk BinMet count coord1
                fix  aveMET all ave/time {Simulation.print_frequency} 1 {Simulation.print_frequency} v_nmet c_TempCompute file "metabolite_properties.dat"
                fix  aveMETBin all ave/time {Simulation.print_frequency} 1 {Simulation.print_frequency} c_BinMetX[*] file "metabolite_binning.dat" mode vector
                fix  aveMEMCHANNEL all ave/time {Simulation.print_frequency} 1 {Simulation.print_frequency} v_meminchannel file "membraneCHANNEL.dat"
                """
                
                # add additional fixes
                params['additional_fixes'] = append_lines_to_list(additional_fixes)

                # initialize the gradient commands
                Simulation.initialize_commands_gradient(params)

                # generate file
                path = create_directory_structure(main_path+"/run"+str(int(index_run))+"/", Simulation, varied_parameters, launch_directories)
                print(f"PATH: {path}")
                Simulation.write_function(path+"/launch.py")

    if not gradients:
        # ------------------------ insert loop here of what you want to iterate
        varied_parameters = ['langevin_seed']

        for ss in [123, 456, 789]:
            
            # count number of simulations launched
            total_simulations+=1

            # mesh size
            resolution = 6
            
            # MD variables
            langevin_seed = ss
            total_sim_time = 10000 # units of time
            discret_snapshots = 50 # units of time
            print_program_iterations = 100
            eq_rounds = 1000

            # simulation box
            xlo = -40
            xhi = 40
            ylo = -40
            yhi = 40
            zlo = -40
            zhi = 40

            # membrane properties
            kappa_a = 2.5e5
            kappa_v = 2.5e5
            
            # decide whether to fix ('') or not (#) the COM
            fix_COM = '#'

            # initialize sim class
            Simulation = SimLauncher_Vesicle(langevin_seed=langevin_seed, 
                                        total_sim_time=total_sim_time, 
                                        discret_snapshots=discret_snapshots, 
                                        print_program_iterations=print_program_iterations,
                                        eq_rounds=eq_rounds,
                                        xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi,
                                        kappa_a=kappa_a,
                                        kappa_v=kappa_v,
                                        resolution=resolution,
                                        fix_COM=fix_COM)

            # generate file
            path = create_directory_structure(main_path+"/run"+str(int(index_run))+"/", Simulation, varied_parameters, launch_directories)
            print(f"PATH: {path}")
            Simulation.write_function(path+"/launch.py")

    print(f"You are launching {total_simulations} simulations.")
    
