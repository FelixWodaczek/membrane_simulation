from itertools import product

import numpy as np

import membrane_simulation.mgradient_classes as cg

def main():
    launch_directories = cg.prepare_bash_launch()

    main_path = "control"
    index_run = cg.check_num_runs(main_path)
    total_simulations = 0

    varied_parameters = ['langevin_seed', 'channel_height']
    langevin_seeds = [123]
    channel_heights = np.round(np.arange(0.1, 1.1, 0.1), 2)

    # MD variables
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
    metabolite_type = 3

    for langevin_seed, channel_height in product(langevin_seeds, channel_heights):
        total_simulations += 1
        
        resolution = 5
        
        sim = cg.SimLauncher_HeterogeneousMembrane(
            langevin_seed=langevin_seed, 
            total_sim_time=total_sim_time, 
            discret_snapshots=discret_snapshots, 
            print_program_iterations=print_program_iterations,
            eq_rounds=eq_rounds,
            xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, zlo=zlo, zhi=zhi,
            kappa_a=kappa_a,
            kappa_v=kappa_v,
            resolution=resolution,
            fix_COM=fix_COM,
            channel_height=channel_height,
        )

        # DEFINE GEOMETRY OF THE SYSTEM
        geometry_factor = channel_height # <--- to change, size of the metabolite channel
        dimensions_boxes = sim.radius_membrane*geometry_factor
        xlo_source = xhi -12 # <--- to gauge: position of the source
        xhi_source = xhi -2  # <------------  
        xlo_sink   = -40     # <--- to change: sink should probably be a sphere (to adapt it to the vesicle)
        xhi_sink   = 0       # <---- atm it is also a rectangle (although we don't expect this to be too important)

        # DEFINE CHEMOSTAT PROPERTIES
        N_gcmc = int(10*sim.langevin_damp/sim.step_size) # do not change
        X_gcmc = 100 # can be changed (the bigger, the more metabolites you should have - also, costlier sims)
        seed_gcmc = langevin_seed*2
        mugcmc = 0 # can be changed (the bigger, the more metabolitres)

        # INPUT ALL THE PARAMETERS
        params = {}
        params['interaction_membrane_metabolites'] = 1.0 # keep it at 1
        params['rc_membrane_metabolites'] = 1.5 # keep at 1.5

        gradient_commands = f"""region SOURCE cylinder x 0 0 {dimensions_boxes} {xlo_source} {xhi_source} side in
        region SINK cylinder x 0 0 {dimensions_boxes} {xlo_sink} {xhi_sink} side in
        fix GCMCSOURCE metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {seed_gcmc} {sim.temperature} {mugcmc} 0 region SOURCE
        fix GCMCSINK metabolites gcmc {N_gcmc} {X_gcmc} 0 {metabolite_type} {11*seed_gcmc} {sim.temperature} -100 0 region SINK
        region CONFINEMENT cylinder x 0 0 {dimensions_boxes+2} EDGE EDGE side in
        fix CONFWALL metabolites wall/region CONFINEMENT harmonic 5000 1 1
        region CHANNEL cylinder x 0 0 {dimensions_boxes+2} {sim.edge_vesicle} EDGE side in"""

        params['gradient_commands'] = cg.append_lines_to_list(gradient_commands)
        
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
        fix  aveMET all ave/time {sim.print_frequency} 1 {sim.print_frequency} v_nmet c_TempCompute file "metabolite_properties.dat"
        fix  aveMETBin all ave/time {sim.print_frequency} 1 {sim.print_frequency} c_BinMetX[*] file "metabolite_binning.dat" mode vector
        fix  aveMEMCHANNEL all ave/time {sim.print_frequency} 1 {sim.print_frequency} v_meminchannel file "membraneCHANNEL.dat"
        compute FORCES vertices group/group metabolites pair yes
        fix aveFORCES all ave/time {sim.print_frequency} 1 {sim.print_frequency} c_FORCES c_FORCES[1] c_FORCES[2] c_FORCES[3] file "metabolite_forces.dat"
        """
        
        # add additional fixes
        params['additional_fixes'] = cg.append_lines_to_list(additional_fixes)

        sim.initialize_commands_gradient(params)

        path = cg.create_directory_structure(main_path+"/run"+str(int(index_run))+"/", sim, varied_parameters, launch_directories)
        sim.write_function(path+'/launch.py')

if __name__ == '__main__':
    main()