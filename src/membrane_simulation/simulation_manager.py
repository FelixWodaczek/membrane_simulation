from pathlib import Path
from dataclasses import field

import numpy as np

from trimem.mc.trilmp import TriLmp, Beads
from trimesh import Trimesh

import membrane_simulation.simulation_utils as sutils 

class SimulationManager():
    def __init__(
            self,
            trimem_params: sutils.TrimemParameters = field(default_factory=sutils.TrimemParameters),
            membrane_params: sutils.MembraneParameters = field(default_factory=sutils.MembraneParameters),
            mesh: int | Path = 2, 
        ):
        self.trimem_params: sutils.TrimemParameters = trimem_params
        self.membrane_params: sutils.MembraneParameters = membrane_params

        if isinstance(mesh, int):
            vertices, faces = sutils.icosphere(mesh)
        elif isinstance(mesh, Path):
            import pandas as pd
            vertices = pd.read_csv(Path(mesh).resolve().joinpath('mesh_coordinates_N_5072_.dat'), header = None, index_col = False, sep = ' ')[[1, 2, 3]].to_numpy()
            faces = pd.read_csv(Path(mesh).resolve().joinpath('mesh_faces_N_5072_.dat'), header = None, index_col = False, sep = ' ')[[0, 1, 2]].to_numpy()
        else:
            raise TypeError(f'mesh must be either an int or a Path object but is {mesh} of type {type(mesh)}')
        
        self.mesh = Trimesh(vertices=vertices, faces=faces)
        # rescaling it so that we start from the right distances
        desired_average_distance = 2**(1.0/6.0) * self.membrane_params.sigma_vertex
        current_average_distance = np.mean(self.mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        self.mesh.vertices *= scaling

        self.pair_styles: list[sutils.BasePairStyle] = []
        self.reactions: list[sutils.Reaction] = []

    def init_trilmp(
            self,
            simulation_box: sutils.Box = None,
            debug_mode: bool = False,
            metabolite_masses: list[float] = None,
        ):
        if simulation_box is None:
            self.box = sutils.Box()
        else:
            self.box = simulation_box

        if metabolite_masses is None:
            metabolite_masses = [1.0] * (self.trimem_params.n_groups - 1)
        assert len(metabolite_masses) == self.trimem_params.n_groups - 1, f"metabolite_masses must have length {self.trimem_params.n_groups - 1} but has length {len(metabolite_masses)}"

        self.trilmp = TriLmp(
            initialize=True,                          # use mesh to initialize mesh reference
            debug_mode=debug_mode,

            mesh_points=self.mesh.vertices,           # input mesh vertices
            mesh_faces=self.mesh.faces,               # input of the mesh faces
            kappa_b=self.membrane_params.kappa_b,     # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a=self.membrane_params.kappa_a,     # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v=self.membrane_params.kappa_v,     # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c=self.membrane_params.kappa_c,     # MEMBRANE MECHANICS: constraint on area difference change (understand meaning) (kB T)
            kappa_t=self.membrane_params.kappa_t,     # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r=self.membrane_params.kappa_r,     # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            
            num_particle_types=self.trimem_params.n_groups,                       # how many particle types will there be in the system
            mass_particle_type=[self.membrane_params.vertex_mass]+metabolite_masses,# the mass of the particle per type
            group_particle_type=self.trimem_params.particle_groups,                 # the names of the groups
            n_bond_types=1,

            step_size=self.trimem_params.step_size,                      # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps=self.trimem_params.traj_steps,                    # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio=self.trimem_params.flip_ratio,                    # MC PART SIMULATION: fraction of edges to flip?
            initial_temperature=self.trimem_params.initial_temperature,          # MD PART SIMULATION: temperature of the system
            switch_mode='random',
            # check_neigh_every=1,                      # NEIGHBOUR LISTS
            equilibration_rounds=self.trimem_params.equilibriation_rounds,                   # MEMBRANE EQUILIBRATION ROUNDS
            pure_MD=self.trimem_params.pure_MD,                          # MD PART SIMULATION: accept every MD trajectory?
            
            box=self.box.to_tuple(),                       # MD PART SIMULATION: box dimensions

            info = self.trimem_params.print_frequency,
            thin = self.trimem_params.discrete_snapshots, # OUTPUT: interval of snapshots
            checkpoint_every=100*self.trimem_params.print_frequency,         # OUTPUT: interval of checkpoints (alternating pickles)
            performance_increment=self.trimem_params.print_program_iterations,    # OUTPUT: output performace stats to prefix_performance.dat file
            energy_increment=self.trimem_params.print_program_iterations,         # OUTPUT: output energies to energies.dat file
        )

    @staticmethod
    def equilibriation_thermostat_command(initial_temperature, langevin_damp, langevin_seed: int, fix_gcmc: bool=True):
        commands = []

        if fix_gcmc:
            # necessary to include the nve because we have removed it
            commands.append(f'fix vertexnve vertices nve')
            commands.append(f'fix lvt vertices langevin {initial_temperature} {initial_temperature}  {langevin_damp} {langevin_seed} zero yes tally yes')
            return '\n'.join(commands)
            
        commands.append(f"fix lvt vertices langevin {initial_temperature} {initial_temperature}  {langevin_damp} {langevin_seed} zero yes")
        return '\n'.join(commands)

    @staticmethod
    def langevin_commands(membrane_vertex_mass, initial_temperature: float, langevin_damp, langevin_seed: int, beads: Beads = None, length_scale: float = 1.):
        fix_gcmc = True

        commands = []
        sc0=membrane_vertex_mass/length_scale

        bds_info_lgv = ''

        if beads is not None:
            if beads.n_types>1:
                bds_info_lgv = ''
                for i in range(beads.n_types):
                    bds_info_lgv += f'scale {i + 2} {(beads.masses[i] / beads.bead_sizes[i])/sc0} '
            else:
                bds_info_lgv=f'scale 2 {(beads.masses / beads.bead_sizes)/sc0}'

            # FIX GCMC (scale the impact of the thermostat on the metabolites - beads masses is by default 1)
            if fix_gcmc:
                bds_info_lgv+= f'scale 2 {(beads.masses / beads.bead_sizes[0])/sc0}'
            

        commands.append("fix mynve all nve")
        commands.append(f"fix lvt all langevin {initial_temperature} {initial_temperature} {langevin_damp} {langevin_seed} zero yes {bds_info_lgv}")
        return '\n'.join(commands)

    @staticmethod
    def walls_command(ylo, yhi, zlo, zhi):
        return f"fix ConfinementMet metabolites wall/reflect xhi EDGE ylo {ylo} yhi {yhi} zlo {zlo} zhi {zhi}"

    def get_pair_style_commands(self):
        # This command tells lamps all the styles that are to be expected
        style_command = f"pair_style hybrid/overlay"
        coeff_and_modify_commands = []
        for pair_style in self.pair_styles:
            # add parameter and initial parameters
            style_command = ' '.join([style_command, pair_style.get_init_string()])
            coeff_and_modify_commands += pair_style.parse_coeff_commands() + pair_style.parse_modify_commands()

        return '\n'.join([style_command]+coeff_and_modify_commands)
            
    def get_chemistry_commands(self):
        molecule_template_commands = []
        chemistry_strings = ["fix freact all bond/react reset_mol_ids no"]
        for reaction in self.reactions:
            molecule_template_commands += reaction.get_molecule_commands()
            chemistry_strings.append(reaction.get_reaction_string())

        return '\n'.join(molecule_template_commands + [' '.join(chemistry_strings)])
        
def main():
    smanager = SimulationManager()
    smanager.init_trilmp()
    # out = smanager.run()

if __name__ == '__main__':
    main()