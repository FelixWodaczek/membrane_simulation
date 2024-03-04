from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from trimem.mc.trilmp import TriLmp, Beads
from trimesh import Trimesh

import py_src.simulation_utils as simulation_utils

@dataclass
class MembraneParams():
    sigma_vertex: float = 1.0
    vertex_mass: float = 1.0
    kappa_b: float = 20.0
    kappa_a: float = 2.5e5
    kappa_v: float = 2.5e5
    kappa_c: float = 0.0
    kappa_t: float = 1.0e4
    kappa_r: float = 1.0e4

@dataclass
class Box():
    xlo: float = -50
    xhi: float = 50
    ylo: float = -50
    yhi: float = 50
    zlo: float = -50
    zhi: float = 50

    def to_tuple(self):
        return (self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi)
    
    def to_list(self):
        return [self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi]
    
    def __str__(self):
        return ' '.join([str(x) for x in self.to_tuple()])

@dataclass
class GCMCRegion():
    name: str
    N: int
    maxp: int
    X: int = 100
    T: float = 1.0
    mu: float = 0
    type: int = 2
    seed: int = 123
    box: Box = field(default_factory=Box)

    def region_command(self):
        return f"region gcmc_region_{self.name} block {self.box} side in"

    def fix_gcmc_command(self):
        return f"fix mygcmc_{self.name} metabolites gcmc {self.N} {self.X} 0 {self.type} {self.seed} {self.T} {self.mu} 0 region gcmc_region_{self.name} max {self.maxp}"

class SimulationManager():
    def __init__(self, resolution=2, membrane_params: MembraneParams = MembraneParams()):
        self.membrane_params = membrane_params

        self.resolution = resolution

        vertices, faces = simulation_utils.icosphere(resolution)
        self.mesh = Trimesh(vertices=vertices, faces=faces)
        # rescaling it so that we start from the right distances
        desired_average_distance = 2**(1.0/6.0) * self.membrane_params.sigma_vertex
        current_average_distance = np.mean(self.mesh.edges_unique_length)
        scaling = desired_average_distance/current_average_distance
        self.mesh.vertices *= scaling
        self.postequilibration_lammps_commands = []

        self.pair_styles: list[simulation_utils.BasePairStyle] = []
        self.reactions: list[simulation_utils.Reaction] = []
        self.template_path = Path(__file__).resolve().parent.joinpath('reaction_templates/')

    def init_trilmp(self, box: Box = Box(xlo=-50, xhi=50, ylo=-50, yhi=50, zlo=-50, zhi=50)):
        # trimem parameters
        self.traj_steps=50
        self.flip_ratio=0.1
        self.step_size=0.001
        self.total_sim_time = 1000 # 00  # in time units
        self.discrete_snapshots = 10   # in time units
        self.print_frequency = int(self.discrete_snapshots/(self.step_size*self.traj_steps))

        self.initial_temperature=1.0                    # MD PART SIMULATION: temperature of the system
        pure_MD=False,                             # MD PART SIMULATION: accept every MD trajectory?

        self.langevin_damp=1.0
        self.langevin_seed=123

        self.box = box
        switch_mode = 'random'

        # Interaction parameters
        self.interaction_range = 1.5 # rc_mm
        self.interaction_strength = 10
        self.sigma_metabolites = 1.0
        self.sigma_tilde_membrane_metabolites = 0.5*(self.sigma_metabolites+self.membrane_params.sigma_vertex)
        self.interaction_range_tilde = self.interaction_range*self.sigma_tilde_membrane_metabolites

        self.trilmp = TriLmp(
            initialize=True,                          # use mesh to initialize mesh reference
            mesh_points=self.mesh.vertices,           # input mesh vertices
            mesh_faces=self.mesh.faces,               # input of the mesh faces
            kappa_b=self.membrane_params.kappa_b,                          # MEMBRANE MECHANICS: bending modulus (kB T)
            kappa_a=self.membrane_params.kappa_a,                          # MEMBRANE MECHANICS: constraint on area change from target value (kB T)
            kappa_v=self.membrane_params.kappa_v,                          # MEMBRANE MECHANICS: constraint on volume change from target value (kB T)
            kappa_c=self.membrane_params.kappa_c,                          # MEMBRANE MECHANICS: constraint on area difference change (understand meaning) (kB T)
            kappa_t=self.membrane_params.kappa_t,                          # MEMBRANE MECHANICS: tethering potential to constrain edge length (kB T)
            kappa_r=self.membrane_params.kappa_r,                          # MEMBRANE MECHANICS: repulsive potential to prevent surface intersection (kB T)
            
            num_particle_types=3,                       # how many particle types will there be in the system
            mass_particle_type=[self.membrane_params.vertex_mass, 1., 1.],# the mass of the particle per type
            group_particle_type=['vertices', 'metabolites', 'waste'],

            step_size=self.step_size,                      # FLUIDITY ---- MD PART SIMULATION: timestep of the simulation
            traj_steps=self.traj_steps,                    # FLUIDITY ---- MD PART SIMULATION: number of MD steps before bond flipping
            flip_ratio=self.flip_ratio,                    # MC PART SIMULATION: fraction of edges to flip?
            # check_neigh_every=1,                      # NEIGHBOUR LISTS
            equilibration_rounds=1,                   # MEMBRANE EQUILIBRATION ROUNDS

            box=self.box.to_tuple(),                       # MD PART SIMULATION: box dimensions

            output_prefix='data/data',                # OUTPUT: prefix for output filenames
            restart_prefix='data/data',               # OUTPUT: name for checkpoint files
            checkpoint_every=self.print_frequency,         # OUTPUT: interval of checkpoints (alternating pickles)
            output_format='lammps_txt',               # OUTPUT: choose different formats for 'lammps_txt', 'lammps_txt_folder' or 'h5_custom'
            output_counter=0,                         # OUTPUT: initialize trajectory number in writer class
            performance_increment=self.print_frequency,    # OUTPUT: output performace stats to prefix_performance.dat file
            energy_increment=self.print_frequency,         # OUTPUT: output energies to energies.dat file
            pure_MD=pure_MD,                          # MD PART SIMULATION: accept every MD trajectory?
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