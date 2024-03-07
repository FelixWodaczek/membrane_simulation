from pathlib import Path
from dataclasses import dataclass, field

import trimesh

def icosphere(n):
    """Icosphere from trimesh."""
    s = trimesh.creation.icosphere(n)
    return s.vertices, s.faces

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


class BasePairStyle():
    def __init__(self, name: str, coeff_commands: list[str] = [], modify_commands: list[str] = []):
        self.name = name
        self.coeff_commands = coeff_commands
        self.modify_commands = modify_commands

    @property
    def init_params(self) -> list:
        return []

    def get_init_string(self) -> str:
        return ' '.join([self.name]+self.init_params)
    
    def parse_coeff_commands(self):
        return ['pair_coeff '+ coeff_command for coeff_command in self.coeff_commands]
    
    def parse_modify_commands(self):
        return ['pair_modify '+ modify_command for modify_command in self.modify_commands]

class TablePairStyle(BasePairStyle):
    def __init__(self, style: str = 'linear', N: int = 2000, *args, **kwargs):
        super().__init__(name='table', *args, **kwargs)
        self.style = style
        self.N = N

    @property
    def init_params(self) -> list:
        return [self.style, str(self.N)]
    
class LJCutPairStyle(BasePairStyle):
    def __init__(self, cutoff: float=2.5, *args, **kwargs):
        super().__init__(name='lj/cut', *args, **kwargs)
        self.cutoff = cutoff

    @property
    def init_params(self) -> list:
        return [str(self.cutoff)]
    
    def set_membrane_attraction(self, n_types, interaction_strength, sigma_tilde, interaction_range):
        self.coeff_commands = ['* * lj/cut 0 0 0']
        for i_regular in range(1, n_types):
            self.coeff_commands.append(f'1 {i_regular+1} lj/cut {interaction_strength} {sigma_tilde} {interaction_range}')
        
        self.modify_commands = ['pair lj/cut shift yes']

class HarmonicCutPairStyle(BasePairStyle):
    def __init__(self, *args, **kwargs):
        super().__init__(name='harmonic/cut', *args, **kwargs)

    def set_metabolite_repulsive_commands(self, n_types, sigma_metabolites):
        self.coeff_commands = ['* * harmonic/cut 0 0']
        # Add interactions between regular particles
        for i_regular in range(1, n_types):
            for j_regular in range(i_regular, n_types):
                self.coeff_commands.append(f"{i_regular+1} {j_regular+1} harmonic/cut 1000 {sigma_metabolites}")

    def set_all_repulsive_commands(self, n_types, sigma_metabolites, sigma_special):
        self.set_metabolite_repulsive_commands(n_types=n_types, sigma_metabolites=sigma_metabolites)

        # Add interactions with special sphere
        for i_regular in range(1, n_types):
            self.coeff_commands.append(f"1 {i_regular+1} harmonic/cut 1000 {sigma_special}")

@dataclass
class BondCreation():
    name: str
    Nevery: int
    Rmin: float
    itype: int = 1
    jtype: int = 2
    bondtype: int = 2
    target_class: str = 'all'
    add_args: dict = field(default_factory=dict)

    def bond_creation_command(self) -> str:
        return f"fix {self.name} {self.target_class} bond/create {self.Nevery} {self.itype} {self.jtype} {self.Rmin} {self.bondtype} {' '.join([f"{k} {v}" for k, v in self.add_args.items()])}"

@dataclass
class BondDeletion():
    name: str
    Nevery: int
    Rmax: float = 0.
    bondtype: int = 2
    prob: float = 1.
    seed: int = 123
    add_args: dict = field(default_factory=dict)

    def bond_deletion_command(self) -> str:
        return f"fix {self.name} all bond/break {self.Nevery} {self.bondtype} {self.Rmax} prob {self.prob} {self.seed} {' '.join([f"{k} {v}" for k, v in self.add_args.items()])}"

class Reaction():
    def __init__(
        self, 
        name: str, 
        pretransform_template: Path, posttransform_template: Path, map_template: Path,
        Nevery: float, Rmin: float, Rmax: float, prob:float, seed: int=123
    ):
        self.name = name
        self.pretransform_template = Path(pretransform_template)
        self.posttransform_template = Path(posttransform_template)
        self.map_template = Path(map_template)
        self.Nevery = Nevery
        self.Rmin = Rmin
        self.Rmax = Rmax
        self.prob = prob
        self.seed = seed
    
    def get_molecule_commands(self) -> list[str]:
        return [
            f"molecule mpre{self.name} {self.pretransform_template}",
            f"molecule mpost{self.name} {self.posttransform_template}",
        ]
    
    def get_reaction_string(self) -> str:
        return f"react {self.name} all {self.Nevery} {self.Rmin} {self.Rmax} mpre{self.name} mpost{self.name} {self.map_template} prob {self.prob} {self.seed}"