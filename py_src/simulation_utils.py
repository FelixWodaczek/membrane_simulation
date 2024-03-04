from pathlib import Path
import trimesh

def icosphere(n):
    """Icosphere from trimesh."""
    s = trimesh.creation.icosphere(n)
    return s.vertices, s.faces

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