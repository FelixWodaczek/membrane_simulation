import argparse
from pathlib import Path

class ChemistryWriter():
    def __init__(self, target_dir: Path):
        self.target_dir = Path(target_dir).resolve()

    def write_single_chemistry(self, name: str, n_cluster: int, type_cluster: int, start_type_single: int, end_type_single: int = None, delete_mode: bool = False):
        if end_type_single is None:
            end_type_single = start_type_single + 1

        fname = f'{name}_nneigh{n_cluster-1}.txt'
        
        with open(self.target_dir.joinpath('pre_'+fname), 'w') as f:
            f.write(f'# Pre-Reaction template for a network of one central atom with {n_cluster-1} neighbours\n\n')
            f.write(f'{n_cluster+1} atoms\n')
            f.write(f'{n_cluster-1} bonds\n\n')
            
            f.write('Types\n\n')
            # write single atom
            f.write(f'1 {start_type_single}\n')
            # write central and surrounding atoms
            f.write('\n'.join([f'{i_type+2} {type_cluster}' for i_type in range(n_cluster)]))
            f.write('\n\n')

            # write bonds for all neighbours of central atom (with index 2)
            f.write('Bonds\n\n')
            f.write('\n'.join(f'{i_bond+1} 1 2 {i_bond+3}' for i_bond in range(n_cluster-1)))

            f.close()
        
        with open(self.target_dir.joinpath('map_'+fname), 'w') as f:
            f.write(f"Map for reaction of changing interacting atom with central atom with {n_cluster-1} neighbours to other type\n\n")
            f.write(f'{n_cluster+1} equivalences\n')
            f.write(f'{n_cluster-1} edgeIDs\n')
            if delete_mode:
                f.write('1 deleteIDs\n')
            f.write('\n')

            f.write('InitiatorIDs\n\n')
            f.write('1\n')
            f.write('2\n\n')

            # Set that all other atoms in cluster have more connections (rest of membrane)
            f.write('EdgeIDs\n\n')
            f.write('\n'.join(f'{i_cluster+3}' for i_cluster in range(n_cluster-1)))
            f.write('\n\n')

            f.write("Equivalences\n\n")
            f.write('\n'.join(f'{i_atom+1} {i_atom+1}' for i_atom in range(n_cluster+1)))
            
            if delete_mode:
                f.write('\n\n')
                f.write("DeleteIDs\n\n")
                f.write('1')

            f.close()


        with open(self.target_dir.joinpath('post_'+fname), 'w') as f:
            f.write(f'# Post-Reaction template for a network of one central atom with {n_cluster-1} neighbours\n\n')
            f.write(f'{n_cluster+1} atoms\n')
            f.write(f'{n_cluster-1} bonds\n\n')
            
            f.write('Types\n\n')
            # write single atom
            f.write(f'1 {end_type_single}\n')
            # write central and surrounding atoms
            f.write('\n'.join([f'{i_type+2} {type_cluster}' for i_type in range(n_cluster)]))
            f.write('\n\n')

            # write bonds for all neighbours of central atom (with index 2)
            f.write('Bonds\n\n')
            f.write('\n'.join(f'{i_bond+1} 1 2 {i_bond+3}' for i_bond in range(n_cluster-1)))

            f.close()

    def write_all_transform_degrade(self, n_min=2, n_max=10):
        for n_cluster in range(n_min, n_max+1):
            self.write_single_chemistry(
                name='TransformMetabolites',
                n_cluster=n_cluster,
                type_cluster=1,
                start_type_single=2,
                end_type_single=3,
                delete_mode=False
            )
            self.write_single_chemistry(
                name='DegradeWaste',
                n_cluster=n_cluster,
                type_cluster=1,
                start_type_single=3,
                end_type_single=3,
                delete_mode=True
            )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--target-dir', type=str, default='reaction_templates/', help='Target directory for writing template files')
    parser.add_argument('-i', '--n-min', type=int, default=3, help='Minimum value for n_cluster')
    parser.add_argument('-a', '--n-max', type=int, default=11, help='Maximum value for n_cluster')
    
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    if not target_dir.is_dir():
        target_dir.mkdir()
    
    ChemistryWriter(target_dir).write_all_transform_degrade(n_min=args.n_min, n_max=args.n_max)
