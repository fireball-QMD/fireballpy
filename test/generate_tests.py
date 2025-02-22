import json
from os.path import split, join
from ase import Atoms
from ase.io import read
from fireballpy.fdata import load_fdata
from fireballpy import Fireball


FDATA = 'CHON'
CHARGES_METHOD = ['lowdin', 'weighted_lowdin', 'mulliken',
                  'mulliken_dipole', 'mulliken_dipole_preserving']


def generate_test(filename: str, kpts: tuple[int, int, int]) -> dict[str, float | list[float]]:
    xyz = join(split(__file__)[0], filename)
    atoms = read(xyz)
    assert isinstance(atoms, Atoms)

    dipoles = ['legacy']
    if kpts == (1,1,1):
        dipoles += ['improved']

    results = {}
    for cm in CHARGES_METHOD:
        for dp in dipoles:
            key = f'{dp}::{cm}'
            print(key)
            atoms.calc = Fireball(fdata=FDATA, kpts=kpts, lazy=False, charges_method=cm,
                                  dipole_method=dp, correction=False, verbose=True)

            results[key] = {}
            results[key]['potential_energy'] = float(atoms.get_potential_energy())
            results[key]['partial_charges'] = atoms.get_charges().tolist()
            results[key]['shell_charges'] = atoms.calc.get_shell_charges().tolist()
            results[key]['eigenvalues'] = atoms.calc.get_eigenvalues().tolist()
            results[key]['forces'] = atoms.get_forces().tolist()
    return results


def save_test(filename: str, results: dict[str, float | list[float]]) -> None:
    with open(join(split(__file__)[0], filename), 'w', encoding='utf-8') as fp:
        json.dump(results, fp, ensure_ascii=False, indent=4)


if __name__ == '__main__':
    load_fdata(FDATA)

    # Molecule
    resmolecule = generate_test('molecule.xyz', kpts=(1,1,1))
    save_test('molecule.json', resmolecule)

    # Periodic
    resperiodic = generate_test('periodic.xyz', kpts=(4,4,4))
    save_test('periodic.json', resperiodic)
