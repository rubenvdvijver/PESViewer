import rdkit
from rdkit import Chem


def num_rad_elecs(molec: rdkit.Chem.Mol) -> int:
    """Count the total number of radical electrons in a molecule.

    @param molec: rdkit Mol object of the target molecule.
    """
    return sum([a.GetNumRadicalElectrons() for a in molec.GetAtoms()])


def gen_bond_combs(n_bonds: int, n_conns: int) -> list:
    """Generates all possible combinations of single double and triple bonds.

    @param n_bonds: Number of shared pairs of electrons: A double bond counts
        as two.
    @param n_conns: Number of connections between atoms: A double bond counts
        as one.

    The only condition to fulfill at this step is that the sum of each of
    the combinations has to be equal to n_bonds. This accepts configurations
    breaking the octet.
    """
    import itertools
    iters = [[1, 2, 3]] * n_conns
    combs = []
    for comb in itertools.product(*iters):
        if sum(comb) == n_bonds:
            combs.append(comb)
    return combs


def filter_valid_structs(mol: rdkit.Chem.Mol, combs: list, hvy_bond_ids: list) -> list:
    """Removes all configurations exceeding the octet rule, radicals are allowed.

    @param mol: rdkit Mol object of the target molecule.
    @param combs: Possible combinations of bond configuration for the molecule.
    @param hvy_bond_ids: List of bonds connecting non-H atoms (C-C, C-O, etc.)
    """
    from copy import deepcopy
    from rdkit.Chem.rdchem import BondType
    from rdkit import rdBase, RDLogger
    # Disable logging
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.ERROR)
    rdBase.DisableLog('rdApp.error')

    atomic_valences = {'C': 4, 'N': 3, 'O': 2, 'H': 1}
    valid_mols = []
    for comb in combs:
        new_mol = deepcopy(mol)
        for i, bo in enumerate(comb):
            bond = new_mol.GetBondWithIdx(hvy_bond_ids[i])
            bond.SetBondType(BondType(bo))
        try:
            Chem.SanitizeMol(new_mol, Chem.SANITIZE_ALL)
        except:  # TODO Catch exact exception if possible (which I think is not)
            continue

        # Correct radical centers
        for a in new_mol.GetAtoms():
            if a.GetExplicitValence() < atomic_valences[a.GetSymbol()]:
                n_rad_elecs = atomic_valences[a.GetSymbol()] \
                                - a.GetExplicitValence()
            else:
                n_rad_elecs = 0
            a.SetNumRadicalElectrons(n_rad_elecs)
        if not any([Chem.MolToSmiles(new_mol) == Chem.MolToSmiles(vmol)
                    for vmol in valid_mols]):
            valid_mols.append(new_mol)
    return valid_mols


def gen_reso_structs(smi: str, min_rads=True) -> list:  # C(=C\\1/[C]C1)\\[CH2]
    """Generate all possible resonant structures of a given molecule.

    @param smi: SMILES code of the given molecule
    @param min_rads: Whether to minimize the number of radical electrons or not.
    """
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    hvy_bond_ids = [b.GetIdx() for b in mol.GetBonds()
                    if b.GetBeginAtom().GetSymbol() != 'H'
                    and b.GetEndAtom().GetSymbol() != 'H']
    num_bonds = int(sum([mol.GetBondWithIdx(i).GetBondTypeAsDouble()
                         for i in hvy_bond_ids]))
    num_conns = len(hvy_bond_ids)
    radic_elecs = num_rad_elecs(mol)
    max_bonds = num_bonds + radic_elecs // 2
    valid_structs = []
    for new_bonds in range(num_bonds, max_bonds + 1):
        combs = gen_bond_combs(new_bonds, num_conns)
        valid_structs += filter_valid_structs(mol, combs, hvy_bond_ids)

    # Remove Structures with larger number of radical electrons
    if min_rads:
        min_num_rads = min([num_rad_elecs(struct) for struct in valid_structs])
        valid_structs = [struct for struct in valid_structs
                         if num_rad_elecs(struct) == min_num_rads]

    # Remove explicit Hydrogen atoms
    for s, struct in enumerate(valid_structs):
        valid_structs[s] = Chem.RemoveHs(struct)

    return valid_structs


if __name__ == '__main__':
    import sys
    from rdkit.Chem import Draw
    smiles = sys.argv[1]
    valid_strs = gen_reso_structs(smiles)
    for j, m in enumerate(valid_strs):
        Draw.MolToFile(m, f'aaa{j}.png')
