from rich import print as rprint
import gemmi
import numpy as np

from ligand_neighbourhood_alignment import dt
from src.ligand_neighbourhood_alignment import alignment_core
from ligand_neighbourhood_alignment.alignment_landmarks import structure_to_landmarks

def _derive_alignment_heirarchy(assemblies: dict[str, dt.Assembly], debug=False):
    # The Alignment hierarchy is the graph of alignments one must perform in order to get from
    # a ligand canonical site to the Reference Assembly Frame

    # In order to calculate the assembly the following steps are performed:
    # 1. Determine the Assembly priority
    # 2. Determine the Chain priority
    # 3. Find each assembly's reference
    # 4. Check per-chain RMSDs and warn if any are high

    # 1. Determine the Assembly priority
    assembly_priority = {_assembly_name: _j for _j, _assembly_name in enumerate(assemblies)}

    # 2. Determine the Chain priority and map assembly names to chains
    chain_priority = {}
    assembly_chains = {}
    chain_priority_count = 0
    for _assembly_name, _ in assembly_priority.items():
        assembly = assemblies[_assembly_name]
        assembly_chains[_assembly_name] = []
        for _generator in assembly.generators:
            _biological_chain_name = _generator.biomol
            assembly_chains[_assembly_name].append(_biological_chain_name)
            if _biological_chain_name not in chain_priority:
                chain_priority[_biological_chain_name] = chain_priority_count
                chain_priority_count += 1

    if debug:
        rprint(f'Assembly priority')
        rprint(assembly_priority)
        rprint(f'Chain priority')
        rprint(chain_priority)

    # 3. Find each assembly's reference
    reference_assemblies = {}
    for _assembly_name, _assembly in assemblies.items():
        # Get the highest priority chain
        reference_chain = min(
            [_generator.biomol for _generator in _assembly.generators], key=lambda _x: chain_priority[_x]
        )

        # Get the highest priority assembly in which it occurs
        reference_assembly = min(
            [_assembly_name for _assembly_name, chains in assembly_chains.items() if reference_chain in chains],
            key=lambda _x: assembly_priority[_x],
        )
        reference_assemblies[_assembly_name] = (reference_assembly, reference_chain)

    if debug:
        rprint(f'Reference Assemblies')
        rprint(reference_assemblies)

    # 4. Check per-chain RMSDs and warn if any are high
    # TODO

    return reference_assemblies, chain_priority


def _chain_to_biochain(chain_name, xtalform: dt.XtalForm, assemblies: dict[str, dt.Assembly]) -> str:
    for _xtal_assembly_name, _xtal_assembly in xtalform.assemblies.items():
        for _j, _chain_name in enumerate(_xtal_assembly.chains):
            if chain_name == _chain_name:
                return assemblies[_xtal_assembly.assembly].generators[_j].biomol
    raise Exception(f'No biochain found for chain {chain_name}')


def _get_assembly_st(as1, as1_ref):
    # Setup new structure to add biochains to
    new_st = gemmi.Structure()
    new_model = gemmi.Model("0")
    new_st.add_model(new_model)

    # Iterate over chain, biochain, transform tuples in the assembly
    for generator in as1.generators:
        # Generate the symmetry operation
        op = gemmi.Op(generator.triplet)

        # Create a clone of base chain to transform
        new_chain = as1_ref[0][generator.chain].clone()

        # Transform the residues
        for residue in new_chain:
            for atom in residue:
                atom_frac = as1_ref.cell.fractionalize(atom.pos)
                new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                new_pos_orth = as1_ref.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                atom.pos = gemmi.Position(*new_pos_orth)
        new_chain.name = generator.biomol
        new_st[0].add_chain(new_chain)

    new_st.add_model(new_model)
    return new_st





def _get_atoms(st):
    atoms = {}
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    atoms[(chain, res.name, res.seqid.num, atom.name)] = atom

    return atoms


def _calculate_assembly_transform(ref=None, mov=None, chain=None, debug=False):
    # Convert to gemmi structures to use superposition algorithm there
    ref_to_mov = alignment_core.calculate_insertion_matching_from_landmarks(ref, mov, chain, chain)  # Keys are ref atomids
    mov_insertion_to_res = {mov_atom_id[1][1]: mov_atom_id[1][0] for mov_atom_id in mov if mov_atom_id[0] == chain}

    ref_poss = []
    mov_poss = []
    for ref_atom_id in ref:
        # Only consider the chain
        if ref_atom_id[0] != chain:
            continue

        # Don't consider unmatched items
        if ref_atom_id[1][1] not in ref_to_mov:
            continue

        # Only align on CAs
        if ref_atom_id[2] != 'CA': 
            continue

        # Get the ref pos
        x, y, z = ref[ref_atom_id]
        pos = gemmi.Position(x, y, z)
        ref_poss.append(pos)

        # Determine the corresponding mov atom id
        mov_insertion = ref_to_mov[ref_atom_id[1][1]]
        mov_res = mov_insertion_to_res[mov_insertion]
        mov_atom_id = (ref_atom_id[0], (mov_res, mov_insertion), ref_atom_id[2])
        
        # Get the mov
        x, y, z = mov[mov_atom_id]
        pos = gemmi.Position(x, y, z)
        mov_poss.append(pos)


    assert (
        len(ref_poss) > 0
    ), f"""
    There are no valid reference positions to align. You may want to check residues numbers are the same between your assembly reference and datasets.
    {ref_to_mov}
    """
    assert (
        len(mov_poss) > 0
    ), f"""
    There are no valid dataset positions to align. You may want to check residues numbers are the same between your assembly reference and datasets.
    {ref_to_mov}
    """

    sup = gemmi.superpose_positions(ref_poss, mov_poss)
    transform = sup.transform

    # transform to interchangable format and return
    assert not np.isnan(np.array(transform.vec.tolist())).any()
    assert not np.isnan(np.array(transform.mat.tolist())).any()

    return {
        'vec': transform.vec.tolist(), 
        'mat': transform.mat.tolist(), 
        'rmsd': sup.rmsd, 
        "count": sup.count,
        'matching': ref_to_mov,
        }


def _calculate_assembly_sequence(
    hierarchy,
    mov_assembly,
):
    # Get the assembly sequence
    assembly_sequence = []
    running_assembly = mov_assembly
    while True:
        next_assembly, alignment_chain = hierarchy[running_assembly]
        # if len(assembly_sequence) != 0:
        if running_assembly == next_assembly:
            break
        else:
            assembly_sequence.append((next_assembly, alignment_chain))
        running_assembly = next_assembly

    return assembly_sequence


def _transform_to_gemmi(transform):
    tr = gemmi.Transform()
    tr.vec.fromlist(transform['vec'])
    tr.mat.fromlist(transform['mat'])
    return tr


def _calculate_assembly_transform_sequence(hierarchy, mov_assembly, assembly_landmarks, debug=False):
    # Get sequence of transforms
    sequence = _calculate_assembly_sequence(hierarchy, mov_assembly)

    # Get the sequence of transforms
    transforms = []
    moving = mov_assembly
    for assembly, chain in sequence:
        transform = _calculate_assembly_transform(
            ref=assembly_landmarks[assembly], 
            mov=assembly_landmarks[moving], 
            chain=chain
        )
        if debug:
            rprint(f'{moving} -> {assembly}')
            rprint(transform)
        moving = assembly
        transforms.append(_transform_to_gemmi(transform))

    # Sum the transforms
    tr = gemmi.Transform()
    tr.vec.fromlist([0.0, 0.0, 0.0])
    tr.mat.fromlist(np.eye(3).tolist())
    for transform in transforms:
        tr = transform.combine(tr)
    return {'vec': tr.vec.tolist(), 'mat': tr.mat.tolist(), 'rmsd': 0.0, "count": 0}


def _chain_to_xtalform_assembly(chain, xtalform):
    for assembly_name, assembly in xtalform.assemblies.items():
        if chain in [x for x in assembly.chains]:
            return assembly_name

    all_chains = {
        assembly_name: [x for x in assembly.chains] for assembly_name, assembly in xtalform.assemblies.items()
    }
    raise Exception(f'Chain {chain} not found in assembly chains: {all_chains}')


def _generate_assembly_from_xtalform(st, xtalform_assembly: dt.XtalFormAssembly, assembly, _debug=False):
    # Setup new structure to add biochains to
    new_st = gemmi.Structure()
    new_model = gemmi.Model("0")
    new_st.add_model(new_model)

    # Iterate over chain, biochain, transform tuples in the assembly
    for biomol, chain, transform in zip(
        [x.biomol for x in assembly.generators],
        xtalform_assembly.chains,
        xtalform_assembly.transforms,
    ):
        if _debug:
            rprint(f'Biomol: {biomol}')
            rprint(f'Chain: {chain}')
        # Generate the symmetry operation
        op = gemmi.Op(transform)

        # Create a clone of base chain to transform
        new_chain = st[0][chain].clone()

        # Transform the residues
        for residue in new_chain:
            for atom in residue:
                atom_frac = st.cell.fractionalize(atom.pos)
                new_pos_frac = op.apply_to_xyz([atom_frac.x, atom_frac.y, atom_frac.z])
                new_pos_orth = st.cell.orthogonalize(gemmi.Fractional(*new_pos_frac))
                atom.pos = gemmi.Position(*new_pos_orth)
        new_chain.name = biomol
        new_st[0].add_chain(new_chain)

    # new_st.add_model(new_model)
    return new_st


def _get_structure_chain_to_assembly_transform(st, chain, xtalform, assemblies, assembly_landmarks, debug=False):
    """
    Get the transform from a chain into a structure onto its biological assembly chain.

    This involves:
    1. Determine the assembly chain that corresponds to the structure chain
    2. Generating the biological assembly from THIS structure
    3. Transforming to landmarks for alignment
    4. Aligning the moving assembly generated from this chain against the canonical version of 
       that biological assembly
    """
    # Map the chain to an xtalform assembly
    xtalform_assembly_name = _chain_to_xtalform_assembly(chain, xtalform)
    xtalform_assembly = xtalform.assemblies[xtalform_assembly_name]
    if debug:
        rprint('Xtalform assembly for generating the structure to be aligned')
        rprint(xtalform_assembly_name)

    # Generate the assembly from the dataset structure
    assembly_st = _generate_assembly_from_xtalform(st, xtalform_assembly, assemblies[xtalform_assembly.assembly])
    assembly_st_chains = [chain.name for model in assembly_st for chain in model]
    if debug:
        rprint('Chains in assembly')
        rprint(assembly_st_chains)

    # Get the landmarks of the structure
    mov_lm = structure_to_landmarks(assembly_st)

    assert (
        len(mov_lm) != 0
    ), f"Should always have more than 0 move landmarks. Chain was {chain}, xtalform assembly name was {xtalform_assembly_name} and assembly chains are {assembly_st_chains}"

    # Align the structure assembly to the reference assembly
    if debug:
        rprint('Assembly landmarks')
        rprint(assembly_landmarks[xtalform_assembly.assembly])
        rprint('Moving landmarks')
        rprint(mov_lm)
    
    tr = _calculate_assembly_transform(
        ref=assembly_landmarks[xtalform_assembly.assembly],
        mov=mov_lm,
        chain=[x.biomol for x in assemblies[xtalform_assembly.assembly].generators][0],
    )

    if debug:
        rprint(tr)

    return tr




def chain_to_assembly_transforms_to_dict(chain_to_assembly_transforms: dict[tuple[str, str]]):
    dic = {}
    for k, v in chain_to_assembly_transforms.items():
        key = "~".join([k[0], k[1]])
        dic[key] = v

    return dic


def dict_to_chain_to_assembly_transforms(dic):
    obj = {}
    for k, v in dic.items():
        obj[(x for x in '~'.split(k))] = v

    return obj


def get_canonical_site_biochain(
    site_reference_ligand_id,
    site_reference_ligand_xtalform_id,
    site_reference_ligand_xtalform,
    xtalform_sites, 
    canonical_site_id, 
    assemblies,
    ):


    for xsid, _xtalform_site in xtalform_sites.items():
        _xtalform_id = _xtalform_site.xtalform_id
        if (
            (_xtalform_id == site_reference_ligand_xtalform_id)
            & (_xtalform_site.canonical_site_id == canonical_site_id)
            & (site_reference_ligand_id in _xtalform_site.members)
        ):
            xtalform_site = _xtalform_site
    site_chain = xtalform_site.crystallographic_chain
    canonical_site_biochain = _chain_to_biochain(
        site_chain, site_reference_ligand_xtalform, assemblies
    )

    return canonical_site_biochain, site_reference_ligand_xtalform
