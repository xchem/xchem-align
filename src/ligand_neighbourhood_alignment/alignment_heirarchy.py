from rich import print as rprint
import gemmi
import numpy as np
import yaml

from ligand_neighbourhood_alignment import dt, constants

AlignmentHeirarchy = dict[str, tuple[str, str]]


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


StructureLandmarks = dict[tuple[str, str, str], tuple[float, float, float]]


def structure_to_landmarks(st):
    landmarks = {}
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.name not in constants.RESIDUE_NAMES:
                    continue
                for atom in residue:
                    pos = atom.pos
                    landmarks[(chain.name, (residue.name, str(residue.seqid.num)), atom.name)] = (pos.x, pos.y, pos.z)

    return landmarks


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


def _landmark_to_structure(lm):
    st = gemmi.Structure()
    model = gemmi.Model("0")
    st.add_model(model)

    used_chains = []
    used_ress = []
    for (chain, (res_name, seqid), atom), (x, y, z) in lm.items():
        if chain not in used_chains:
            st[0].add_chain(gemmi.Chain(chain))
            used_chains.append(chain)

        if (chain, seqid) not in used_ress:
            new_residue = gemmi.Residue()
            new_residue.name = res_name
            new_residue.seqid = gemmi.SeqId(str(seqid))
            st[0][chain].add_residue(new_residue)
            used_ress.append((chain, seqid))

        new_atom = gemmi.Atom()
        new_atom.name = atom
        pos = gemmi.Position(x, y, z)
        new_atom.pos = pos
        st[0][chain][seqid][0].add_atom(new_atom)

    st.setup_entities()
    return st


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
    ref_matches = [atom_id for atom_id in ref if atom_id in mov]
    mov_matches = [atom_id for atom_id in mov if atom_id in ref]
    chain_matches_ref = [atom_id for atom_id in ref if atom_id[0] == chain]
    chain_matches_mov = [atom_id for atom_id in ref if atom_id[0] == chain]
    ref_poss = [
        gemmi.Position(x, y, z)
        for atom_id, (x, y, z) in ref.items()
        if (atom_id[0] == chain) & (atom_id in mov) & (atom_id[2] == 'CA')
    ]
    mov_poss = [
        gemmi.Position(x, y, z)
        for atom_id, (x, y, z) in mov.items()
        if (atom_id[0] == chain) & (atom_id in ref) & (atom_id[2] == 'CA')
    ]

    assert (
        len(ref_poss) > 0
    ), "There are no valid reference positions to align. You may want to check residues numbers are the same between your assembly reference and datasets."
    assert (
        len(mov_poss) > 0
    ), "There are no valid dataset positions to align. You may want to check residues numbers are the same between your assembly reference and datasets."

    sup = gemmi.superpose_positions(ref_poss, mov_poss)
    transform = sup.transform

    # transform to interchangable format and return
    assert not np.isnan(np.array(transform.vec.tolist())).any()
    assert not np.isnan(np.array(transform.mat.tolist())).any()

    return {'vec': transform.vec.tolist(), 'mat': transform.mat.tolist(), 'rmsd': sup.rmsd, "count": sup.count}


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
            ref=assembly_landmarks[assembly], mov=assembly_landmarks[moving], chain=chain
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


def _generate_assembly_from_xtalform(st, xtalform_assembly: dt.XtalFormAssembly, assembly):
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
    # Map the chain to an xtalform assembly
    xtalform_assembly_name = _chain_to_xtalform_assembly(chain, xtalform)
    xtalform_assembly = xtalform.assemblies[xtalform_assembly_name]

    # Generate the assembly from the dataset structure
    assembly_st = _generate_assembly_from_xtalform(st, xtalform_assembly, assemblies[xtalform_assembly.assembly])
    assembly_st_chains = [chain.name for model in assembly_st for chain in model]

    # Get the landmarks of the structure
    mov_lm = structure_to_landmarks(assembly_st)

    assert (
        len(mov_lm) != 0
    ), f"Should always have more than 0 move landmarks. Chain was {chain}, xtalform assembly name was {xtalform_assembly_name} and assembly chains are {assembly_st_chains}"

    # Align the structure assembly to the reference assembly
    tr = _calculate_assembly_transform(
        ref=assembly_landmarks[xtalform_assembly.assembly],
        mov=mov_lm,
        chain=[x.biomol for x in assemblies[xtalform_assembly.assembly].generators][0],
    )

    return tr


def save_yaml(path, obj, obj_to_dict):
    with open(path, "w") as f:
        dic = obj_to_dict(obj)
        yaml.safe_dump(dic, f, sort_keys=False)


def load_yaml(path, dict_to_obj):
    with open(path, "r") as f:
        dic = yaml.safe_load(f)

    return dict_to_obj(dic)


def assembly_landmarks_to_dict(assembly_landmarks: dict[tuple[str, tuple[str, str], str], tuple[float, float, float]]):
    dic = {}
    for assembly, data in assembly_landmarks.items():
        dic[assembly] = {}
        for k, v in data.items():
            key = "~".join([k[0], k[1][0], k[1][1], k[2]])
            dic[assembly][key] = v

    return dic


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


def dict_to_assembly_landmarks(dic):
    obj = {}
    for assembly, data in dic.items():
        obj[assembly] = {}
        for k, v in data.items():
            # key = [x for x in '~'.split(k)]
            chain, rname, rid, aname = [x for x in k.split('~')]
            obj[assembly][(chain, (rname, rid), aname)] = v

    return obj
