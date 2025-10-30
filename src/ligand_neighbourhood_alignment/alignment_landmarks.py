import gemmi

from ligand_neighbourhood_alignment import constants


def assembly_landmarks_to_dict(assembly_landmarks: dict[tuple[str, tuple[str, str], str], tuple[float, float, float]]):
    dic = {}
    for assembly, data in assembly_landmarks.items():
        dic[assembly] = {}
        for k, v in data.items():
            key = "~".join([k[0], k[1][0], k[1][1], k[2]])
            dic[assembly][key] = v

    return dic

def dict_to_assembly_landmarks(dic):
    obj = {}
    for assembly, data in dic.items():
        obj[assembly] = {}
        for k, v in data.items():
            # key = [x for x in '~'.split(k)]
            chain, rname, rid, aname = [x for x in k.split('~')]
            obj[assembly][(chain, (rname, rid), aname)] = v

    return obj

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