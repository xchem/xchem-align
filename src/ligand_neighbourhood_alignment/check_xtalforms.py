import numpy as np

def _check_xtalforms(xtalforms, structures):
    """
    For every crystalform, check every other crystalform, and raise an exception if the unit cell is within tolderance
    and it shares the same spacegroup
    """
    identical_crystalforms = []
    all_deltas = {}
    for ref_xtalform_id, ref_xtalform in xtalforms.items():
        ref_structure = structures[ref_xtalform.reference]
        ref_spacegroup = ref_structure.spacegroup_hm
        ref_structure_cell = ref_structure.cell
        for mov_xtalform_id, mov_xtalform in xtalforms.items():
            mov_structure = structures[mov_xtalform.reference]
            mov_spacegroup = mov_structure.spacegroup_hm
            mov_structure_cell = mov_structure.cell

            if ref_xtalform_id == mov_xtalform_id:
                continue

            if mov_spacegroup != ref_spacegroup:
                continue

            # Check deltas
            deltas = np.round(np.array(
            [
                mov_structure_cell.a / ref_structure_cell.a,
                mov_structure_cell.b / ref_structure_cell.b,
                mov_structure_cell.c / ref_structure_cell.c,
                mov_structure_cell.alpha / ref_structure_cell.alpha,
                mov_structure_cell.beta / ref_structure_cell.beta,
                mov_structure_cell.gamma / ref_structure_cell.gamma,
            ]
            ), 2)
            if np.any(deltas > 1.1) | np.any(deltas < 0.91):
                continue

            # If this code executes then the spacegroup is the same and the unit cell is close: raise 
            # an Exception!
            identical_crystalform = set([ref_xtalform_id, mov_xtalform_id])
            if identical_crystalform not in identical_crystalforms:
                identical_crystalforms.append(identical_crystalform)
            all_deltas[identical_crystalform] = deltas
        
    if len(identical_crystalforms) != 0:
        exception = f'There were identical crystalforms present in the assemblies.yaml \n'
        for identical_crystalform in identical_crystalforms:
            identical_crystalform_tuple = tuple(identical_crystalform)
            print(identical_crystalform_tuple)
            exception += f'\t{identical_crystalform_tuple[0]} is identical to {identical_crystalform_tuple[1]}: deltas {all_deltas[identical_crystalform]}\n'
        exception += 'Please remove the redundancy and rerun!'
        raise Exception(
                exception
            )