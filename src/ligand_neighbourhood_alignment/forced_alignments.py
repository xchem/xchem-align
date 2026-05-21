from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment import alignment_heirarchy
from ligand_neighbourhood_alignment import alignment_core

# from ligand_neighbourhood_alignment.structure import _get_transform_from_residues


def _update_reference_structure_transforms(
    reference_structure_transforms,
    key,
    structures,
    canonical_site: dt.CanonicalSite,
    conformer_sites: dict[str, dt.ConformerSite],
    assemblies,
    xtalforms,
    dataset_assignments,
    xtalform_sites,
    canonical_site_id,
    debug=False,
):
    """
    Calculate the forced transforms that will superpose the relevant residues of a forced transform dataset
    to the relevant canconical site.

    Yes, the use of the word reference is unbelievably confusing here. Sorry. This should be fixed.
    """

    site_reference_ligand_id = conformer_sites[canonical_site.reference_conformer_site_id].reference_ligand_id
    site_reference_ligand_xtalform_id = dataset_assignments[site_reference_ligand_id[0]]
    site_reference_ligand_xtalform = xtalforms[site_reference_ligand_xtalform_id]

    # # Get the biochain of the canonical site
    canonical_site_biochain, site_reference_ligand_xtalform = alignment_heirarchy.get_canonical_site_biochain(
        site_reference_ligand_id,
        site_reference_ligand_xtalform_id,
        site_reference_ligand_xtalform,
        xtalform_sites,
        canonical_site_id,
        assemblies,
    )

    # Determine whether the biochain is shared, and if not skip
    reference_structure = structures[key[0]]
    reference_structure_xtalform = xtalforms[dataset_assignments[key[0]]]
    xtalform_chains = [
        chain for assembly in reference_structure_xtalform.assemblies.values() for chain in assembly.chains
    ]

    reference_structure_biochains = {
        chain: alignment_heirarchy._chain_to_biochains(chain, reference_structure_xtalform, assemblies)
        for chain in xtalform_chains
    }
    reference_structure_biochains_inv = {
        biochain: chain 
        for chain, biochains 
        in reference_structure_biochains.items()
        for biochain in biochains
        }

    canonical_site_reference_assemblies = set([xtalform_assembly.assembly for xtalform_assembly in site_reference_ligand_xtalform.assemblies.values()])
    reference_assemblies = set([xtalform_assembly.assembly for xtalform_assembly in reference_structure_xtalform.assemblies.values()])
    canonical_site_crystalographic_chain = alignment_heirarchy.get_canonical_site_crystallographic_chain(
            site_reference_ligand_id,
            site_reference_ligand_xtalform_id,
            xtalform_sites,
            canonical_site_id,
        )
    canonical_site_assemby = site_reference_ligand_xtalform.assemblies[
        alignment_heirarchy._chain_to_xtalform_assembly(
            canonical_site_crystalographic_chain, 
            site_reference_ligand_xtalform,
        )
        ].assembly
    
    # print(f'No assembly overlap between {canonical_site.reference_conformer_site_id} and {key}')
    # print(f'Reference assemblies: {reference_assemblies}')
    # print(f'Canonical site assembly: {canonical_site_assemby}')
    # print(f'Canonical site reference assemblies: {canonical_site_reference_assemblies}')
    # print(f'Site reference ligand xtalform id: {site_reference_ligand_xtalform_id}')
    # print(f'Reference structure xtalform: {dataset_assignments[key[0]]}')
    # print(f'Reference structure xtalform chains: {xtalform_chains}')
    # print(f'Reference structure biochains: {reference_structure_biochains}')
    # print(f'Reference structure biochains inverse: {reference_structure_biochains_inv}')


    # Skip if the assemblies aren't shared!
    if canonical_site_assemby not in reference_assemblies:
        print(f'No assembly overlap between {canonical_site.reference_conformer_site_id} and {key}')
        print(f'Reference assemblies: {reference_assemblies}')
        print(f'Canonical site assembly: {canonical_site_assemby}')
        print(f'Canonical site reference assemblies: {canonical_site_reference_assemblies}')
        print(f'Site reference ligand xtalform id: {site_reference_ligand_xtalform_id}')
        print(f'Reference structure xtalform: {dataset_assignments[key[0]]}')
        print(f'Reference structure xtalform chains: {xtalform_chains}')
        print(f'Reference structure biochains: {reference_structure_biochains}')
        print(f'Reference structure biochains inverse: {reference_structure_biochains_inv}')
        return

    # If the assemblies are shared, then get the chain in the reference assembly that maps to the
    # same biochain as the ligand binds to

    # # Align the reference to the biochain reference using the canonical site residues
    # try:
    residue_mapping = alignment_core.get_residue_mapping(  # ref-to-mov mapping
        reference_chain=canonical_site_biochain,
        moving_chain=reference_structure_biochains_inv[canonical_site_biochain],
        reference_structure=structures[site_reference_ligand_id[0]],
        moving_structure=reference_structure,
    )
    # except:
    #     print(f'Error in determining resiude mapping for:')
    #     print(f'{key}')
    #     print(f'Canonical site biochain: {canonical_site_biochain}')
    #     print(f'Biochain-to-reference chain: {reference_structure_biochains_inv}')
    #     print(f'Reference site ligand: {site_reference_ligand_id}')
    #     raise Exception
    
    # Calculate the transform from the residue mapping
    transform = alignment_core.align_on_residues(
        residue_mapping=residue_mapping,
        reference_structure=structures[site_reference_ligand_id[0]],
        moving_structure=reference_structure,
    )

    reference_structure_transforms[key] = transform

    return transform
