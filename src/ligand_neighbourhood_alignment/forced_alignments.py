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
        chain: alignment_heirarchy._chain_to_biochain(chain, reference_structure_xtalform, assemblies)
        for chain in xtalform_chains
    }
    reference_structure_biochains_inv = {v: k for k, v in reference_structure_biochains.items()}

    # # Align the reference to the biochain reference using the canonical site residues
    residue_mapping = alignment_core.get_residue_mapping(  # ref-to-mov mapping
        reference_chain=canonical_site_biochain,
        moving_chain=reference_structure_biochains_inv[canonical_site_biochain],
        reference_structure=structures[site_reference_ligand_id[0]],
        moving_structure=reference_structure,
    )

    # Calculate the transform from the residue mapping
    transform = alignment_core.align_on_residues(
        residue_mapping=residue_mapping,
        reference_structure=structures[site_reference_ligand_id[0]],
        moving_structure=reference_structure,
    )

    reference_structure_transforms[key] = transform

    return transform
