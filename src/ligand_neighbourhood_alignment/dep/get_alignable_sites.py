from ligand_neighbourhood_alignment.data import AlignableSite


def match_site(containing_site, contained_site):
    if len([ligand_id for ligand_id in contained_site if ligand_id in containing_site]) == len(contained_site):
        return True
    else:
        return False


# def get_alignable_sites_()


def get_alignable_sites(connected_components, alignable_sites: list[AlignableSite] | None):
    alignable_site_num: int = 0
    if alignable_sites:
        alignable_site_num = len(alignable_sites)
    else:
        alignable_sites = []

    for site_ligand_ids in connected_components:
        site_match: bool = False
        for alignable_site in alignable_sites:
            site_match = match_site(
                containing_site=site_ligand_ids,
                contained_site=alignable_site.ligand_ids,
            )

            # Update site if matched
            if site_match:
                for _ligand_id in site_ligand_ids:
                    if _ligand_id not in alignable_site.ligand_ids:
                        alignable_site.ligand_ids.append(_ligand_id)

        # If not site has matched create a new alignable site
        if not site_match:
            alignable_sites.append(AlignableSite(id=alignable_site_num, name="", ligand_ids=site_ligand_ids))
            alignable_site_num += 1

    return alignable_sites
