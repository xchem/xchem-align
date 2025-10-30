import numpy as np

from ligand_neighbourhood_alignment import dt
from ligand_neighbourhood_alignment.ligand_neighbourhoods import _get_centroid_res

def _get_dist(pos_1, pos_2):
    return np.linalg.norm(np.array([pos_1.x, pos_1.y, pos_1.z]) - np.array([pos_2.x, pos_2.y, pos_2.z]))


def _crystalform_incremental_cluster(observation_centroid_residues, xtalform_sites, neighbourhoods, cutoff=10.0):
    # Check if each observation centroid is within 10A of an xtalform site centroid, and if so assign that dataset to it
    # If not, select the remainer with the most members within 10A and form a new "site" from it
    # Continue until all observations are assigned

    # Get CA positions
    centroid_ca_positions = {
        observation_id: neighbourhoods[observation_id].atoms[(centroid_res[0], centroid_res[1], 'CA')]
        for observation_id, centroid_res in observation_centroid_residues.items()
    }

    # Identify current xtalform centre residues (for the considered canonical site)
    centre_residues_positions = {
        xtalform_site_id: centroid_ca_positions[tuple([dt.string_to_altloc(x) for x in xtalform_site_id.split('/')])]
        for xtalform_site_id in xtalform_sites
    }

    # While observations remain to be assigned
    assignments = {
        xtalform_site_id: [_x for _x in xtalform_sites[xtalform_site_id].members]
        for xtalform_site_id in xtalform_sites
    }
    assigned_observations = [
        observation_id for xtalform_site_id in assignments for observation_id in assignments[xtalform_site_id]
    ]
    observations_to_assign = [observation_id for observation_id in centroid_ca_positions]

    while len(observations_to_assign) > 0:
        for observation_id in observations_to_assign:
            # Assign observations near current xtalform sites
            for xtalform_site_id, xtalform_site_pos in centre_residues_positions.items():
                if _get_dist(centroid_ca_positions[observation_id], xtalform_site_pos) < cutoff:
                    assignments[xtalform_site_id].append(observation_id)
                    assigned_observations.append(observation_id)
                    break

        # Identify observation with most remaining items within 10A
        remaining_observations = [_k for _k in observations_to_assign if _k not in assigned_observations]
        if len(remaining_observations) == 0:
            break
        num_near = {
            observation_id: len(
                [
                    other_observation_id
                    for other_observation_id in remaining_observations
                    if (
                        _get_dist(centroid_ca_positions[observation_id], centroid_ca_positions[other_observation_id])
                        < cutoff
                    )
                    & (observation_id != other_observation_id)
                ]
            )
            for observation_id in remaining_observations
        }
        new_centroid_observation_id = min(num_near, key=lambda _x: num_near[_x])

        # Add it to centre residue list
        centre_residues_positions[new_centroid_observation_id] = centroid_ca_positions[new_centroid_observation_id]
        assignments[new_centroid_observation_id] = []

        # Update observations to assign
        # assigned_observations.append(new_centroid_observation_id)
        observations_to_assign = [_k for _k in observations_to_assign if _k not in assigned_observations]

    return assignments





def _update_xtalform_sites(
    xtalform_sites: dict[str, dt.XtalFormSite],
    canonical_site: dt.CanonicalSite,
    canonical_site_id: str,
    dataset_assignments: dict[str, str],
    conformer_sites: dict[str, dt.ConformerSite],
    neighbourhoods,
    debug=False,
):
    # Iterate over Canonical Sites, collecting Observations that are in the same crystalform
    # Then get their centroid CA positions,
    # spatially cluster on these with a reasonably broad cutoff, including old items
    # Then name these collections according to their centremost member

    # Partition by canonical site and xtalform
    crystalform_observations = {
        xtalform_name: [
            member
            for conformer_site_id in canonical_site.conformer_site_ids
            for member in conformer_sites[conformer_site_id].members
            if dataset_assignments[member[0]] == xtalform_name
        ]
        for xtalform_name in set(dataset_assignments.values())
    }

    # Get Observation centroid CA names and positions
    crystalform_observation_centroids = {
        xtalform_name: {
            member: _get_centroid_res(
                [_x for _x in set([(aid[0], aid[1]) for aid in neighbourhoods[member].atoms])], neighbourhoods[member]
            )
            for member in crystalform_observations[xtalform_name]
        }
        for xtalform_name in crystalform_observations
    }

    # Spatially cluster
    crystalform_observation_cluster_assignments = {
        xtalform_name: _crystalform_incremental_cluster(
            crystalform_observation_centroids[xtalform_name],
            {
                xid: xs
                for xid, xs in xtalform_sites.items()
                if (xs.xtalform_id == xtalform_name) & (xs.canonical_site_id == canonical_site_id)
            },
            neighbourhoods,
        )
        for xtalform_name in crystalform_observations
    }

    # if debug:
    #     raise Exception

    # Create the xtalforms or assign new observations
    for xtalform_name in crystalform_observation_cluster_assignments:
        for centroid_residue, asigned_observation_ids in crystalform_observation_cluster_assignments[
            xtalform_name
        ].items():
            # If the centroid is known, assign any new observations
            if centroid_residue in xtalform_sites:
                for asigned_observation_id in asigned_observation_ids:
                    if asigned_observation_id not in xtalform_sites[centroid_residue].members:
                        xtalform_sites[centroid_residue].members.append(asigned_observation_id)

            # Otherwise create a new crystalform site
            else:
                xtalform_site_id = "/".join(centroid_residue)
                xtalform_site = dt.XtalFormSite(
                    xtalform_name,
                    crystalform_observation_centroids[xtalform_name][centroid_residue][0],
                    canonical_site_id,
                    asigned_observation_ids,
                )
                xtalform_sites[xtalform_site_id] = xtalform_site

    # Otherwise if not matched create a new xtalform site
    ...
