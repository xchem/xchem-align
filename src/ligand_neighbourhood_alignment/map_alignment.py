import itertools
from math import ceil, floor
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
from loguru import logger

from ligand_neighbourhood_alignment import alignment_heirarchy
from ligand_neighbourhood_alignment.dt import (
    Block,
    CanonicalSites,
    ConformerSites,
    Dataset,
    DatasetID,
    LigandBindingEvent,
    LigandID,
    LigandNeighbourhood,
    LigandNeighbourhoods,
    Output,
    SiteTransforms,
    SystemData,
    Transform,
    Transforms,
    gemmi_to_transform,
    read_xmap,
    transform_to_gemmi,
    write_xmap,
)

from ligand_neighbourhood_alignment import dt


def get_coord_array(neighbourhood: LigandNeighbourhood):
    coords = []
    for a in neighbourhood.atoms:
        coords.append([a.x, a.y, a.z])
    for a in neighbourhood.artefact_atoms:
        coords.append([a.x, a.y, a.z])

    return np.array(coords)


def _get_coord_array(neighbourhood: dt.Neighbourhood):
    coords = []
    for aid, a in neighbourhood.atoms.items():
        coords.append([a.x, a.y, a.z])
    for aid, a in neighbourhood.artefact_atoms.items():
        coords.append([a.x, a.y, a.z])

    return np.array(coords)


def get_bounds(coord_array):
    lb = np.min(coord_array, axis=0)
    up = np.max(coord_array, axis=0)

    return lb, up


def get_transformed_bounds(rlb, rub):
    arr = np.array([[rlb.x, rlb.y, rlb.z], [rub.x, rub.y, rub.z]])

    lb, ub = get_bounds(arr)
    return lb, ub


def get_grid_bounds(tlb, tub, xmap):
    cell = xmap.unit_cell

    tlbf = cell.fractionalize(gemmi.Position(*tlb))
    tubf = cell.fractionalize(gemmi.Position(*tub))

    tlbg = (
        floor(xmap.nu * tlbf.x),
        floor(xmap.nv * tlbf.y),
        floor(xmap.nw * tlbf.z),
    )
    tubg = (
        ceil(xmap.nu * tubf.x),
        ceil(xmap.nv * tubf.y),
        ceil(xmap.nw * tubf.z),
    )

    return tlbg, tubg


def get_subblocks(xrm, xr):
    xb = []
    xbt = []
    xit = xr[0]
    xbi = []
    for x, xm in zip(xr, xrm):
        if (xm == 0) & (len(xbt) > 0):
            xb.append(xbt)
            xbi.append(xit)
            xit = x
            xbt = []

        else:
            xbt.append(x)
    if len(xbt) > 0:
        xb.append(xbt)
        xbi.append(xit)
    return xb, xbi


def get_blocks(rglb, rgub, xmap):
    cell = xmap.unit_cell

    xr = np.arange(rglb[0], rgub[0] + 1)
    yr = np.arange(rglb[1], rgub[1] + 1)
    zr = np.arange(rglb[2], rgub[2] + 1)
    logger.debug(f"X range: {xr[0]} : {xr[-1]}")

    xrm = np.mod(xr, xmap.nu)
    yrm = np.mod(yr, xmap.nv)
    zrm = np.mod(zr, xmap.nw)
    logger.debug(f"X mod range: {xrm[0]} : {xrm[-1]}")

    xb, xbi = get_subblocks(xrm, xr)
    yb, ybi = get_subblocks(yrm, yr)
    zb, zbi = get_subblocks(zrm, zr)
    logger.debug(f"X subblock 0 range: {xb[0][0]} : {xb[0][-1]}")
    logger.debug(f"X subblock initial: {xbi[0]}")

    blocks = []
    for xsb, xsbi in zip(xb, xbi):
        for ysb, ysbi in zip(yb, ybi):
            for zsb, zsbi in zip(zb, zbi):
                # Transform jkl -> fractional -> orthogonal -> initial point
                # Fractional pos of block 0,0,0
                transform_vec = xmap.unit_cell.orthogonalize(
                    gemmi.Fractional(
                        xsbi / xmap.nu,
                        ysbi / xmap.nv,
                        zsbi / xmap.nw,
                    )
                )
                # Transform mat: fractional to orth
                orth_arr = np.array(cell.orth.mat.tolist())
                # Transform mat: grid to frac
                frac_arr = np.diag([1 / xmap.nu, 1 / xmap.nv, 1 / xmap.nw])

                block = Block(
                    xi=xsb[0],
                    yi=ysb[0],
                    zi=zsb[0],
                    xmi=np.mod(xsbi, xmap.nu),
                    ymi=np.mod(ysbi, xmap.nv),
                    zmi=np.mod(zsbi, xmap.nw),
                    dx=len(xsb),
                    dy=len(ysb),
                    dz=len(zsb),
                    transform=Transform(
                        vec=[
                            transform_vec.x,
                            transform_vec.y,
                            transform_vec.z,
                        ],
                        mat=(orth_arr @ frac_arr).tolist(),
                        alignable_ids=[]
                    ),
                )

                blocks.append(block)

    return blocks


def get_interpolation_range(neighbourhood: LigandNeighbourhood, transform, reference_xmap):
    # Get the gemmi transform
    transform_gemmi = transform

    # Coord array
    coord_arr = get_coord_array(neighbourhood)

    # Find the bounds
    lb, ub = get_bounds(coord_arr)
    logger.debug(f"Bounds of interpolation range in native frame are: {lb} : {ub}")

    # Transform the bounds
    corners = []
    for xb, yb, zb in itertools.product([lb[0], ub[0]], [lb[1], ub[1]], [lb[2], ub[2]]):
        corner = transform_gemmi.apply(gemmi.Position(xb, yb, zb))
        corners.append([corner.x, corner.y, corner.z])
    corner_array = np.array(corners)
    logger.debug(f"Transformed Corner array: {corner_array}")

    # Get the new bounds
    tlb, tub = get_bounds(corner_array)
    logger.debug(f"Bounds of interpolation range bounds in reference frame are: {tlb} : {tub}")

    # Get grid bounds
    rglb, rgub = get_grid_bounds(tlb, tub, reference_xmap)
    logger.debug(f"Grid bounds are: {rglb} : {rgub}")

    # Get the blocks
    blocks: list[Block] = get_blocks(rglb, rgub, reference_xmap)
    logger.debug(f"Num blocks: {blocks}")
    for b in blocks:
        s = f"Block: xi {b.xi} yi {b.yi} zi {b.zi} dx {b.dx} dy {b.dy} dz {b.dz}"
        logger.debug(s)
        logger.debug(b.transform)

    return blocks


def _get_interpolation_range(neighbourhood: dt.Neighbourhood, transform, reference_xmap):
    # Get the gemmi transform
    transform_gemmi = transform

    # Coord array
    coord_arr = _get_coord_array(neighbourhood)

    # Find the bounds
    lb, ub = get_bounds(coord_arr)
    logger.debug(f"Bounds of interpolation range in native frame are: {lb} : {ub}")

    # Transform the bounds
    corners = []
    for xb, yb, zb in itertools.product([lb[0], ub[0]], [lb[1], ub[1]], [lb[2], ub[2]]):
        corner = transform_gemmi.apply(gemmi.Position(xb, yb, zb))
        corners.append([corner.x, corner.y, corner.z])
    corner_array = np.array(corners)
    logger.debug(f"Transformed Corner array: {corner_array}")

    # Get the new bounds
    tlb, tub = get_bounds(corner_array)
    logger.debug(f"Bounds of interpolation range bounds in reference frame are: {tlb} : {tub}")

    # Get grid bounds
    rglb, rgub = get_grid_bounds(tlb, tub, reference_xmap)
    logger.debug(f"Grid bounds are: {rglb} : {rgub}")

    # Get the blocks
    blocks: list[Block] = get_blocks(rglb, rgub, reference_xmap)
    logger.debug(f"Num blocks: {blocks}")
    for b in blocks:
        s = f"Block: xi {b.xi} yi {b.yi} zi {b.zi} dx {b.dx} dy {b.dy} dz {b.dz}"
        logger.debug(s)
        logger.debug(b.transform)

    return blocks


def interpolate_range(
    reference_xmap,
    xmap,
    interpolation_ranges: list[Block],
    transform,
):
    # Make a xmap on reference template
    new_xmap = gemmi.FloatGrid(reference_xmap.nu, reference_xmap.nv, reference_xmap.nw)
    new_xmap.set_unit_cell(reference_xmap.unit_cell)
    new_xmap.spacegroup = gemmi.find_spacegroup_by_name("P1")
    grid_np = np.array(new_xmap, copy=False)
    logger.debug(f"Xmap shape is: {grid_np.shape}")

    # Interpolate values
    for interpolation_range in interpolation_ranges:
        arr = np.zeros(
            (
                interpolation_range.dx,
                interpolation_range.dy,
                interpolation_range.dz,
            ),
            dtype=np.float32,
        )
        logger.debug(f"Sample grid shape: {arr.shape}")

        range_transform = gemmi.Transform()
        range_transform.vec.fromlist(interpolation_range.transform.vec)
        range_transform.mat.fromlist(interpolation_range.transform.mat)

        xmap.interpolate_values(
            arr,
            transform.combine(transform_to_gemmi(interpolation_range.transform)),
        )

        # Assign to blocks
        rxi = interpolation_range.xmi
        rxf = interpolation_range.xmi + interpolation_range.dx
        ryi = interpolation_range.ymi
        ryf = interpolation_range.ymi + interpolation_range.dy
        rzi = interpolation_range.zmi
        rzf = interpolation_range.zmi + interpolation_range.dz
        logger.debug(f"Block X Range in output xmap: {rxi} : {rxf}")
        logger.debug(f"Block Y Range in output xmap: {ryi} : {ryf}")
        logger.debug(f"Block Z Range in output xmap: {rzi} : {rzf}")

        grid_np[
            rxi:rxf,
            ryi:ryf,
            rzi:rzf,
        ] = arr

    return new_xmap


def align_xmap(
    neighbourhoods: LigandNeighbourhoods,
    g,
    transforms: Transforms,
    site_transforms: SiteTransforms,
    reference_xmap,
    subsite_reference_id: LigandID,
    site_id: int,
    subsite_id: int,
    lid: LigandID,
    xmap,
    output_path: Path,
):
    # Get the ligand neighbourhood
    neighbourhood: LigandNeighbourhood = neighbourhoods.get_neighbourhood(lid)

    # Get the xmap

    # Get the Transform to reference
    running_transform = gemmi.Transform()
    shortest_path = nx.shortest_path(g, lid, subsite_reference_id)
    logger.debug(f"Shortest path: {shortest_path}")

    previous_ligand_id = lid
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = transforms.get_transform(
                (
                    next_ligand_id,
                    previous_ligand_id,
                ),
            )
            running_transform = transform.combine(running_transform)
        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    logger.debug(
        f"Transform from native frame to subsite frame to site frame is: {gemmi_to_transform(running_transform)}"
    )

    # Define the interpolation range
    interpolation_range = get_interpolation_range(neighbourhood, running_transform, reference_xmap)

    # Interpolate
    new_xmap = interpolate_range(
        reference_xmap,
        xmap,
        interpolation_range,
        running_transform.inverse(),
    )

    # Output the xmap
    write_xmap(
        new_xmap,
        output_path,
        neighbourhood,
        running_transform,
    )


def _get_box(neighbourhood: dt.Neighbourhood, xmap, transform):
    # transform_gemmi = transform_to_gemmi(transform)
    transform_gemmi = transform

    box = gemmi.FractionalBox()
    for atom_id, atom in neighbourhood.atoms.items():
        transformed_pos = transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
        box.extend(
            xmap.unit_cell.fractionalize(gemmi.Position(transformed_pos.x, transformed_pos.y, transformed_pos.z))
        )

    for atom_id, atom in neighbourhood.artefact_atoms.items():
        transformed_pos = transform_gemmi.apply(gemmi.Position(atom.x, atom.y, atom.z))
        box.extend(
            xmap.unit_cell.fractionalize(gemmi.Position(transformed_pos.x, transformed_pos.y, transformed_pos.z))
        )
    return box


def _write_xmap_from_ccp4(ccp4, path):
    ccp4.write_ccp4_map(str(path))


def _write_cut_xmap(new_xmap, aligned_res, crystallographic_output_path):
    new_st = gemmi.Structure()
    new_model = gemmi.Model("0")
    new_st.add_model(new_model)
    chain = gemmi.Chain("A")

    new_st.add_model(new_model)
    new_st[0].add_chain(chain)
    new_st[0][0].add_residue(aligned_res)

    new_st.cell = new_xmap.unit_cell
    new_st.spacegroup_hm = gemmi.find_spacegroup_by_name('P 1').hm

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = new_xmap
    ccp4.grid.spacegroup = gemmi.SpaceGroup("P1")
    ccp4.update_ccp4_header()
    ccp4.set_extent(new_st.calculate_fractional_box(margin=7.5))
    ccp4.write_ccp4_map(str(crystallographic_output_path))


def _write_xmap(xmap, path: Path, neighbourhood: dt.Neighbourhood, transform):
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = xmap
    ccp4.setup(float("nan"))
    ccp4.update_ccp4_header()

    if transform:
        box = _get_box(neighbourhood, xmap, transform)
        box_min = box.minimum
        box_max = box.maximum
        box_min_str = f"{round(box_min.x, 2)} {round(box_min.y, 2)} {round(box_min.z, 2)}"
        box_max_str = f"{round(box_max.x, 2)} {round(box_max.y, 2)} {round(box_max.z, 2)}"
        # print(f"Box Extent is: min {box_min_str} : max {box_max_str}")

        ccp4.set_extent(box)

    ccp4.setup(float("nan"))
    ccp4.update_ccp4_header()

    ccp4.write_ccp4_map(str(path))


def get_ligand_bounds(lig, border):
    poss = []
    for atom in lig:
        pos = atom.pos
        poss.append([pos.x, pos.y, pos.z])
    pos_array = np.array(poss)
    min_pos = np.min(pos_array, axis=0)
    max_pos = np.max(pos_array, axis=0)
    return min_pos, max_pos


def get_frame_bounds(ligand_lower_bound, ligand_upper_bound, border, step):
    min_border = ligand_lower_bound - border
    max_border = min_border + (np.ceil(((ligand_upper_bound + 5) - min_border) / step) * step)
    return min_border, max_border


def get_frame_array(frame_lower_bound, frame_upper_bound, step):
    interval = np.round((frame_upper_bound - frame_lower_bound) / step)
    return np.zeros((int(interval[0]), int(interval[1]), int(interval[2])), dtype=np.float32)


def get_frame_transform(frame_lower_bound, frame_array, step):
    tr = gemmi.Transform()
    tr.vec.fromlist([x for x in frame_lower_bound])
    tr.mat.fromlist((np.eye(3) * step).tolist())
    return tr


def get_cell(frame_array, step):
    shape = frame_array.shape
    cell = gemmi.UnitCell(shape[0] * step, shape[1] * step, shape[2] * step, 90.0, 90.0, 90.0)
    return cell


def get_new_map(cell, sample, frame_min, step):
    shape = sample.shape
    grid = gemmi.FloatGrid(shape[0], shape[1], shape[2])
    grid_array = np.array(grid, copy=False)
    grid_array[:, :, :] = sample[:, :, :]
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.grid.set_unit_cell(cell)
    ccp4.grid.spacegroup = gemmi.SpaceGroup("P1")
    ccp4.update_ccp4_header()
    ccp4.set_header_float(50, frame_min[0])
    ccp4.set_header_float(51, frame_min[1])
    ccp4.set_header_float(52, frame_min[2])
    return ccp4


def resample_xmap(new_xmap, aligned_res):
    step = 0.5
    border = 7.5
    m = new_xmap
    lig = aligned_res

    ligand_lower_bound, ligand_upper_bound = get_ligand_bounds(lig, border)
    frame_lower_bound, frame_upper_bound = get_frame_bounds(ligand_lower_bound, ligand_upper_bound, border, step)

    frame_array = get_frame_array(frame_lower_bound, frame_upper_bound, step)

    frame_transform = get_frame_transform(frame_lower_bound, frame_array, step)
    m.interpolate_values(frame_array, frame_transform)
    cell = get_cell(frame_array, step)
    new_map = get_new_map(cell, frame_array, frame_lower_bound, step)
    return new_map


def __align_xmap(
    neighbourhood: dt.Neighbourhood,
    g,
    ligand_neighbourhood_transforms: dict[tuple[tuple[str, str, str], tuple[str, str, str]], dt.Transform],
    reference_xmap,
    subsite_reference_id: tuple[str, str, str],
    lid: tuple[str, str, str],
    xmap,
    output_path: Path,
    crystallographic_output_path,
    aligned_res,
    chain_to_assembly_transform,
    assembly_transform,
):
    # Get the xmap

    # Get the Transform to reference
    running_transform = gemmi.Transform()
    shortest_path = nx.shortest_path(g, lid, subsite_reference_id)
    logger.debug(f"SHORTEST PATH: {lid}: {shortest_path}")

    previous_ligand_id = lid
    for next_ligand_id in shortest_path:
        # Get the transform from previous frame to new one
        # Transform is 2 onto 1
        if next_ligand_id != previous_ligand_id:
            transform = ligand_neighbourhood_transforms[
                (
                    next_ligand_id,
                    previous_ligand_id,
                )
            ]
            running_transform = transform_to_gemmi(transform).combine(running_transform)
        # Apply the translation to the new frame
        previous_ligand_id = next_ligand_id

    # Update the transform with the assembly alignment
    # # Get the xtalform to assembly transform
    running_transform = alignment_heirarchy._transform_to_gemmi(chain_to_assembly_transform).combine(running_transform)

    # # Get the assembly alignment transform
    running_transform = alignment_heirarchy._transform_to_gemmi(assembly_transform).combine(running_transform)

    logger.debug(
        f"Transform from native frame to subsite frame to site frame is: {gemmi_to_transform(running_transform)}"
    )

    # Define the interpolation range
    interpolation_range = _get_interpolation_range(neighbourhood, running_transform, reference_xmap)

    # Interpolate
    new_xmap = interpolate_range(
        reference_xmap,
        xmap,
        interpolation_range,
        running_transform.inverse(),
    )

    # Resample the xmap to the aligned structure frame
    resampled_xmap = resample_xmap(new_xmap, aligned_res)

    _write_xmap_from_ccp4(
        resampled_xmap,
        output_path,
    )

    # Write a restricted crystallographic map
    _write_cut_xmap(new_xmap, aligned_res, crystallographic_output_path)


def read_xmap_from_mtz(
    mtz_path: Path,
    map_type="2Fo-Fc",
):
    mtz = gemmi.read_mtz_file(str(mtz_path))

    if map_type == "2Fo-Fc":
        try:
            grid = mtz.transform_f_phi_to_map("2FOFCWT", "PH2FOFCWT", sample_rate=4)
            return grid
        except Exception:
            logger.info("Trying FWT PHWT")
        try:
            grid = mtz.transform_f_phi_to_map("FWT", "PHWT", sample_rate=4)
            return grid
        except Exception:
            logger.warning("Failed to find structure factors!")

        raise Exception(f"Couldn't open {mtz_path}!")
    if map_type == "Fo-Fc":
        try:
            grid = mtz.transform_f_phi_to_map("FOFCWT", "PHFOFCWT", sample_rate=4)
            return grid
        except Exception:
            logger.info("Trying DELFWT DELPHWT")
        try:
            grid = mtz.transform_f_phi_to_map("DELFWT", "PHDELWT", sample_rate=4)
            return grid
        except Exception:
            logger.warning("Failed to find structure factors!")

        raise Exception(f"Couldn't open {mtz_path}!")


def _align_xmaps(
    system_data: SystemData,
    structures,
    canonical_sites: CanonicalSites,
    conformer_sites: ConformerSites,
    neighbourhoods: LigandNeighbourhoods,
    g,
    transforms: Transforms,
    site_transforms: SiteTransforms,
    output: Output,
):
    for canonical_site_id, canonical_site in canonical_sites.iter():
        logger.debug(f"Aligning site: {canonical_site_id}")

        reference_lid: LigandID = canonical_site.reference_ligand_id
        referance_ds: Dataset = system_data.get_dataset(DatasetID(dtag=reference_lid.dtag))
        logger.debug(f"PDB: {referance_ds.pdb}")
        reference_mtz_path = Path(referance_ds.mtz)
        reference_xmap = read_xmap_from_mtz(reference_mtz_path)

        for conformer_site_id, conformer_site in conformer_sites.iter():
            if conformer_site_id not in canonical_site.subsite_ids:
                continue
            logger.debug(f"Aligning subsite: {conformer_site_id}")
            # Get the site reference
            # TODO: Make work
            conformer_site_reference_id = conformer_site.reference_ligand_id

            # For each ligand neighbourhood, find the xmap and transform
            for lid in conformer_site.members:
                logger.debug(f"Aligning xmap: {lid}")

                dtag, chain, residue = (
                    lid.dtag,
                    lid.chain,
                    lid.residue,
                )

                logger.debug(f"Reference pdb: {referance_ds.pdb}")
                logger.debug(f"Moving pdb: {system_data.get_dataset(DatasetID(dtag=lid.dtag)).pdb}")

                # Get the ligand binding event
                dataset = system_data.get_dataset(DatasetID(dtag=lid.dtag))
                lbe: LigandBindingEvent = dataset.ligand_binding_events[lid]

                # Get the xmap path
                #

                if lbe.xmap != "None":
                    xmap_path: Path = Path(lbe.xmap)
                    logger.debug(f"Xmap path: {xmap_path}")

                    xmap = read_xmap(xmap_path)
                else:
                    mtz_path: Path = Path(dataset.mtz)
                    logger.debug(f"Mtz path is: {mtz_path}")

                    xmap = read_xmap_from_mtz(mtz_path)

                aop = Path(output.source_dir)

                # Align the event map
                output_path = aop / output.dataset_output[dtag][chain][residue].aligned_event_maps[canonical_site_id]
                logger.debug("Aligning event map...")
                align_xmap(
                    neighbourhoods,
                    g,
                    transforms,
                    site_transforms,
                    reference_xmap,
                    conformer_site_reference_id,
                    canonical_site_id,
                    conformer_site_id,
                    lid,
                    xmap,
                    output_path,
                )
                # Align the refined map
                output_path = aop / output.dataset_output[dtag][chain][residue].aligned_xmaps[canonical_site_id]

                if dataset.mtz:
                    logger.debug("Aligning 2Fo-Fc map...")

                    xmap = read_xmap_from_mtz(Path(dataset.mtz))

                    align_xmap(
                        neighbourhoods,
                        g,
                        transforms,
                        site_transforms,
                        reference_xmap,
                        conformer_site_reference_id,
                        canonical_site_id,
                        conformer_site_id,
                        lid,
                        xmap,
                        output_path,
                    )


# def _align_xmap(
#     system_data: SystemData,
#     canonical_sites: CanonicalSites,
#     conformer_sites: ConformerSites,
#     neighbourhoods: LigandNeighbourhoods,
#     g,
#     transforms: Transforms,
#     site_transforms: SiteTransforms,
#     output: Output,
# ):
#     reference_lid: LigandID = canonical_site.reference_ligand_id
#     referance_ds: Dataset = system_data.get_dataset(DatasetID(dtag=reference_lid.dtag))
#     logger.debug(f"PDB: {referance_ds.pdb}")
#     reference_mtz_path = Path(referance_ds.mtz)
#     reference_xmap = read_xmap_from_mtz(reference_mtz_path)

#     # Get the site reference
#     # TODO: Make work
#     conformer_site_reference_id = conformer_site.reference_ligand_id

#     dtag, chain, residue = (
#         lid.dtag,
#         lid.chain,
#         lid.residue,
#     )
#     # if dtag != "Mpro-J0055":
#     #     continue

#     logger.debug(f"Reference pdb: {referance_ds.pdb}")

#     # Get the ligand binding event
#     dataset = system_data.get_dataset(DatasetID(dtag=lid.dtag))
#     lbe: LigandBindingEvent = dataset.ligand_binding_events[lid]

#     # Get the xmap path
#     #

#     if lbe.xmap != "None":
#         xmap_path: Path = Path(lbe.xmap)
#         logger.debug(f"Xmap path: {xmap_path}")

#         xmap = read_xmap(xmap_path)
#     else:
#         mtz_path: Path = Path(dataset.mtz)
#         logger.debug(f"Mtz path is: {mtz_path}")

#         xmap = read_xmap_from_mtz(mtz_path)

#     aop = Path(output.source_dir)

#     # Align the event map
#     output_path = aop / output.dataset_output[dtag][chain][residue].aligned_event_maps[canonical_site_id]
#     logger.debug("Aligning event map...")
#     align_xmap(
#         neighbourhoods,
#         g,
#         transforms,
#         site_transforms,
#         reference_xmap,
#         conformer_site_reference_id,
#         canonical_site_id,
#         conformer_site_id,
#         lid,
#         xmap,
#         output_path,
#     )
#     # Align the refined map
#     output_path = aop / output.dataset_output[dtag][chain][residue].aligned_xmaps[canonical_site_id]

#     if dataset.mtz:
#         logger.debug("Aligning 2Fo-Fc map...")

#         xmap = read_xmap_from_mtz(Path(dataset.mtz))

#         align_xmap(
#             neighbourhoods,
#             g,
#             transforms,
#             site_transforms,
#             reference_xmap,
#             conformer_site_reference_id,
#             canonical_site_id,
#             conformer_site_id,
#             lid,
#             xmap,
#             output_path,
#         )
