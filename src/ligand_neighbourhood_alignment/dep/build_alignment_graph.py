from pathlib import Path

from loguru import logger

from ligand_neighbourhood_alignment import constants
from ligand_neighbourhood_alignment.data import (
    LigandNeighbourhoods,
    SystemData,
    read_assigned_xtalforms,
    read_xtalforms,
)
from ligand_neighbourhood_alignment.get_alignability import get_alignability
from ligand_neighbourhood_alignment.get_graph import get_graph
from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from ligand_neighbourhood_alignment.save_graph import save_graph
from ligand_neighbourhood_alignment.save_neighbourhoods import save_neighbourhoods
from ligand_neighbourhood_alignment.save_transforms import save_transforms


def build_alignment_graph(output_dir: Path):
    # logger.info(f"System sites json path is: {system_sites_json_path}")

    data_json_path: Path = output_dir / constants.DATA_JSON_PATH
    logger.info(f"Data json path is: {data_json_path}")

    # # Load the input data
    system_data: SystemData = SystemData.parse_file(str(data_json_path))

    num_input_datasets = len(system_data.datasets)

    logger.debug(f"Got {num_input_datasets} datasets")

    xtalforms = read_xtalforms(output_dir)
    assigned_xtalforms = read_assigned_xtalforms(output_dir)

    # Identify sites in the input data
    ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
        system_data,
        xtalforms,
        assigned_xtalforms,
    )

    num_neighbourhoods = len(ligand_neighbourhoods.ligand_neighbourhoods)
    logger.info(f"Found {num_neighbourhoods} ligand neighbourhoods")

    # Save the neighbourhoods
    logger.info("Saving neighbourhoods!")

    save_neighbourhoods(ligand_neighbourhoods, output_dir / "neighbourhoods.json")
    logger.info("Saved neighbourhoods!")

    # Get alignability
    logger.info("Getting alignbaility matrix...!")
    alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)
    logger.info("Got alignability matrix!")

    # logger.debug(alignability_matrix)
    logger.debug("Alignability matrix shape: {alignability_matrix.shape}")

    # Generate the graph
    logger.info("Getting alignability graph...")
    g = get_graph(alignability_matrix, ligand_neighbourhoods)
    logger.info("Got alignability graph!")

    # Write the graph
    logger.info("Saving Graph...")
    save_graph(g, output_dir / "alignability.gml")
    logger.info("Saved Graph!")

    # Generate the transforms
    logger.info("Getting transforms...")
    # transforms = get_transforms(ligand_neighbourhoods, g)
    logger.info("Got transforms!")

    # Save the transforms
    logger.info("Saving transforms...")
    save_transforms(transforms, output_dir / "transforms.json")
    logger.info("Saved transforms!")
