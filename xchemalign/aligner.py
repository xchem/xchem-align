# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse, yaml
from pathlib import Path

from . import utils
from .utils import Constants

# Local alignment imports
from ligand_neighbourhood_alignment import constants as lna_constants
from ligand_neighbourhood_alignment.align_xmaps import _align_xmaps
from ligand_neighbourhood_alignment.data import (
    Assemblies,
    AssignedXtalForms,
    CanonicalSites,
    ChainOutput,
    ConformerSites,
    Dataset,
    DatasetID,
    DatasetOutput,
    Datasource,
    LigandBindingEvent,
    LigandBindingEvents,
    LigandID,
    LigandNeighbourhoods,
    LigandOutput,
    Options,
    Output,
    PanDDA,
    SiteTransforms,
    SystemData,
    Transforms,
    XtalForms,
    read_assigned_xtalforms,
    read_canonical_sites,
    read_graph,
    read_neighbourhoods,
    read_output,
    read_site_transforms,
    read_structures,
    read_system_data,
    read_transforms,
    read_xtalforms,
    save_assigned_xtalforms,
    save_canonical_sites,
    save_data,
    save_output,
)
from ligand_neighbourhood_alignment.generate_aligned_structures import _align_structures_from_sites
from ligand_neighbourhood_alignment.generate_sites_from_components import (
    _generate_sites_from_components,
    get_components,
    get_conformer_sites_from_components,
    get_site_transforms,
    get_sites_from_conformer_sites,
    get_structures,
    get_subsite_transforms,
)
from ligand_neighbourhood_alignment.get_alignability import get_alignability
from ligand_neighbourhood_alignment.get_graph import get_graph
from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from ligand_neighbourhood_alignment.cli import _add_model_building_dir_to_system_data, _add_manual_dir_to_system_data, \
    _add_pandda_to_system_data, _add_data_to_system_data, _get_assigned_xtalforms


class Aligner():

    def __init__(self, version_dir, metadata_file, xtalforms_file, logger=None):
        self.version_dir = Path(version_dir)
        self.aligned_dir = self.version_dir / Constants.META_ALIGNED_FILES
        self.metadata_file = metadata_file
        self.xtalforms_file = xtalforms_file
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()
        self.errors = []
        self.warnings = []

    def _log_error(self, msg):
        self.logger.error(msg)
        self.errors.append(msg)

    def validate(self):

        if not self.version_dir.exists():
            self._log_error('version dir {} does not exist'.format(self.version_dir))
        elif not self.version_dir.is_dir():
            self._log_error('version dir {} is not a directory'.format(self.version_dir))
        else:
            p = self.version_dir / self.metadata_file
            if not p.exists():
                self._log_error('metadata_file {} does not exist'.format(p))
            if not p.is_file():
                self._log_error('metadata_file {} is not a file'.format(p))

        p = Path(self.xtalforms_file)
        if not p.exists():
            self._log_error('xtalforms_file {} does not exist'.format(p))
        elif not p.is_file():
            self._log_error('xtalforms_file {} is not a file'.format(p))

        return len(self.errors), len(self.warnings)

    def run(self):
        self.logger.info('Running aligner...')

        meta_path = self.version_dir / self.metadata_file
        meta = utils.read_config_file(str(meta_path))

        if not self.aligned_dir.is_dir():
            self.aligned_dir.mkdir()
            self.logger.info('created aligned directory', self.aligned_dir)

        new_meta = self._perform_alignments(meta)
        # TODO - should aligned metadata be written to its own file?
        with open(meta_path, 'w') as stream:
            yaml.dump(new_meta, stream, sort_keys=False)

    def _perform_alignments(self, meta):
        self.logger.info('Performing alignments')

        # Initialize the output directory and create empty
        # jsons in it

        # Add the datasources in the options json and add them to
        # the datasource json
        # visits = meta[lna_constants.META_INPUT]
        crystals = meta[Constants.META_XTALS]
        output_path = Path(meta[Constants.META_OUTPUT_DIR])

        # Assert that
        if len(crystals) == 0:
            self.logger.error(f"Did not find any crystals in metadata. Exiting.")
            raise Exception

        dataset_ids = [DatasetID(dtag=dtag) for dtag in crystals]
        datasets = [
            Dataset(
                dtag=dtag,
                pdb=str(output_path / crystal[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB][Constants.META_FILE]),
                xmap='',
                mtz=str(output_path / crystal[Constants.META_XTAL_FILES].get(Constants.META_XTAL_MTZ, {}).get(Constants.META_FILE)),
                ligand_binding_events=LigandBindingEvents(
                    ligand_ids=[
                        LigandID(
                            dtag=dtag,
                            chain=binding_event.get(Constants.META_PROT_CHAIN),
                            residue=binding_event.get(Constants.META_PROT_RES),
                        )
                        for binding_event
                        in crystal[Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})
                    ],
                    ligand_binding_events=[
                        LigandBindingEvent(
                            id=0,
                            dtag=dtag,
                            chain=binding_event.get(Constants.META_PROT_CHAIN),
                            residue=binding_event.get(Constants.META_PROT_RES),
                            xmap=str(output_path / binding_event.get(Constants.META_FILE)),
                        )
                        for binding_event
                        in crystal[Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})
                    ]
                ),
            )
            for dtag, crystal
            in crystals.items()
        ]
        if len(dataset_ids) != len(datasets):
            self.logger.error(f"Number of dataset ids found in metadata not equal to number of datasets. Exiting.")
            raise Exception

        if (len(datasets) == 0) or (len(datasets) == 0):
            self.logger.error(f"Did not find any datasets in metadata. Exiting.")
            raise Exception

        self.logger.info(f"Ligand binding events in datasets:")
        for _dataset_id, _dataset in zip(dataset_ids, datasets):
            _num_binding_events = len(_dataset.ligand_binding_events.ligand_binding_events)
            self.logger.info(f"\t{_dataset_id.dtag} : Num Ligand binding events: {_num_binding_events}")

        system_data = SystemData(
            datasources=[],
            panddas=[],
            dataset_ids=dataset_ids,
            datasets=datasets
        )

        # Copy the xtalform json into the source directory (checking validity)
        # xtalforms = XtalForms.read(Path(meta[lna_constants.XTALFORM_JSON]))
        xtalforms = XtalForms.read(self.xtalforms_file)

        # Parse the data sources and PanDDAs, matching ligands up to events
        # system_data = _add_data_to_system_data(system_data)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)
        self.logger.info(f"Assigned xtalforms are:")
        for _dataset_id, _assigned_xtalform in assigned_xtalforms:
            self.logger.info(f"\t{_dataset_id} : {_assigned_xtalform}")

        # Build the alignment graph
        ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
            system_data,
            xtalforms,
            assigned_xtalforms,
        )
        self.logger.info(f"Found {len(ligand_neighbourhoods.ligand_ids)} ligand neighbourhoods.")
        self.logger.info(f"Ligand neighbourhoods are:")
        for _ligand_id, _ligand_neighbourhood in zip(ligand_neighbourhoods.ligand_ids, ligand_neighbourhoods.ligand_neighbourhoods):
            _dtag, _chain, _residue = _ligand_id.dtag, _ligand_id.chain, _ligand_id.residue
            _num_atoms, _num_art_atoms = len(_ligand_neighbourhood.atoms), len(_ligand_neighbourhood.artefact_atoms)
            self.logger.info(f"\t{_dtag} {_chain} {_residue} : Num atoms: {_num_atoms} : Num artefact atoms: {_num_art_atoms}")

        # Get alignability
        alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)
        _x, _y = alignability_matrix.shape
        _z = len(ligand_neighbourhoods.ligand_ids)
        if (_x != _y) or (_x != _z) or (_y != _z):
            self.logger.error(f"Alignability matrix should be of shape: {_z} x {_z}, however is {_x} x {_y}")
            raise Exception

        # Generate the graph
        g = get_graph(alignability_matrix, ligand_neighbourhoods)

        # Generate canonical, conformer and xtalform sites from the
        # alignment graph

        # Get the connected components
        connected_components = get_components(g)

        # Get the subsites from the connected components with overlap
        conformer_sites: ConformerSites = get_conformer_sites_from_components(
            connected_components, ligand_neighbourhoods
        )

        # Merge the connected components with shared residues into sites
        _sites = get_sites_from_conformer_sites(conformer_sites, ligand_neighbourhoods)
        if len(_sites) == 0:
            self.logger.error(f"Number of sites is 0: this should be impossible. Exiting.")
            raise Exception

        canonical_sites: CanonicalSites = CanonicalSites(
            site_ids=[s.id for s in _sites],
            sites=_sites,
            reference_site=_sites[0],
            reference_site_id=_sites[0].id,
        )

        # Get the xtalform sites
        # xtalform_sites = get_xtalform_sites_from_canonical_sites(
        #     canonical_sites,
        #     assigned_xtalforms,
        #     xtalforms,
        #     # assemblies,
        # )

        # Get the subsite transforms
        structures = get_structures(system_data)
        subsite_transforms = get_subsite_transforms(canonical_sites, structures)

        # Get the site transforms
        site_transforms = get_site_transforms(canonical_sites, structures)
        site_transforms = SiteTransforms(
            canonical_site_transform_ids=[key for key in site_transforms.keys()],
            canonical_site_transforms=[tr for tr in site_transforms.values()],
            conformer_site_transform_ids=[key for key in subsite_transforms.keys()],
            conformer_site_transforms=[tr for tr in subsite_transforms.values()],
        )

        # Fully specify the output now that the sites are known
        output = Output(source_dir=str(self.version_dir),
        system_data='',
        xtalforms='',
        assigned_xtalforms='',
        neighbourhoods='',
        graph='',
        transforms='',
        sites='',
        site_transforms='',
        aligned_dir=Constants.META_ALIGNED_FILES,
        dataset_output={})

        dataset_output_dict = {}
        for ligand_id in ligand_neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            if dtag not in dataset_output_dict:
                dataset_output = DatasetOutput(aligned_chain_output={})
                dataset_output_dict[dtag] = dataset_output
            else:
                dataset_output = dataset_output_dict[dtag]

            if chain not in dataset_output.aligned_chain_output:
                chain_output = ChainOutput(
                    aligned_ligands={},
                )
                dataset_output_dict[dtag].aligned_chain_output[chain] = chain_output
            else:
                chain_output = dataset_output_dict[dtag].aligned_chain_output[chain]

            chain_output.aligned_ligands[residue] = LigandOutput(
                aligned_structures={}, aligned_artefacts={}, aligned_xmaps={}, aligned_event_maps={}
            )

            # Add output for each canonical site that the ligand is aligned to
            for site_id, site in canonical_sites.iter():
                if ligand_id not in site.members:
                    continue

                chain_output.aligned_ligands[residue].aligned_structures[site_id] = (
                        output.aligned_dir
                        + "/"
                        + lna_constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                        output.aligned_dir
                        + "/"
                        + lna_constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + lna_constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue, site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                        output.aligned_dir
                        + "/"
                        + lna_constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

        # Save the output file
        output.dataset_output = dataset_output_dict

        # Align structures to each canonical site
        _align_structures_from_sites(
            structures,
            canonical_sites,
            conformer_sites,
            transforms,
            ligand_neighbourhoods,
            xtalforms,
            assigned_xtalforms,
            g,
            site_transforms,
            output,
        )

        # Align xmaps to each canonical site
        _align_xmaps(
            system_data,
            structures,
            canonical_sites,
            conformer_sites,
            ligand_neighbourhoods,
            g,
            transforms,
            site_transforms,
            output,
        )

        # Update the metadata with aligned file locations and site information
        for dtag, crystal in crystals.items():

            # Skip if no output for this dataset
            if dtag not in dataset_output_dict:
                continue

            # Otherwise iterate the output data structure, adding the aligned structure,
            # artefacts, xmaps and event maps to the metadata
            aligned_output = crystal[Constants.META_ALIGNED_FILES] = {}
            dataset_output = dataset_output_dict[dtag]
            for chain_name, chain_output in dataset_output.aligned_chain_output.items():
                aligned_chain_output = aligned_output[chain_name] = {}
                for ligand_residue, ligand_output in chain_output.aligned_ligands.items():
                    aligned_ligand_output = aligned_chain_output[ligand_residue] = {}
                    for site_id, aligned_structure_path in ligand_output.aligned_structures.items():
                        aligned_artefacts_path = ligand_output.aligned_artefacts[site_id]
                        aligned_event_map_path = ligand_output.aligned_event_maps[site_id]
                        aligned_xmap_path = ligand_output.aligned_xmaps[site_id]
                        aligned_ligand_output[site_id] = {
                            Constants.META_AIGNED_STRUCTURE: aligned_structure_path,
                            Constants.META_AIGNED_ARTEFACTS: aligned_artefacts_path,
                            Constants.META_AIGNED_EVENT_MAP: aligned_event_map_path,
                            Constants.META_AIGNED_X_MAP: aligned_xmap_path
                        }

        return meta


def main():
    parser = argparse.ArgumentParser(description='aligner')

    parser.add_argument('-v', '--version-dir', required=True, help="Path to version dir")
    parser.add_argument('-d', '--metadata', default='metadata.yaml', help="Metadata YAML file")
    parser.add_argument('-x', '--xtalforms', default='xtalforms.json', help="Crystal forms JSON file")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("aligner: ", args)

    logger = utils.Logger(logfile=args.log_file, level=args.log_level)

    a = Aligner(args.version_dir, args.metadata, args.xtalforms, logger=logger)
    num_errors, num_warnings = a.validate()

    if not args.validate:
        if num_errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            a.run()


if __name__ == "__main__":
    main()
