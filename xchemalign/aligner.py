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

import argparse
from pathlib import Path

from . import utils, processor

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


class Aligner(processor.Processor):

    def validate(self):
        v = validator.Validator(self.base_dir, self.input_dirs, self.output_dir, self.target_name, logger=self.logger)
        warnings, errors = v.validate_paths()
        return warnings, errors

    def run(self):
        self.logger.info('Running aligner...')
        v_dir = self.read_versions()
        if not v_dir:
            self.logger.error('Error with version dir. Please fix and try again.')
            return None
        self.logger.info('Using version dir {}'.format(v_dir))

        meta = self.read_metadata(v_dir)

        self._perform_alignments(self.config, meta, self.meta_history, [], self.output_dir)

    def _perform_alignments(self, config, meta, meta_history, panddas_dirs, output_dir):
        self.logger.info('Performing alignments (well, not really)')

        # Initialize the output directory and create empty
        # jsons in it

        # Add the datasources in the options json and add them to
        # the datasource json
        # visits = meta[constants.META_INPUT]
        crystals = meta[constants.META_CRYSTALS]

        # datasources = [
        #     Datasource(
        #         path=visit[constants.META_DATA_DIR],
        #         datasource_type=visit[constants.META_DATA_DIR_TYPE],
        #     )
        #     for visit
        #     in visits
        # ]
        # panddas = [
        #     PanDDA(
        #         path=visit[constants.META_PANDDAS][constants.META_PANDDAS_PATH],
        #         event_table_path=visit[constants.META_PANDDAS][constants.META_PANDDAS_EVENT_TABLE_PATH],
        #     )
        #     for visit
        #     in visits
        # ]
        dataset_ids = [DatasetID(dtag=dtag) for dtag in crystals]
        datasets = [
            Dataset(
                dtag=dtag,
                pdb=crystal[constants.META_XTAL_PDB][constants.META_XTAL_PDB_FILE],
                xmap=None,
                mtz=crystal[constants.META_XTAL_MTZ][constants.META_XTAL_MTZ_FILE],
                ligand_binding_events=LigandBindingEvents(
                    ligand_ids=[
                        LigandID(
                            dtag=dtag,
                            chain=binding_event[constants.META_BINDING_EVENT_CHAIN],
                            residue=binding_event[constants.META_BINDING_EVENT_RES],
                        )
                        for binding_event
                        in crystal[constants.META_BINDING_EVENT]
                    ],
                    ligand_binding_events=[
                        LigandBindingEvent(
                            id=0,
                            dtag=dtag,
                            chain=binding_event[constants.META_BINDING_EVENT_CHAIN],
                            residue=binding_event[constants.META_BINDING_EVENT_RES],
                            xmap=binding_event[constants.META_BINDING_EVENT_EVENT_MAP],
                        )
                        for binding_event
                        in crystal[constants.META_BINDING_EVENT]
                    ]
                ),
            )
            for dtag, crystal
            in crystals.items()
        ]
        system_data = SystemData(
            # datasources=datasources,
            # panddas=panddas,
            datasources=[],
            panddas=[],
            dataset_ids=dataset_ids,
            datasets=datasets
        )

        # Copy the xtalform json into the source directory (checking validity)
        xtalforms = XtalForms.read(Path(meta[constants.XTALFORM_JSON]))

        # Parse the data sources and PanDDAs, matching ligands up to events
        # system_data = _add_data_to_system_data(system_data)

        # Assign each dataset to the clsoest xtalform and fail if this
        # is not possible
        assigned_xtalforms = _get_assigned_xtalforms(system_data, xtalforms)

        # Build the alignment graph
        ligand_neighbourhoods: LigandNeighbourhoods = get_ligand_neighbourhoods(
            system_data,
            xtalforms,
            assigned_xtalforms,
        )

        # Get alignability
        alignability_matrix, transforms = get_alignability(ligand_neighbourhoods, system_data)

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
        output = read_output(Path(meta[constants.ALIGNED_DIR]))
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
            aligned_output = crystal[constants.META_ALIGNED] = {}

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
                            constants.META_ALIGNED_STRUCTURE: aligned_structure_path,
                            constants.META_ALIGNED_ARTEFACTS: aligned_artefacts_path,
                            constants.META_ALIGNED_EVENT_MAP: aligned_event_map_path,
                            constants.META_ALIGNED_XMAP: aligned_xmap_path
                        }

        return meta


def main():
    parser = argparse.ArgumentParser(description='aligner')

    parser.add_argument('-c', '--config-file', default='config.yaml', help="Configuration file")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("aligner: ", args)

    logger = utils.Logger(logfile=args.log_file, level=args.log_level)

    a = Aligner(args.config_file, logger=logger)

    warnings, errors = a.validate()

    if not args.validate:
        if errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            a.run()


if __name__ == "__main__":
    main()
