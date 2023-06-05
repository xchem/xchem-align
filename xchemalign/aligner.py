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
import json
import os
from pathlib import Path

import gemmi
import yaml

# Local alignment imports
from ligand_neighbourhood_alignment import constants as lna_constants
from ligand_neighbourhood_alignment.align_xmaps import _align_xmaps
from ligand_neighbourhood_alignment.data import (
    CanonicalSites,
    ChainOutput,
    ConformerSites,
    Dataset,
    DatasetID,
    DatasetOutput,
    LigandBindingEvent,
    LigandBindingEvents,
    LigandID,
    LigandNeighbourhoods,
    LigandOutput,
    Output,
    SiteTransforms,
    SystemData,
    XtalForms,
)
from ligand_neighbourhood_alignment.generate_aligned_structures import _align_structures_from_sites
from ligand_neighbourhood_alignment.generate_sites_from_components import (
    get_components,
    get_conformer_sites_from_components,
    get_site_transforms,
    get_sites_from_conformer_sites,
    get_structures,
    get_subsite_transforms,
    get_xtalform_sites_from_canonical_sites,
)
from ligand_neighbourhood_alignment.get_alignability import get_alignability
from ligand_neighbourhood_alignment.get_graph import get_graph
from ligand_neighbourhood_alignment.get_ligand_neighbourhoods import get_ligand_neighbourhoods
from ligand_neighbourhood_alignment.cli import _get_assigned_xtalforms

from . import utils
from .utils import Constants

def try_make(path):
    if not Path(path).exists():
        os.mkdir(path)

class Aligner():

    def __init__(self, version_dir, metadata, xtalforms, logger=None):
        self.version_dir = Path(version_dir)                                # e.g. path/to/upload_1
        self.base_dir = self.version_dir.parent                             # e.g. path/to
        self.aligned_dir = self.version_dir / Constants.META_ALIGNED_FILES  # e.g. path/to/upload_1/aligned_files
        self.xtal_dir = self.version_dir / Constants.META_XTAL_FILES        # e.g. path/to/upload_1/crystallographic_files
        self.metadata_file = self.version_dir / metadata                    # e.g. path/to/upload_1/metadata.yaml
        self.xtalforms_file = self.base_dir.parent / xtalforms              # e.g. path/to/xtalforms.yaml
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
            p = self.metadata_file
            if not p.exists():
                self._log_error('metadata file {} does not exist'.format(p))
            if not p.is_file():
                self._log_error('metadata file {} is not a file'.format(p))

            p = Path(self.xtalforms_file)
            if not p.exists():
                self._log_error('xtalforms file {} does not exist'.format(p))
            elif not p.is_file():
                self._log_error('xtalforms file {} is not a file'.format(p))

        return len(self.errors), len(self.warnings)

    def run(self):
        self.logger.info('Running aligner...')

        meta = utils.read_config_file(str(self.metadata_file))

        if not self.aligned_dir.is_dir():
            self.aligned_dir.mkdir()
            self.logger.info('created aligned directory', self.aligned_dir)

        new_meta = self._perform_alignments(meta)

        # TODO - should aligned metadata_file be written to its own file?
        with open(self.version_dir / Constants.METADATA_ALIGN_FILENAME, 'w') as stream:
            yaml.dump(new_meta, stream, sort_keys=False, default_flow_style=None)

    def _perform_alignments(self, meta):
        self.logger.info('Performing alignments')

        # Initialize the output directory and create empty
        # jsons in it

        # Add the datasources in the options json and add them to
        # the datasource json
        # visits = meta[lna_constants.META_INPUT]
        crystals = meta[Constants.META_XTALS]
        output_path = Path(meta[Constants.CONFIG_OUTPUT_DIR])

        # Assert that
        if len(crystals) == 0:
            self.logger.error(f"Did not find any crystals in metadata file. Exiting.")
            raise Exception

        dataset_ids = [DatasetID(dtag=dtag) for dtag in crystals]
        # paths to files will be defined like this: upload_1/crystallographic_files/8dz1/8dz1.pdb
        # this is relative to the output_path variable that is defined from the metadata_file.yaml
        datasets = [
            Dataset(
                dtag=dtag,
                pdb=str(output_path / crystal[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB][Constants.META_FILE]),
                xmap='',
                mtz=str(output_path / crystal[Constants.META_XTAL_FILES].get(Constants.META_XTAL_MTZ, {}).get(
                    Constants.META_FILE)),
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
            self.logger.error(f"Number of dataset ids found in metadata_file not equal to number of datasets. Exiting.")
            raise Exception

        if (len(datasets) == 0) or (len(datasets) == 0):
            self.logger.error(f"Did not find any datasets in metadata_file. Exiting.")
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

        # Convert the xtalform yaml into a json and read into an Xtalforms object
        # xtalforms = XtalForms.read(Path(meta[lna_constants.XTALFORM_JSON]))
        with open(self.xtalforms_file, "r") as f:
            xtalform_dict = yaml.safe_load(f)
        xtalform_path = self.version_dir / 'xtalforms.json'
        with open(xtalform_path, "w") as f:
            json.dump(xtalform_dict, f)
        xtalforms = XtalForms.read(xtalform_path)

        # Parse the data sources and PanDDAs, matching ligands up to events
        # system_data = _add_data_to_system_data(system_data)

        # Assign each dataset to the closest xtalform and fail if this
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
        for _ligand_id, _ligand_neighbourhood in zip(ligand_neighbourhoods.ligand_ids,
                                                     ligand_neighbourhoods.ligand_neighbourhoods):
            _dtag, _chain, _residue = _ligand_id.dtag, _ligand_id.chain, _ligand_id.residue
            _num_atoms, _num_art_atoms = len(_ligand_neighbourhood.atoms), len(_ligand_neighbourhood.artefact_atoms)
            self.logger.info(
                f"\t{_dtag} {_chain} {_residue} : Num atoms: {_num_atoms} : Num artefact atoms: {_num_art_atoms}")

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
        xtalform_sites = get_xtalform_sites_from_canonical_sites(
            canonical_sites,
            assigned_xtalforms,
            xtalforms,
            # assemblies,
        )

        # Get the subsite transforms
        structures = get_structures(system_data)
        subsite_transforms = get_subsite_transforms(canonical_sites, structures)

        # Get the site transforms
        site_transforms = get_site_transforms(canonical_sites, structures)
        site_transforms = SiteTransforms(
            # pylint: disable=consider-iterating-dictionary
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
                        aligned_dir=str(self.aligned_dir),
                        dataset_output={})

        # Create the aligned dir
        try_make(output.aligned_dir)

        dataset_output_dict = {}
        for ligand_id in ligand_neighbourhoods.ligand_ids:
            dtag, chain, residue = (
                ligand_id.dtag,
                ligand_id.chain,
                ligand_id.residue,
            )

            # Create output dataset dir if not already exists
            dataset_output_dir = output.aligned_dir + "/" + f"{dtag}"
            try_make(dataset_output_dir)

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
                        Constants.META_ALIGNED_FILES + '/' + dtag + '/' +
                        lna_constants.ALIGNED_STRUCTURE_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_artefacts[site_id] = (
                        Constants.META_ALIGNED_FILES + '/' + dtag + '/' +
                        lna_constants.ALIGNED_STRUCTURE_ARTEFACTS_TEMPLATE.format(
                    dtag=dtag, chain=chain, residue=residue, site=site_id
                )
                )

                chain_output.aligned_ligands[residue].aligned_xmaps[site_id] = (
                        Constants.META_ALIGNED_FILES + '/' + dtag + '/' +
                        lna_constants.ALIGNED_XMAP_TEMPLATE.format(dtag=dtag, chain=chain, residue=residue,
                                                                     site=site_id)
                )

                chain_output.aligned_ligands[residue].aligned_event_maps[site_id] = (
                        Constants.META_ALIGNED_FILES + '/' + dtag + '/' +
                        lna_constants.ALIGNED_EVENT_MAP_TEMPLATE.format(
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

        # Update the metadata_file with aligned file locations and site information
        new_meta = {}

        # Add the xtalform information
        meta_xtalforms = {}
        meta_assemblies = {}
        for xtalform_id, xtalform in xtalforms.xtalforms.items():
            xtalform_reference = xtalform.reference
            reference_structure = gemmi.read_structure(system_data.get_dataset(xtalform_reference).pdb)
            reference_spacegroup = reference_structure.spacegroup_hm
            reference_unit_cell = reference_structure.cell

            meta_xtalforms[xtalform_id] = {
                Constants.META_XTALFORM_REFERENCE: xtalform_reference.dtag,
                Constants.META_XTALFORM_SPACEGROUP: reference_spacegroup,
                Constants.META_XTALFORM_CELL: {
                    "a":reference_unit_cell.a,
                    "b": reference_unit_cell.b,
                    "c": reference_unit_cell.c,
                    "alpha": reference_unit_cell.alpha,
                    "beta": reference_unit_cell.beta,
                    "gamma": reference_unit_cell.gamma,
                },
            }
            for assembly_id, assembly in xtalform.assemblies.items():

                assembly_ref_chains = []
                for generator_id, generator in assembly.generators.items():
                    ref_chain, chain, triplet = generator.reference_chain, generator.chain, generator.triplet
                    assembly_ref_chains.append(ref_chain)

                assembly_ref_chains_tup = tuple(assembly_ref_chains)

                # Create an assembly or add one
                if assembly_ref_chains_tup not in meta_assemblies:
                    meta_assemblies[assembly_ref_chains_tup] = {
                        Constants.META_ASSEMBLIES_XTALFORMS: [xtalform_id,]
                    }
                else:
                    meta_assemblies[assembly_ref_chains_tup][Constants.META_ASSEMBLIES_XTALFORMS].append(xtalform_id)

        new_meta[Constants.META_XTALFORMS] = meta_xtalforms
        new_meta[Constants.META_ASSEMBLIES] = meta_assemblies

        print(conformer_sites)
        print("##############################")
        print(canonical_sites)
        print("##############################")
        print(xtalform_sites)

        # Add the conformer sites
        conformer_sites_meta = new_meta[Constants.META_CONFORMER_SITES] = {}
        for conformer_site_id, conformer_site in conformer_sites.conformer_sites.items():
            conformer_sites_meta[conformer_site_id] = {
                Constants.META_CONFORMER_SITE_NAME: None,
                Constants.META_CONFORMER_SITE_REFERENCE_LIG: {

                    Constants.META_DTAG: conformer_site.reference_ligand_id.dtag,
                    Constants.META_CHAIN: conformer_site.reference_ligand_id.chain,
                    Constants.META_RESIDUE: conformer_site.reference_ligand_id.residue,
                },
                Constants.META_CONFORMER_SITE_RESIDUES: {

                        Constants.META_CHAIN: [res.chain for res in conformer_site.residues],
                        Constants.META_RESIDUE: [res.residue for res in conformer_site.residues]

                    },
                Constants.META_CONFORMER_SITE_MEMBERS: {
                        Constants.META_DTAG: [lid.dtag for lid in conformer_site.members],
                     Constants.META_CHAIN: [lid.chain for lid in conformer_site.members],
                     Constants.META_RESIDUE: [lid.residue for lid in conformer_site.members]
                     }

            }

        # Add the canonical sites
        canonical_sites_meta = new_meta[Constants.META_CANONICAL_SITES] = {}
        for canonical_site_id, canonical_site in zip(canonical_sites.site_ids, canonical_sites.sites):
            canonical_sites_meta[canonical_site_id]= {
                Constants.META_CANONICAL_SITE_REF_SUBSITE: canonical_site.reference_subsite_id,
                Constants.META_CANONICAL_SITE_CONFORMER_SITES: canonical_site.subsite_ids,
                Constants.META_CANONICAL_SITE_RESIDUES: {

                        Constants.META_CHAIN: [res.chain for res in canonical_site.residues],
                        Constants.META_RESIDUE: [res.residue for res in canonical_site.residues]

                    },
                Constants.META_CANONICAL_SITE_MEMBERS: {
                        Constants.META_DTAG: [lid.dtag for lid in canonical_site.members],
                     Constants.META_CHAIN: [lid.chain for lid in canonical_site.members],
                     Constants.META_RESIDUE: [lid.residue for lid in canonical_site.members]
                     }
            }

        # Add the xtalform sites - note the chain is that of the original crystal structure, NOT the assembly
        xtalform_sites_meta = new_meta[Constants.META_XTALFORM_SITES] = {}
        for xtalform_site_id, xtalform_site in xtalform_sites.xtalform_sites.items():
            xtalform_sites_meta[xtalform_site_id] = {
                Constants.META_XTALFORM_SITE_XTALFORM_ID: xtalform_site.xtalform_id,
                Constants.META_XTALFORM_SITE_CANONICAL_SITE_ID: xtalform_site.site_id,
                Constants.META_XTALFORM_SITE_LIGAND_CHAIN: xtalform_site.crystallographic_chain,
                Constants.META_XTALFORM_SITE_MEMBERS: {
                        Constants.META_DTAG: [lid.dtag for lid in xtalform_site.members],
                     Constants.META_CHAIN: [lid.chain for lid in xtalform_site.members],
                     Constants.META_RESIDUE: [lid.residue for lid in xtalform_site.members]
                     }
            }

        # Add the output aligned files
        for dtag, crystal in crystals.items():

            # Skip if no output for this dataset
            if dtag not in dataset_output_dict:
                continue

            crystal_output = new_meta[dtag] = {}

            # Otherwise iterate the output data structure, adding the aligned structure,
            # artefacts, xmaps and event maps to the metadata_file
            assigned_xtalform = assigned_xtalforms.get_xtalform_id(dtag)
            crystal_output[Constants.META_ASSIGNED_XTALFORM] = assigned_xtalform

            aligned_output = crystal_output[Constants.META_ALIGNED_FILES] = {}
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

        return new_meta


def main():
    parser = argparse.ArgumentParser(description='aligner')

    parser.add_argument('-d', '--version-dir', required=True, help="Path to version dir")
    parser.add_argument('-m', '--metadata_file', default=Constants.METADATA_XTAL_FILENAME, help="Metadata YAML file")
    parser.add_argument('-x', '--xtalforms', default=Constants.XTALFORMS_FILENAME, help="Crystal forms JSON file")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("aligner: ", args)

    logger = utils.Logger(logfile=args.log_file, level=args.log_level)

    a = Aligner(args.version_dir, args.metadata_file, args.xtalforms, logger=logger)
    num_errors, num_warnings = a.validate()

    if not args.validate:
        if num_errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            a.run()


if __name__ == "__main__":
    main()
