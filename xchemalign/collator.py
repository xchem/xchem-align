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

import argparse, os, shutil
import yaml

from . import utils, processor

from .utils import Constants

import pandas as pd
from pathlib import Path
import numpy as np
import gemmi


class Collator(processor.Processor):

    def __init__(self, config_file, logger=None):
        super(Collator, self).__init__(config_file, logger=logger)
        self.all_xtals = None
        self.new_or_updated_xtals = None

    def validate(self):

        self.logger.info('validating paths')
        num_errors, num_warnings = self.validate_paths()

        if num_errors == 0:
            self.logger.info('validating data')
            meta = self.validate_data()
        else:
            meta = {}

        return meta, len(self.errors), len(self.warnings)

    def run(self, meta):
        self.logger.info('Running collator...')
        v_dir = self.read_versions()
        if not v_dir:
            self.logger.error('Error with version dir. Please fix and try again.')
            return None
        self.logger.info('Using version dir {}'.format(v_dir))
        self.logger.info('Coping files ...')
        new_meta = self._copy_files(meta)
        self.logger.info('Munging the history ...')
        all_xtals, new_xtals = self._munge_history(meta)
        self.logger.info('Writing metadata ...')
        self._write_metadata(new_meta, all_xtals, new_xtals)
        self.logger.info('Copying config ...')
        self._copy_config()
        self.logger.info('Run complete')
        return new_meta

    def read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = Constants.VERSION_DIR_PREFIX + str(version)
            if (self.output_path / v_dir).is_dir():
                version += 1
            else:
                break
        if version == 1:
            self.logger.error('No version directory found. Please create one named upload_1')
            return None

        # the working version dir is one less than the current value
        version -= 1
        self.logger.info('Version is {}'.format(version))
        v_dir = Constants.VERSION_DIR_PREFIX + str(version)

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(1, version):
                self.logger.info('Reading metadata for version {}'.format(v))
                meta_file = self.output_path / (Constants.VERSION_DIR_PREFIX + str(v)) / Constants.METADATA_XTAL_FILENAME
                if meta_file.is_file():
                    with open(meta_file, 'r') as stream:
                        meta = yaml.safe_load(stream)
                        self.meta_history.append(meta)
                else:
                    self.logger.error('Metadata file {} not found'.format(meta_file))
                    return None

        self.version_dir = Path(v_dir)

        num_old_metas = len(self.meta_history)
        if num_old_metas:
            self.logger.info('Found {} metadata files from previous versions'.format(num_old_metas))

        return v_dir

    def _copy_files(self, meta):
        cryst_path = self.version_dir / Constants.META_XTAL_FILES
        ext_cryst_path = self.output_path / cryst_path
        self.logger.info('Using cryst_dir of', ext_cryst_path)
        if ext_cryst_path.exists():
            self.logger.info('removing old cryst_dir')
            shutil.rmtree(ext_cryst_path)
        self.logger.info('creating cryst_dir')
        os.makedirs(ext_cryst_path)

        num_event_maps = 0

        for xtal_name, data in meta[Constants.META_XTALS].items():
            dir = cryst_path / xtal_name
            os.makedirs(self.output_path / dir)

            xtal_files = data[Constants.META_XTAL_FILES]

            # handle the PDB file
            pdb = xtal_files[Constants.META_XTAL_PDB]
            if pdb:
                pdb_input = self.base_path / pdb
                pdb_name = xtal_name + '.pdb'
                pdb_output = dir / pdb_name
                f = shutil.copy2(pdb_input, self.output_path / pdb_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy PDB file {} to {}'.format(pdb_input, self.output_path / pdb_output))
                    return None
                digest = utils.gen_sha256(self.output_path / pdb_output)
                xtal_files[Constants.META_XTAL_PDB] = {Constants.META_FILE: str(pdb_output), Constants.META_SHA256: digest}
            else:
                self.logger.error('PDB entry missing for {}'.format(xtal_name))
                return meta

            # handle the MTZ file
            mtz = xtal_files.get(Constants.META_XTAL_MTZ)
            if mtz:
                mtz_input = self.base_path / mtz
                mtz_name = xtal_name + '.mtz'
                mtz_output = dir / mtz_name
                f = shutil.copy2(mtz_input, self.output_path / mtz_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy MTZ file {} to {}'.format(mtz_input, self.output_path / mtz_output))
                    return None
                digest = utils.gen_sha256(self.output_path / mtz_output)
                xtal_files[Constants.META_XTAL_MTZ] = {Constants.META_FILE: str(mtz_output), Constants.META_SHA256: digest}
            else:
                self.logger.warn('MTZ entry missing for {}'.format(xtal_name))

            # handle the CIF file
            cif = xtal_files.get(Constants.META_XTAL_CIF)
            if cif:
                cif_input = self.base_path / cif
                cif_name = xtal_name + '.cif'
                cif_output = dir / cif_name
                f = shutil.copy2(cif_input, self.output_path / cif_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy CIF file {} to {}'.format(cif_input, self.output_path / cif_output))
                    return None
                digest = utils.gen_sha256(self.output_path / cif_output)
                xtal_files[Constants.META_XTAL_CIF] = {Constants.META_FILE: str(cif_output), Constants.META_SHA256: digest}

                # # convert ligand PDB to SDF
                # # The ligand CIF file does not seem to be readable using OpenBabel so we resort to using the PDB
                # # that also seems to be generated but is not referenced in the database
                # sdf_file = os.path.join(dir, xtal_name + '.sdf')
                # ligand_pdb = cif_input[:-4] + '.pdb'
                # if os.path.isfile(ligand_pdb):
                #     count = obabel_utils.convert_molecules(ligand_pdb, 'pdb', sdf_file, 'sdf')
                #     if count:
                #         digest = utils.gen_sha256(sdf_file)
                #         xtal_files['ligand_sdf'] = {'file': sdf_file, 'sha256': digest}
                #     else:
                #         self.logger.warn('Ligand SDF file was not generated')
                # else:
                #     self.logger.warn('Ligand PDB file {} not found'.format(ligand_pdb))
            else:
                self.logger.warn('CIF entry missing for {}'.format(xtal_name))

            if pdb:
                event_tables = {}
                for input in self.inputs:
                    for panddas_file in input.panddas_event_file_paths:
                        pfp = panddas_file
                        # self.logger.info('Reading CSV:', pfp)
                        df = pd.read_csv(input.get_input_dir_path() / pfp)
                        event_tables[input.get_input_dir_path() / pfp] = df

                # find the best event maps, then copy them to the standard location and update the metadata
                if xtal_name in meta[Constants.META_XTALS]:
                    best_event_map_paths = self.get_dataset_event_maps(xtal_name, pdb_input, event_tables)
                    if best_event_map_paths:
                        p_paths = []
                        for k, tup in best_event_map_paths.items():
                            ccp4_file_path = tup[0]
                            # print('handling', xtal_name, ccp4_file_path)
                            ccp4_output = cryst_path / xtal_name / '{}_{}_{}.ccp4'.format(k[0], k[1], k[2])
                            self.logger.info('copying CCP4 file {} to {}'.format(ccp4_file_path, self.output_path / ccp4_output))
                            f = shutil.copy2(ccp4_file_path, self.output_path / ccp4_output, follow_symlinks=True)
                            if f:
                                digest = utils.gen_sha256(self.output_path / ccp4_output)
                                p_paths.append({
                                    Constants.META_FILE: str(ccp4_output),
                                    Constants.META_SHA256: digest,
                                    Constants.META_PROT_MODEL: int(k[0]),
                                    Constants.META_PROT_CHAIN: k[1],
                                    Constants.META_PROT_RES: k[2],
                                    Constants.META_PROT_INDEX: tup[1],
                                    Constants.META_PROT_BDC: tup[2]})
                                num_event_maps += 1
                            else:
                                self.logger.error('Failed to copy CCP4 file {} to {}'.format(ccp4_file_path, self.output_path / ccp4_output))

                        if p_paths:
                            meta[Constants.META_XTALS][xtal_name][Constants.META_XTAL_FILES][Constants.META_BINDING_EVENT] = p_paths
                else:
                    self.logger.warn('crystal {} not found in metadata - strange!'.format(xtal_name))

        self.logger.info('found {} event map files'.format(num_event_maps))

        return meta

    def _munge_history(self, meta):
        all_xtals = {}
        new_or_updated_xtals = {}

        # get any user defined overrides
        overrides = self.config.get('overrides', {})

        count = 0
        for metad in self.meta_history:
            count += 1
            self.logger.info('Munging metadata {}'.format(count))
            xtals = metad['crystals']
            total = 0
            for xtal_name, xtal_data in xtals.items():
                total += 1
                all_xtals[xtal_name] = xtal_data
            self.logger.info('Metadata {} has {} items'.format(count, total))

        count += 1
        self.logger.info('Munging current metadata')
        xtals = meta['crystals']
        total = 0
        for xtal_name, xtal_data in xtals.items():
            total += 1
            if xtal_name in all_xtals:
                old_xtal_data = all_xtals[xtal_name]
                status = self._get_xtal_status(xtal_name, old_xtal_data, xtal_data)
                xtal_data[Constants.META_STATUS] = status
                if status != Constants.META_STATUS_UNCHANGED:
                    new_or_updated_xtals[xtal_name] = xtal_data
            else:
                xtal_data[Constants.META_STATUS] = Constants.META_STATUS_NEW
                new_or_updated_xtals[xtal_name] = xtal_data
            all_xtals[xtal_name] = xtal_data

            # look for any deprecations
            xtals_overrides = overrides.get(Constants.META_XTALS, {})
            xtal_override = overrides.get(Constants.META_XTALS, {}).get(xtal_name)
            if xtal_override:
                status_override = xtal_override.get(Constants.META_STATUS)
                if status_override:
                    status = status_override.get(Constants.META_STATUS)
                    if status:
                        xtal_data[Constants.META_STATUS] = status
                        reason = status_override.get(Constants.META_REASON)
                        self.logger.info('status for xtal {} is overridden to be {}'.format(xtal_name, status))
                        if reason:
                            xtal_data[Constants.META_REASON] = reason
                        else:
                            self.logger.warn('status is overridden, but no reason was given')
                    else:
                        self.logger.warn('status is declared to be overridden, but no new status was given')

        self.logger.info('Metadata {} has {} items'.format(count, total))
        self.logger.info('Munging resulted in {} total xtals, {} are new or updated'.format(
            len(all_xtals), len(new_or_updated_xtals)))

        self.all_xtals = all_xtals
        self.new_or_updated_xtals = new_or_updated_xtals
        return all_xtals, new_or_updated_xtals

    def _get_xtal_status(self, xtal_name, old_data, new_data):
        """
        Compare status using the SHA256 digest of the PDB file.
        This should be reliable, but will be subject to saying the entry is updated even if there is a trivial
        change.

        :param xtal_name:
        :param old_data:
        :param new_data:
        :return:
        """
        old_pdb_data = old_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        new_pdb_data = new_data[Constants.META_XTAL_FILES][Constants.META_XTAL_PDB]
        old_sha256 = old_pdb_data[Constants.META_SHA256]
        new_sha256 = new_pdb_data[Constants.META_SHA256]
        if old_sha256 == new_sha256:
            return Constants.META_STATUS_UNCHANGED
        else:
            return Constants.META_STATUS_SUPERSEDES

    def _get_xtal_status_old(self, xtal_name, old_data, new_data):
        """
        Compare status using the last_updated property.
        Even if this can be relied on, it is not present for manual entries.

        :param xtal_name:
        :param old_data:
        :param new_data:
        :return:
        """
        old_date = old_data.get(Constants.META_LAST_UPDATED)
        if old_date:
            new_date = new_data.get(Constants.META_LAST_UPDATED)
        if not old_date or not new_date:
            self.logger.warn('Dates not defined for {}, must assume xtal is updated'.format(xtal_name))
            return Constants.META_STATUS_SUPERSEDES
        elif utils.to_datetime(new_date) > utils.to_datetime(old_date):
            self.logger.info('Xtal {} is updated'.format(xtal_name))
            return Constants.META_STATUS_SUPERSEDES
        else:
            # self.logger.info('Xtal {} is unchanged'.format(xtal_name))
            return Constants.META_STATUS_UNCHANGED

    def _write_metadata(self, meta, all_xtals, new_xtals):
        f = self.output_path / self.version_dir / Constants.METADATA_XTAL_FILENAME
        with open(f, 'w') as stream:
            yaml.dump(meta, stream, sort_keys=False)
        f = self.output_path / self.version_dir / 'all_xtals.yaml'
        with open(f, 'w') as stream:
            yaml.dump(all_xtals, stream, sort_keys=False)
            f = self.output_path / self.version_dir / 'new_xtals.yaml'
        with open(f, 'w') as stream:
            yaml.dump(new_xtals, stream, sort_keys=False)

    def _copy_config(self):
        f = shutil.copy2(self.config_file, self.output_path / self.version_dir)
        if not f:
            print('Failed to copy config file to {}'.format((self.output_path / self.version_dir)))
            return False
        return True

    def get_ligand_coords(self, structure: gemmi.Structure, ) -> dict[tuple[str, str, str], np.array]:
        ligand_coords = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name in Constants.LIGAND_NAMES:

                        poss = []
                        for atom in residue:
                            pos = atom.pos
                            poss.append([pos.x, pos.y, pos.z])

                        arr = np.array(poss)
                        mean = np.mean(arr, axis=0)
                        ligand_coords[(model.name, chain.name, residue.seqid.num)] = mean

        return ligand_coords

    def get_closest_event_map(
            self,
            xtal_name: str,
            ligand_coord: np.array,
            event_tables: dict[Path, pd.DataFrame],
    ) -> Path:
        distances = {}
        data = {}
        for pandda_path, event_table in event_tables.items():
            # print('Processing', xtal_name, pandda_path)
            dataset_events = event_table[event_table[Constants.EVENT_TABLE_DTAG] == xtal_name]
            for idx, row in dataset_events.iterrows():
                event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
                bdc = row[Constants.EVENT_TABLE_BDC]
                x,y,z = row[Constants.EVENT_TABLE_X], row[Constants.EVENT_TABLE_Y], row[Constants.EVENT_TABLE_Z]
                distance = np.linalg.norm(np.array([x,y,z]).flatten() - ligand_coord.flatten())
                # print('Distance:', distance)
                for template in Constants.EVENT_MAP_TEMPLATES:
                    event_map_path = pandda_path.parent.parent / Constants.PROCESSED_DATASETS_DIR / xtal_name / template.format(
                        dtag=xtal_name,
                        event_idx=event_idx,
                        bdc=bdc
                    )
                    if event_map_path.exists() and event_map_path.is_file():
                        distances[event_map_path] = distance
                        data[event_map_path] = (event_idx, bdc)
                        continue

        if distances:
            k = min(distances, key=lambda _key: distances[_key])
            return k, data[k]
        else:
            return None, None

    def get_dataset_event_maps(self, xtal_name: str, pdb_file: Path, event_tables: dict[Path, pd.DataFrame]
                               ) -> dict[tuple[str, str, str], Path]:
        # Get the relevant structure
        # self.logger.info('Reading', xtal_name, pdb_file)
        # self.logger.info('Using {} event tables'.format(len(event_tables)))

        structure = gemmi.read_structure(str(pdb_file))

        # Get the coordinates of ligands
        ligand_coords = self.get_ligand_coords(structure)

        # Get the closest events within some reasonable radius
        closest_event_maps = {}
        for ligand_key, ligand_coord in ligand_coords.items():
            # print('coord:', ligand_coord)
            closest_event_map, data = self.get_closest_event_map(xtal_name, ligand_coord, event_tables)
            if closest_event_map:
                # print('closest:', closest_event_map)
                closest_event_maps[ligand_key] = (closest_event_map, data[0], data[1])

        return closest_event_maps


def main():

    parser = argparse.ArgumentParser(description='collator')

    parser.add_argument('-c', '--config-file', default='config.yaml', help="Configuration file")
    parser.add_argument('-l', '--log-file', help="File to write logs to")
    parser.add_argument('--log-level', type=int, default=0, help="Logging level")
    parser.add_argument('-v', '--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("collator: ", args)

    c = Collator(args.config_file, logger=logger)

    meta, num_errors, num_warnings = c.validate()

    if not args.validate:
        if num_errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            c.run(meta)

if __name__ == "__main__":
    main()
