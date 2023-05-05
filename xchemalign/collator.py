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


class Collator(processor.Processor):

    def __init__(self, config_file, logger=None):
        super(Collator, self).__init__(config_file, logger=logger)
        self.all_xtals = None
        self.new_or_updated_xtals = None

    def validate(self):

        num_errors, num_warnings = self.validate_paths()

        if num_errors == 0:
            meta = self.validate_soakdb_data()
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
            v_dir = self.output_path / (processor.VERSION_DIR_PREFIX + str(version))
            if v_dir.is_dir():
                version += 1
            else:
                break
        if version == 1:
            self.logger.error('No version directory found. Please create one named upload_1')
            return None

        # the working version dir is one less than the current value
        version -= 1
        self.logger.info('Version is {}'.format(version))
        v_dir = self.output_path / (processor.VERSION_DIR_PREFIX + str(version))

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(1, version):
                self.logger.info('Reading metadata for version {}'.format(v))
                meta_file = self.output_path / (processor.VERSION_DIR_PREFIX + str(v)) / processor.METADATA_FILENAME
                if meta_file.is_file():
                    with open(meta_file, 'r') as stream:
                        meta = yaml.safe_load(stream)
                        self.meta_history.append(meta)
                else:
                    self.logger.error('Metadata file {} not found'.format(meta_file))
                    return None

        self.version_dir = v_dir

        num_old_metas = len(self.meta_history)
        if num_old_metas:
            self.logger.info('Found {} metadata files from previous versions'.format(num_old_metas))

        return v_dir

    def _copy_files(self, meta):
        cryst_dir = self.version_dir / 'crystallographic'
        self.logger.info('Using cryst_dir of', cryst_dir)
        if cryst_dir.exists():
            self.logger.info('removing old cryst_dir')
            shutil.rmtree(cryst_dir)
        self.logger.info('creating cryst_dir')
        os.makedirs(cryst_dir)

        for name, data in meta['crystals'].items():
            dir = os.path.join(cryst_dir, name)
            os.makedirs(dir)

            xtal_files = data['crystallographic_files']

            # handle the PDB file
            pdb = xtal_files['xtal_pdb']
            if pdb:
                pdb_input = self.base_path / pdb
                pdb_output = os.path.join(dir, name + '.pdb')
                f = shutil.copy2(pdb_input, pdb_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy PDB file {} to {}'.format(pdb_input, pdb_output))
                    return None
                digest = utils.gen_sha256(pdb_output)
                xtal_files['xtal_pdb'] = {'file': pdb_output, 'sha256': digest}
            else:
                self.logger.warn('PDB entry missing for {}'.format(name))

            # handle the MTZ file
            mtz = xtal_files['xtal_mtz']
            if mtz:
                mtz_input = self.base_path / mtz
                mtz_output = os.path.join(dir, name + '.mtz')
                f = shutil.copy2(mtz_input, mtz_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy MTZ file {} to {}'.format(mtz_input, mtz_output))
                    return None
                digest = utils.gen_sha256(mtz_output)
                xtal_files['xtal_mtz'] = {'file': mtz_output, 'sha256': digest}
            else:
                self.logger.warn('MTZ entry missing for {}'.format(name))

            # handle the CIF file
            cif = xtal_files['ligand_cif']
            if cif:
                cif_input = self.base_path / cif
                cif_output = os.path.join(dir, name + '.cif')
                f = shutil.copy2(cif_input, cif_output, follow_symlinks=True)
                if not f:
                    self.logger.error('Failed to copy CIF file {} to {}'.format(cif_input, cif_output))
                    return None
                digest = utils.gen_sha256(cif_output)
                xtal_files['ligand_cif'] = {'file': cif_output, 'sha256': digest}

                # # convert ligand PDB to SDF
                # # The ligand CIF file does not seem to be readable using OpenBabel so we resort to using the PDB
                # # that also seems to be generated but is not referenced in the database
                # sdf_file = os.path.join(dir, name + '.sdf')
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
                self.logger.warn('CIF entry missing for {}'.format(name))

        return meta

    def _munge_history(self, meta):
        all_xtals = {}
        new_or_updated_xtals = {}

        # handle any user defined deprecations
        if 'overrides' in self.config and 'deprecations' in self.config['overrides']:
            deprecations = self.config['overrides']['deprecations']
        else:
            deprecations = {}
        self.logger.info('{} deprecations were defined'.format(len(deprecations)))

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
                old_date = old_xtal_data['last_updated']
                new_date = xtal_data['last_updated']
                if not old_date or not new_date:
                    self.logger.warn('Dates not defined for {}, must assume xtal is updated {} {}'.format(xtal_name, ))
                    xtal_data['status'] = 'supersedes'
                    new_or_updated_xtals[xtal_name] = xtal_data
                elif utils.to_datetime(new_date) > utils.to_datetime(old_date):
                    self.logger.info('Xtal {} is updated'.format(xtal_name))
                    xtal_data['status'] = 'supersedes'
                    new_or_updated_xtals[xtal_name] = xtal_data
                else:
                    # self.logger.info('Xtal {} is unchanged'.format(xtal_name))
                    xtal_data['status'] = 'unchanged'
            else:
                xtal_data['status'] = 'new'
                new_or_updated_xtals[xtal_name] = xtal_data
            all_xtals[xtal_name] = xtal_data

            # look for any deprecations
            if xtal_name in deprecations:
                xtal_data['status'] = 'deprecated'
                xtal_data['reason'] = deprecations[xtal_name]
                self.logger.info('Deprecating xtal {}'.format(xtal_name))

        self.logger.info('Metadata {} has {} items'.format(count, total))
        self.logger.info('Munging resulted in {} total xtals, {} are new or updated'.format(
            len(all_xtals), len(new_or_updated_xtals)))

        self.all_xtals = all_xtals
        self.new_or_updated_xtals = new_or_updated_xtals
        return all_xtals, new_or_updated_xtals

    def _write_metadata(self, meta, all_xtals, new_xtals):
        f = os.path.join(self.version_dir, processor.METADATA_FILENAME)
        with open(f, 'w') as stream:
            yaml.dump(meta, stream, sort_keys=False)
        f = os.path.join(self.version_dir, 'all_xtals.yaml')
        with open(f, 'w') as stream:
            yaml.dump(all_xtals, stream, sort_keys=False)
            f = os.path.join(self.version_dir, 'new_xtals.yaml')
        with open(f, 'w') as stream:
            yaml.dump(new_xtals, stream, sort_keys=False)

    def _copy_config(self):
        f = shutil.copy2(self.config_file, self.version_dir)
        if not f:
            print('Failed to copy config file to {}'.format(self.version_dir))
            return False
        return True


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
