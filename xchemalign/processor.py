import argparse, os, shutil
import yaml

from . import utils, obabel_utils
from .validator import Validator

_METADATA_FILENAME = 'metadata.yaml'
_CONFIG_FILENAME = 'config.yaml'
_VERSION_DIR_PREFIX = 'upload_'

class Processor:

    def __init__(self, config_file, logger=None):
        self.config_file = config_file

        with open(config_file, 'r') as stream:
            config = yaml.safe_load(stream)

        self.input_dirs = config['input_dirs']
        self.output_dir = config['output_dir']
        self.target_name = config['target_name']
        self.version_dir = None
        self.meta_history = []
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def validate(self):
        v = Validator(self.input_dirs, self.output_dir, self.target_name, logger=self.logger)
        meta = v.validate_all()
        infos, warnings, errors = self.logger.get_num_messages()
        return meta, warnings, errors

    def run(self, meta):
        self.logger.log('Running ...', level=0)
        v_dir = self._read_versions()
        self.logger.log('Using version dir {}'.format(v_dir), level=0)
        self.logger.log('Coping files ...', level=0)
        new_meta = self._copy_files(meta)
        self.logger.log('Writing metadata ...', level=0)
        self._write_metadata(new_meta)
        return new_meta

    def _read_versions(self):
        # find out which version dirs exist
        version = 1
        while True:
            v_dir = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(version))
            if os.path.isdir(v_dir):
                version += 1
            else:
                break

        # the working version dir is one less than the current value
        version -= 1
        v_dir = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(version))

        # read the metadata from the earlier versions
        if version > 1:
            for v in range(version -1)
                meta_file = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(v), _METADATA_FILENAME)
                if os.path.isfile(meta_file):
                    with open(meta_file, 'r') as stream:
                        meta = yaml.safe_load(stream)
                        self.meta_history.append(meta)

        return v_dir


# def _create_version_dir(self):
    #     version = 1
    #     while True:
    #         v_dir = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(version))
    #         if os.path.exists(v_dir):
    #             meta_file = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(version), _METADATA_FILENAME)
    #             with open(meta_file, 'r') as stream:
    #                 try:
    #                     meta = yaml.safe_load(stream)
    #                 except yaml.YAMLError as exc:
    #                     self.logger.log('Failed to read metadata file', meta_file)
    #                     print(exc)
    #                 self.meta_history.append(meta)
    #             version += 1
    #             dir_exists = True
    #         else:
    #             dir_exists = False
    #             break
    #
    #     if version == 1 or self.new_version:
    #         # we don't have any versions so create the first
    #         pass
    #     else:
    #         # no new version so we must remove the latest one
    #         version -= 1
    #         v_dir = os.path.join(self.output_dir, _VERSION_DIR_PREFIX + str(version))
    #         self.logger.log('Removing current version dir {}'.format(v_dir), level=0)
    #         shutil.rmtree(v_dir)
    #         # also remove the last metadata item from the history
    #         self.meta_history.pop()
    #
    #     self.logger.log('Creating current version dir {}'.format(v_dir), level=0)
    #     self.logger.log('Metadata history has {} items'.format(len(self.meta_history)), level=0)
    #
    #     os.makedirs(v_dir)
    #     self.version_dir = v_dir
    #     return v_dir

    def _copy_files(self, meta):
        cryst_dir = os.path.join(self.version_dir, 'crystallographic')
        self.logger.log('Using cryst_dir of', cryst_dir, level=0)
        if os.path.exists(cryst_dir):
            self.logger.log('removing old cryst_dir', level=0)
            shutil.rmtree(cryst_dir)
        self.logger.log('creating cryst_dir', level=0)
        os.makedirs(cryst_dir)

        for name, data in meta['crystals'].items():
            dir = os.path.join(cryst_dir, name)
            os.makedirs(dir)

            xtal_files = data['crystallographic_files']

            # handle the PDB file
            pdb = xtal_files['xtal_pdb']
            pdb_input = self.input_dir + pdb
            pdb_output = os.path.join(dir, name + '.pdb')
            f = shutil.copy2(pdb_input, pdb_output, follow_symlinks=True)
            if not f:
                self.logger.log('Failed to copy PDB file {} to {}'.format(pdb_input, pdb_output), level=2)
                return None
            digest = utils.gen_sha256(pdb_output)
            xtal_files['xtal_pdb'] = {'file': pdb_output, 'sha256': digest}

            # handle the MTZ file
            mtz = xtal_files['xtal_mtz']
            mtz_input = self.input_dir + mtz
            mtz_output = os.path.join(dir, name + '.mtz')
            f = shutil.copy2(mtz_input, mtz_output, follow_symlinks=True)
            if not f:
                self.logger.log('Failed to copy MTZ file {} to {}'.format(mtz_input, mtz_output), level=2)
                return None
            digest = utils.gen_sha256(mtz_output)
            xtal_files['xtal_mtz'] = {'file': mtz_output, 'sha256': digest}

            # handle the CIF file
            cif = xtal_files['ligand_cif']
            cif_input = self.input_dir + cif
            cif_output = os.path.join(dir, name + '.cif')
            f = shutil.copy2(cif_input, cif_output, follow_symlinks=True)
            if not f:
                self.logger.log('Failed to copy CIF file {} to {}'.format(cif_input, cif_output), level=2)
                return None
            digest = utils.gen_sha256(cif_output)
            xtal_files['ligand_cif'] = {'file': cif_output, 'sha256': digest}

            # convert ligand PDB to SDF
            # The ligand CIF file does not seem to be readable using OpenBabel so we resort to using the PDB
            # that also seems to be generated but is not referenced in the database
            sdf_file = os.path.join(dir, name + '.sdf')
            ligand_pdb = cif_input[:-4] + '.pdb'
            if os.path.isfile(ligand_pdb):
                count = obabel_utils.convert_molecules(ligand_pdb, 'pdb', sdf_file, 'sdf')
                if count:
                    digest = utils.gen_sha256(sdf_file)
                    xtal_files['ligand_sdf'] = {'file': sdf_file, 'sha256': digest}
                else:
                    self.logger.log('Ligand SDF file was not generated', level=1)
            else:
                self.logger.log('Ligand PDB file {} not found'.format(ligand_pdb), level=1)

        return meta

    def _perform_alignments(self, meta):
        # do Conor's stuff
        return meta

    def _extract_components(self, meta):
        # extract molfile, apo pdbs etc.
        return meta

    def _write_metadata(self, meta):
        f = os.path.join(self.version_dir, _METADATA_FILENAME)
        with open(f, 'w') as stream:
            yaml.dump(meta, stream, sort_keys=False)


def main():

    parser = argparse.ArgumentParser(description='processor')

    parser.add_argument('-c', '--config-file', required=True, help="Configuration file")
    parser.add_argument('-l', '--log-file', help="Sqlite DB file")
    parser.add_argument('--log-level', type=int, default=1, help="Logging level")
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("processor: ", args)

    if args.log_file:
        logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    else:
        logger = None

    p = Processor(args.config_file, logger=logger)

    meta, warnings, errors = p.validate()

    print(meta)

    if not args.validate:
        if errors:
            print('There are errors, cannot continue')
            exit(1)
        else:
            p.run(meta)


if __name__ == "__main__":
    main()
