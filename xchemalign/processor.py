import argparse, os, shutil
import yaml

from . import utils, obabel_utils
from .validator import Validator

_metadata_filename = 'metadata.yaml'


class Processor:

    def __init__(self, input_dir, output_dir, target_name, dbfile, new_version=True, logger=None):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.target_name = target_name
        self.dbfile = dbfile
        self.new_version = new_version
        self.version_dir = None
        self.meta_history = []
        if logger:
            self.logger = logger
        else:
            self.logger = utils.Logger()

    def validate(self):
        v = Validator(self.input_dir, self.output_dir, self.target_name, self.dbfile, logger=self.logger)
        meta = v.validate_all()
        infos, warnings, errors = self.logger.get_num_messages()
        return meta, warnings, errors

    def run(self, meta):
        self.logger.log('Running ...', level=0)
        v_dir = self._create_version_dir()
        self.logger.log('Using version dir {}'.format(v_dir), level=0)
        self.logger.log('Coping files ...', level=0)
        new_meta = self._copy_files(meta)
        self.logger.log('Performing alignments ...', level=0)
        new_meta = self._perform_alignments(new_meta)
        self.logger.log('Extracting components ...', level=0)
        new_meta = self._extract_components(new_meta)
        self.logger.log('Writing metadata ...', level=0)
        self._write_metadata(new_meta)
        return new_meta

    def _create_version_dir(self):
        version = 1
        while True:
            v_dir = os.path.join(self.output_dir, 'v' + str(version))
            if os.path.exists(v_dir):
                meta_file = os.path.join(self.output_dir, 'v' + str(version), _metadata_filename)
                with open(meta_file, 'r') as stream:
                    try:
                        meta = yaml.safe_load(stream)
                    except yaml.YAMLError as exc:
                        self.logger.log('Failed to read metadata file', meta_file)
                        print(exc)
                    self.meta_history.append(meta)
                version += 1
                dir_exists = True
            else:
                dir_exists = False
                break

        if version == 1 or self.new_version:
            # we don't have any versions so create the first
            pass
        else:
            # no new version so we must remove the latest one
            version -= 1
            v_dir = os.path.join(self.output_dir, 'v' + str(version))
            self.logger.log('Removing current version dir {}'.format(v_dir), level=0)
            shutil.rmtree(v_dir)
            # also remove the last metadata item from the history
            self.meta_history.pop()

        self.logger.log('Creating current version dir {}'.format(v_dir), level=0)
        self.logger.log('Metadata history has {} items'.format(len(self.meta_history)), level=0)

        os.makedirs(v_dir)
        self.version_dir = v_dir
        return v_dir

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
        f = os.path.join(self.version_dir, _metadata_filename)
        with open(f, 'w') as stream:
            yaml.dump(meta, stream, sort_keys=False)


def main():

    parser = argparse.ArgumentParser(description='processor')

    parser.add_argument('-i', '--input-dir', required=True, help="Input directory")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")
    parser.add_argument('-t', '--target-name', required=True, help="Target name")
    parser.add_argument('-d', '--db-file', required=True, help="Sqlite DB file")
    parser.add_argument('-l', '--log-file', help="Sqlite DB file")
    parser.add_argument('--log-level', type=int, default=1, help="Logging level")
    parser.add_argument('--new-version', action='store_true', help='Create a new version')
    parser.add_argument('--validate', action='store_true', help='Only perform validation')

    args = parser.parse_args()
    print("processor: ", args)

    if args.log_file:
        logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    else:
        logger = None

    p = Processor(args.input_dir, args.output_dir, args.target_name, args.db_file,
                  new_version=args.new_version, logger=logger)

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
