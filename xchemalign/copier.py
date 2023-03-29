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

import os, argparse, shutil
from . import dbreader, validator


def copy_files(base_dir, input_dir, output_dir):
    if base_dir:
        dbfile = os.path.join(base_dir, input_dir, 'processing', 'database', 'soakDBDataFile.sqlite')
    else:
        dbfile = os.path.join(input_dir, 'processing', 'database', 'soakDBDataFile.sqlite')
    df = dbreader.filter_dbmeta(dbfile)
    count = 0
    copied = 0
    for index, row in df.iterrows():
        count += 1
        xtal_name = row['CrystalName']
        xtal_dir = validator.generate_xtal_dir(input_dir, xtal_name)
        print('processing {} {}'.format(count, xtal_name))
        for colname in ['RefinementPDB_latest', 'RefinementMTZ_latest', 'RefinementCIF']:
            file = row[colname]
            if file:
                ok = _copy_file(file, base_dir, input_dir, xtal_dir, output_dir)
                if ok:
                    copied += 1
        # for the ligand CIF file also copy the corresponding PDB file
        if file:
            ok = _copy_file(file[:-3] + 'pdb', base_dir, input_dir, xtal_dir, output_dir)
            if ok:
                copied += 1
    print('Copied {} files'.format(copied))


def _copy_file(filepath, base_dir, input_dir, xtal_dir, output_dir):

    # print('handling', filepath)
    inputpath, outputpath = validator.generate_filenames(filepath, xtal_dir, output_dir)
    if base_dir:
        full_inputpath = base_dir + '/' + inputpath
    else:
        full_inputpath = inputpath
    # print('copying', full_inputpath, outputpath)

    if not os.path.isfile(full_inputpath):
        print('File {} not found'.format(full_inputpath))
        return False

    os.makedirs(os.path.dirname(outputpath), exist_ok=True)
    f = shutil.copy2(full_inputpath, outputpath, follow_symlinks=True)
    if not f:
        print('Failed to copy file {} to {}'.format(inputpath, outputpath))
        return False
    return True


def find_panddas(dirname):
    """
    Find the PanDDAs event maps (.ccp4 files) that are relevant to the data contained in this directory
    :param dirname: The directory in which to look e.g. /dls/labxchem/data/2020/lb18145-153
    :return: A list of paths to the .ccp4 files?
    """
    # TODO Conor to implement.

    # Will be used by copier.py to copy those .ccp4 files from Diamond.
    # The collator.py tool will copy those files to a standard location under outputs/<target_name>/upload_n/crystallographic
    # and list those files in the metadata.yaml file.

    # Question: will the returned value be a list of .ccp4 files for the whole set of data, or will it be crystal specific?

    pass


def main():

    parser = argparse.ArgumentParser(description='copier')

    parser.add_argument('-b', '--base-dir', help="Base directory")
    parser.add_argument('-i', '--input-dir', required=True, help="Input directory (below base-dir)")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")

    args = parser.parse_args()
    print("copier: ", args)

    copy_files(args.base_dir, args.input_dir, args.output_dir)

    # _copy_file('compound/foo.cif', 'dls/labxchem/data/lb18145/lb18145-216', 'data', 'data/outputs')
    # _copy_file('/dls/labxchem/data/lb18145/lb18145-216/compound/foo.cif', 'dls/labxchem/data/lb18145/lb18145-216', 'data', 'data/outputs')


if __name__ == "__main__":
    main()
