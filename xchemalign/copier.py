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
        for colname in ['RefinementCIF', 'RefinementPDB_latest', 'RefinementMTZ_latest']:
            file = row[colname]
            if file:
                ok = _copy_file(file, base_dir, input_dir, xtal_dir, output_dir)
                if ok:
                    copied += 1
    print('Copied {} files'.format(copied))


def _copy_file(filepath, base_dir, input_dir, xtal_dir, output_dir):

    # print('handling', filepath)
    inputpath, outputpath = validator.generate_filenames(filepath, base_dir, xtal_dir, output_dir)
    # print('copying', inputpath, outputpath)

    if not os.path.isfile(inputpath):
        print('File {} not found'.format(inputpath))
        return False


    os.makedirs(os.path.dirname(outputpath), exist_ok=True)
    f = shutil.copy2(inputpath, outputpath, follow_symlinks=True)
    if not f:
        print('Failed to copy file {} to {}'.format(inputpath, outputpath))
        return False
    #print('Copied file {} to {}'.format(path, outputpath))
    return True


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
