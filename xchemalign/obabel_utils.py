from openbabel import pybel
ob = pybel.ob


def convert_molecules(in_file, in_format, out_file, out_format):
    print('converting {} to {}'.format(in_file, out_file))
    mols = pybel.readfile(in_format, in_file)
    count = 0
    with pybel.Outputfile(out_format, out_file, overwrite=True) as writer:
        for mol in mols:
            count += 1
            writer.write(mol)
        print('Converted {} mols'.format(count))
    return count