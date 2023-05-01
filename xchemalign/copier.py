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
from pathlib import Path

import pandas as pd

from . import dbreader, validator


def _generate_path(base_dir, input_dir, file):
    if input_dir:
        if base_dir:
            return base_dir + '/' + input_dir + '/' + file
        else:
            return input_dir + '/' + file
    else:
        if base_dir:
            return base_dir + file
        else:
            return file


def copy_files(base_dir, input_dir, output_dir, soakdb_file, panddas_files):

    dbfile = _generate_path(base_dir, input_dir, soakdb_file)

    print('Reading soakdb file', dbfile)
    df = dbreader.filter_dbmeta(dbfile)
    count = 0
    copied = 0
    datasets = {}
    for index, row in df.iterrows():
        count += 1
        xtal_name = row['CrystalName']
        xtal_dir = validator.generate_xtal_dir(input_dir, xtal_name)
        print('processing {} {}'.format(count, xtal_name))

        file = row['RefinementPDB_latest']
        if file:
            ok = _copy_file(file, base_dir, input_dir, xtal_dir, output_dir)
            if ok:
                copied += 1
                datasets[xtal_name] = file
                # if PDB is OK then continue with the other files
                file = row['RefinementMTZ_latest']
                if file:
                    ok = _copy_file(file, base_dir, input_dir, xtal_dir, output_dir)
                    if ok:
                        copied += 1
                file = row['RefinementCIF']
                if file:
                    ok = _copy_file(file, base_dir, input_dir, xtal_dir, output_dir)
                    if ok:
                        copied += 1
                        # for the ligand CIF file also copy the corresponding PDB file
                        ok = _copy_file(file[:-3] + 'pdb', base_dir, input_dir, xtal_dir, output_dir)
                        if ok:
                            copied += 1

    # copy the specified csv files with the panddas info
    print('Copying panddas csv files')
    for panddas_file in panddas_files:
        if base_dir:
            f = base_dir + input_dir + '/' + panddas_file
        else:
            f = input_dir + '/' + panddas_file
        ok = _copy_csv(panddas_file, base_dir, input_dir, output_dir)

    # copy the relevant panddas event map files
    _copy_panddas(datasets, base_dir, input_dir, output_dir, panddas_files)

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


def _copy_csv(filepath, base_dir, input_dir, output_dir):
    if base_dir:
        csv_src = base_dir + '/' + input_dir + '/' + filepath
    else:
        csv_src = input_dir + '/' + filepath
    if not os.path.isfile(csv_src):
        print('File {} not found'.format(csv_src))
        return False
    dest_file = output_dir + '/' + input_dir + '/' + filepath
    dest_dir = os.path.dirname(dest_file)
    print('Copying', csv_src, 'to', dest_dir)
    os.makedirs(dest_dir, exist_ok=True)
    f = shutil.copy2(csv_src, dest_dir + '/', follow_symlinks=True)
    if not f:
        print('Failed to copy file {} to {}'.format(dest_dir, dest_dir))
        return False
    return True


def _copy_panddas(datasets, base_dir, input_dir, output_dir, panddas_files):
    """
    Find the PanDDAs event maps (.ccp4 files) that are relevant to the data contained in this directory
    :param dirname: The directory in which to look e.g. /dls/labxchem/data/2020/lb18145-153
    :return: A list of paths to the .ccp4 files?
    """

    if not panddas_files:
        return

    print(len(datasets), 'datasets')

    panddas_dict = {}
    for panddas_file in panddas_files:
        f = _generate_path(base_dir, input_dir, panddas_file)
        print('Reading CSV:', f)
        df = pd.read_csv(f)
        print('Data frame shape:', df.shape)
        panddas_dict[Path(f)] = df

        for xtal_name, pdb_path in datasets.items():

            dataset_events = df[df[Constants.EVENT_TABLE_DTAG] == xtal_name]
            for idx, row in dataset_events.iterrows():
                print('  event', idx, xtal_name)
                event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
                bdc = row[Constants.EVENT_TABLE_BDC]
                event_map_path = Path(panddas_file).parent / Constants.PROCESSED_DATASETS_DIR / xtal_name / Constants.EVENT_MAP_TEMPLATE.format(
                    dtag=xtal_name,
                    event_idx=event_idx,
                    bdc=bdc
                )
                print(event_map_path)


from typing import Protocol
from pathlib import Path
import numpy as np
import gemmi


class DatasetInterface(Protocol):
    dtag: str
    pdb: Path


class Constants:
    EVENT_TABLE_DTAG = "dtag"
    EVENT_TABLE_EVENT_IDX = "event_idx"
    EVENT_TABLE_X = "x"
    EVENT_TABLE_Y = "y"
    EVENT_TABLE_Z = "z"
    EVENT_TABLE_BDC = "1-BDC"
    LIGAND_NAMES = ["LIG", "XXX"]
    PROCESSED_DATASETS_DIR = "processed_datasets"
    EVENT_MAP_TEMPLATE = "{dtag}-event_{event_idx}_1-BDC_{bdc}_map.native.ccp4"


def get_ligand_coords(structure: gemmi.Structure, ) -> dict[tuple[str, str, str], np.array]:
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
                    ligand_coords[(model.name, chain.name, residue.name)] = mean

    return ligand_coords


def get_closest_event_map(
        xtal_name: str,
        ligand_coord: np.array,
        event_tables: dict[Path, pd.DataFrame],
) -> Path:
    distances = {}
    for pandda_path, event_table in event_tables.items():
        print('Processing', pandda_path)
        dataset_events = event_table[event_table[Constants.EVENT_TABLE_DTAG] == xtal_name]
        for idx, row in dataset_events.iterrows():
            event_idx = row[Constants.EVENT_TABLE_EVENT_IDX]
            bdc = row[Constants.EVENT_TABLE_BDC]
            x,y,z = row[Constants.EVENT_TABLE_X], row[Constants.EVENT_TABLE_Y], row[Constants.EVENT_TABLE_Z]
            distance = np.linalg.norm(np.array([x,y,z]).flatten() - ligand_coord.flatten())
            print('Distance:', distance)
            event_map_path = pandda_path / Constants.PROCESSED_DATASETS_DIR / xtal_name / Constants.EVENT_MAP_TEMPLATE.format(
                dtag=xtal_name,
                event_idx=event_idx,
                bdc=bdc
            )
            distances[event_map_path] = distance

    return min(distances, key=lambda _key: distances[_key])


def get_dataset_event_maps(
        xtal_name: str,
        pdb_file: Path,
        base_dir: str,
        input_dir: str,
        event_tables: dict[Path, pd.DataFrame],
) -> dict[tuple[str, str, str], Path]:
    # Get the relevant structure
    p = _generate_path(base_dir, None, pdb_file)
    print('Reading', xtal_name, p)
    structure = gemmi.read_structure(p)

    # Get the coordinates of ligands
    ligand_coords = get_ligand_coords(structure)

    # Get the closest events within some reasonable radius
    closest_event_maps = {}
    for ligand_key, ligand_coord in ligand_coords.items():
        print('coord:', ligand_coord)
        closest_event_map = get_closest_event_map(xtal_name, ligand_coord, event_tables)
        closest_event_maps[ligand_key] = closest_event_map

    return closest_event_maps


def test_function_get_all_event_maps(datasets: dict[str, DatasetInterface], pandda_event_tables: dict[Path, pd.DataFrame]) -> dict[str, dict[tuple[str, str, str], Path]]:

    datasets_event_maps = {}
    for dtag, dataset in datasets.items():
        dataset_event_maps = get_dataset_event_maps(
            dataset,
            pandda_event_tables
        )
        datasets_event_maps[dtag] = dataset_event_maps

    return datasets_event_maps





def main():
    parser = argparse.ArgumentParser(description='copier')

    parser.add_argument('-b', '--base-dir', help="Base directory")
    parser.add_argument('-i', '--input-dir', required=True, help="Input directory (below base-dir)")
    parser.add_argument('-s', '--soakdb-file', default='processing/database/soakDBDataFile.sqlite',
                        help="Relative path to soakdb file")
    parser.add_argument('-p', '--panddas-files', nargs='*', help="Relative path to CSV files with panddas data")
    parser.add_argument('-o', '--output-dir', required=True, help="Output directory")

    args = parser.parse_args()
    print("copier: ", args)

    copy_files(args.base_dir, args.input_dir, args.output_dir, args.soakdb_file, args.panddas_files)

    # _copy_file('compound/foo.cif', 'dls/labxchem/data/lb18145/lb18145-216', 'data', 'data/outputs')
    # _copy_file('/dls/labxchem/data/lb18145/lb18145-216/compound/foo.cif', 'dls/labxchem/data/lb18145/lb18145-216', 'data', 'data/outputs')


if __name__ == "__main__":
    main()
