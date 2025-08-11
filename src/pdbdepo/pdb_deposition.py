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
import shutil
from pathlib import Path

import gemmi
from gemmi import cif

from xchemalign import dbreader, utils
from xchemalign.utils import Constants
from pdbdepo import merge_sf


def process_input(
    base_dir: Path,
    input_path: Path,
    collator_output_dir: Path,
    soakdb_file: Path,
    meta: dict,
    output_dir: Path,
    logger,
):
    path_to_soakdb_file = base_dir / input_path / soakdb_file
    crystals = meta[Constants.META_XTALS]
    logger.info("reading soakdb file:", path_to_soakdb_file)
    df = dbreader.read_pdb_depo(path_to_soakdb_file)
    logger.info("read {} rows".format(len(df)))
    for index, row in df.iterrows():
        mmcif = row[Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST]
        if mmcif:
            refinement_type = Constants.SOAKDB_VALUE_BUSTER
        else:
            refinement_type = Constants.SOAKDB_VALUE_REFMAC
        xtal_name = row[Constants.SOAKDB_XTAL_NAME]
        xtal_out_dir = output_dir / xtal_name
        if xtal_out_dir.is_dir():
            # delete it
            shutil.rmtree(xtal_out_dir)
        # create the output dir
        xtal_out_dir.mkdir(parents=True)
        print(index, xtal_name, refinement_type)
        meta_data = crystals.get(xtal_name)
        if meta_data is None:
            logger.warn("Crystal {} not found in metadata. This is strange.".format(xtal_name))
        else:
            mtz_latest = row.get(Constants.SOAKDB_COL_MTZ_LATEST)
            mtz_free = row.get(Constants.SOAKDB_COL_MTZ_FREE)
            pdb = row.get(Constants.SOAKDB_COL_PDB)
            crystallographic_files = meta_data.get(Constants.META_XTAL_FILES)
            ligand_binding_events = crystallographic_files.get(Constants.META_BINDING_EVENT)
            ccp4_files = []
            if ligand_binding_events:
                for ligand_binding_event in ligand_binding_events:
                    ccp4_filename = ligand_binding_event.get(Constants.META_FILE)
                    if not ccp4_filename:
                        print("WARNING: no event map file for " + xtal_name)
                    else:
                        ccp4_files.append(str(collator_output_dir.parent / ligand_binding_event[Constants.META_FILE]))

            mtz_free_path = None
            if mtz_free:
                # fix the path of the mtz_free file as it is probably just the filename, not the path
                mtz_free_path = Path(mtz_free)
                if not mtz_free_path.is_absolute():
                    parent = Path(pdb)
                    while True:
                        parent = parent.parent
                        if parent == Path('/'):
                            print("Couldn't find dir ")
                            break
                        if parent.name == xtal_name:
                            mtz_free_path = parent / mtz_free
                            break
                if not (base_dir / utils.make_path_relative(mtz_free_path)).is_file():
                    print("WARNING: couldn't find mtz_free file " + mtz_free)

            print(
                "  mtz_latest: {}\n  mtz_free: {}\n  mmcif: {}\n  pdb: {}\n  ccp4: {}".format(
                    mtz_latest, str(mtz_free_path), mmcif, pdb, ccp4_files
                )
            )

            if refinement_type == Constants.SOAKDB_VALUE_BUSTER:
                cif_doc = read_buster_structure(base_dir / utils.make_path_relative(Path(mmcif)))
            else:
                cif_doc = read_refmac_structure(base_dir / utils.make_path_relative(Path(pdb)))
            cif_doc.write_file(str(xtal_out_dir / 'structure.mmcif'))

            # handle the structure factors
            merge_sf.run(
                str(base_dir / utils.make_path_relative(Path(mtz_latest))),
                str(base_dir / utils.make_path_relative(mtz_free_path)),
                ccp4_files,
                str(xtal_out_dir / 'structure_factors.cif'),
            )


def read_buster_structure(mmcif):
    doc = cif.read(str(mmcif))
    # is this a Structure or a Document?
    return doc


def read_refmac_structure(pdb):
    struc = gemmi.read_pdb(str(pdb))
    return struc.make_mmcif_document()


def run(collator_dir, output_dir, logger):
    collator_path = Path(collator_dir)
    config = utils.read_config_file(collator_path / 'config.yaml')
    meta = utils.read_config_file(collator_path / 'meta_collator.yaml')
    base_dir = config.get(Constants.CONFIG_BASE_DIR)

    inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
    logger.info("found {} inputs".format(len(inputs)))
    for input in inputs:
        if input[Constants.CONFIG_TYPE] == Constants.CONFIG_TYPE_MODEL_BUILDING:
            soakdb_file = utils.find_soakdb_file(input)
            input_path = input[Constants.CONFIG_DIR]
            logger.info("Running for " + input_path)
            process_input(
                Path(base_dir), Path(input_path), collator_path, Path(soakdb_file), meta, Path(output_dir), logger
            )


def main():
    # Example:
    #   python -m pdbdepo.pdb_deposition -w path_to_collator_output -o depo

    parser = argparse.ArgumentParser(description="pdb deposition")

    parser.add_argument("-w", "--collator-dir", required=True, help="collator's output dir")

    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("pdb_deposition: ", args)
    utils.LOG = logger

    run(args.collator_dir, args.output_dir, logger)

    logger.report()
    logger.close()


if __name__ == "__main__":
    main()
