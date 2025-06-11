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

import sys
import argparse
import shutil
from pathlib import Path

import pandas as pd

import gemmi
from gemmi import cif

from xchemalign import dbreader, utils
from .utils import Constants


def run(base_dir: Path, input_path: Path, soakdb_file: Path, meta: dict, output_dir: Path, logger):
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
        name = row[Constants.SOAKDB_XTAL_NAME]
        xtal_out_dir = output_dir / name
        if xtal_out_dir.is_dir():
            # delete it
            shutil.rmtree(xtal_out_dir)
        # create the output dir
        xtal_out_dir.mkdir(parents=True)
        print(index, name, refinement_type)
        meta_data = crystals.get(name)
        if meta_data is None:
            logger.warn("Crystal {} not found in metadata. This is strange.".format(name))
        else:
            # event_map = ?? # how do we find the event map when sometimes there are multiple ones
            mtz_latest = row.get(Constants.SOAKDB_COL_MTZ_LATEST)
            mtz_free = row.get(Constants.SOAKDB_COL_MTZ_FREE)
            pdb = row.get(Constants.SOAKDB_COL_PDB)
            print("  mtz_latest: {}\n  mtz_free: {}\n  mmcif: {}\n  pdb: {}".format(mtz_latest, mtz_free, mmcif, pdb))

            if refinement_type == Constants.SOAKDB_VALUE_BUSTER:
                cif_doc = read_buster_structure(base_dir / utils.make_path_relative(Path(mmcif)))
            else:
                cif_doc = read_refmac_structure(base_dir / utils.make_path_relative(Path(pdb)))
            cif_doc.write_file(str(xtal_out_dir / 'structure.mmcif'))


def read_buster_structure(mmcif):
    doc = cif.read(str(mmcif))
    # is this a Structure or a Document?
    return doc


def read_refmac_structure(pdb):
    struc = gemmi.read_pdb(str(pdb))
    return struc.make_mmcif_document()


def main():
    # Example:
    #   python -m xchemalign.pdb_deposition -c config.yaml -o depo

    parser = argparse.ArgumentParser(description="copier")

    parser.add_argument("-c", "--config-file", required=True, help="config.yaml file")
    parser.add_argument("-m", "--meta-file", required=True, help="meta_aligner.yaml file")

    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    logger = utils.Logger(logfile=args.log_file, level=args.log_level)
    logger.info("pdb_deposition: ", args)
    utils.LOG = logger

    config = utils.read_config_file(args.config_file)
    meta = utils.read_config_file(args.meta_file)
    base_dir = config.get(Constants.CONFIG_BASE_DIR)

    inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
    logger.info("found {} inputs".format(len(inputs)))
    for input in inputs:
        if input[Constants.CONFIG_TYPE] == Constants.CONFIG_TYPE_MODEL_BUILDING:
            soakdb_file = utils.find_soakdb_file(input)
            input_path = input[Constants.CONFIG_DIR]
            logger.info("Running for " + input_path)
            run(Path(base_dir), Path(input_path), Path(soakdb_file), meta, Path(args.output_dir), logger)

    logger.report()
    logger.close()


if __name__ == "__main__":
    main()
