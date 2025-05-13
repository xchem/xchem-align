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
import datetime

import gemmi
from gemmi import cif


def run(latest_mtz, free_mtz, event_map, output):
    today = datetime.date.today()
    formatted_date = today.strftime("%Y-%m-%d")

    doc_final = cif.Document()

    # handle the latest_mtz data
    cif_s = read_mtz(latest_mtz)
    mtz_doc = cif.read_string(cif_s)
    mtz_blk = mtz_doc.sole_block()
    mtz_blk.name = "rxxxxsf"
    mtz_blk.set_pair("_audit.revision_id", "1_0")
    mtz_blk.set_pair("_audit.creation_date", formatted_date)
    mtz_blk.set_pair("_audit.update_record", "Initial release")
    mtz_blk.set_pair("_diffrn.id", "1")
    mtz_blk.set_pair("_diffrn.details", "data from final refinement with ligand, final.mtz")
    doc_final.add_copied_block(mtz_blk)

    # handle the free_mtz data
    cif_s = read_mtz(free_mtz)
    mtz_doc = cif.read_string(cif_s)
    mtz_blk = mtz_doc.sole_block()
    mtz_blk.name = "rxxxxAsf"
    mtz_blk.set_pair("_diffrn.id", "1")
    mtz_blk.set_pair("_diffrn.details", "data from original reflections, data.mtz")
    doc_final.add_copied_block(mtz_blk)

    # handle the event map data
    cif_s = read_ccp4(event_map)
    mtz_doc = cif.read_string(cif_s)
    mtz_blk = mtz_doc.sole_block()
    mtz_blk.name = "rxxxxBsf"
    mtz_blk.set_pair("_diffrn.id", "1")
    mtz_blk.set_pair("_diffrn.details", "data for ligand evidence map (PanDDA event map), event_map_1.mtz")
    mtz_blk.set_pair("_diffrn_radiation_wavelength.id", "1")
    mtz_blk.set_pair("_diffrn_radiation_wavelength.wavelength", "0.92209")

    doc_final.add_copied_block(mtz_blk)

    doc_final.write_file(output)


def read_mtz(file):
    mtz = gemmi.read_mtz_file(file)
    mtz.title = "MMMM"
    print("read mtz")
    to_cif = gemmi.MtzToCif()
    cif_s = to_cif.write_cif_to_string(mtz)
    return cif_s


def read_ccp4(file):
    # see https://gemmi.readthedocs.io/en/latest/hkl.html#id7

    RESOLUTION_LIMIT = 1.5

    map = gemmi.read_ccp4_map(file)
    print("read ccp4")
    sf = gemmi.transform_map_to_f_phi(map.grid, half_l=True)
    data = sf.prepare_asu_data(dmin=RESOLUTION_LIMIT)

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = sf.spacegroup
    mtz.set_cell_for_all(sf.unit_cell)
    mtz.add_dataset("unknown")
    mtz.add_column("FWT", "F")
    mtz.add_column("PHWT", "P")
    mtz.set_data(data)

    to_cif = gemmi.MtzToCif()
    cif_s = to_cif.write_cif_to_string(mtz)
    return cif_s


def main():
    parser = argparse.ArgumentParser(description="merge_sf")

    parser.add_argument("-l", "--latest-mtz", required=True, help="MTZ latest file to merge")
    parser.add_argument("-f", "--free-mtz", required=True, help="MTZ free file to merge")
    parser.add_argument("-e", "--event-map", required=True, help="Panddas event map file to merge")
    parser.add_argument("-o", "--output", required=True, help="Output CIF file")

    parser.add_argument("--log-level", type=int, default=0, help="Logging level")

    args = parser.parse_args()

    run(args.latest_mtz, args.free_mtz, args.event_map, args.output)


if __name__ == "__main__":
    main()
