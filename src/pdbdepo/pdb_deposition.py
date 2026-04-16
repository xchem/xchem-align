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
import bz2
import datetime
import glob
import re
import shutil
from pathlib import Path
from collections import OrderedDict

import gemmi
from gemmi import cif

from mmcif_gen.facilities import xchem

from rdkit import rdBase
from rdkit import Chem

from xchemalign import dbreader, utils
from xchemalign.utils import Constants
from pdbdepo import merge_sf
from pdbdepo import scrape_processing_stats

LOG = utils.Logger()


def info(*args, **kwargs):
    if LOG:
        LOG.info(*args, **kwargs)


def warn(*args, **kwargs):
    if LOG:
        LOG.warn(*args, **kwargs)


def error(*args, **kwargs):
    if LOG:
        LOG.error(*args, **kwargs)


def run_mmcifgen(inv_id, soakdb_sqlite, metadata_csv, out_dir):
    p = Path(out_dir)
    if not p.exists():
        p.mkdir(parents=True)

    xchem.run(inv_id, soakdb_sqlite, metadata_csv, out_dir, 'xchem_operations.json')


def generate_dimple_log_path(xtal_dir_path):
    dimple_log_path = xtal_dir_path / 'dimple/dimple/dimple.log'
    if dimple_log_path.is_file():
        return dimple_log_path
    else:
        warn('Dimple log file not found: ' + str(dimple_log_path))
        return None


def generate_collection_info_path(xtal_dir_path, xtal_name):
    collection_info_path = xtal_dir_path / 'autoprocessing/' / (xtal_name + '_collection_info.cif')
    if collection_info_path.is_file():
        return collection_info_path
    else:
        warn('Collection info CIF file not found: ' + str(collection_info_path))
        return None


def read_software_templates():
    d = {}
    templates = glob.glob('config/pdb-depo/*.cif')
    for t in templates:
        c = cif.read(t)
        n = t[t.rfind('/') + 1 : -4].lower()
        d[n] = c
    info('read software templates', d)
    return d


def process_input(
    base_dir: Path,
    input_path: Path,
    collator_output_dir: Path,
    soakdb_file_p: Path,
    input_config: dict,
    meta_collator: dict,
    mmcifgen_block,
    output_dir: Path,
    cmpd_codes_dict: dict = {},
    debug=False,
):
    software_templates = read_software_templates()

    crystals = meta_collator[Constants.META_XTALS]

    # read the sequences
    default_seq, variants = utils.read_sequences(base_dir, input_config)

    xtals_list = []
    cmpd_moldata = {}
    # grab the ligand smiles of the collator output
    for crystal in crystals:
        mol_data = []
        cmpd_moldata[crystal] = mol_data
        ligand_cif = crystals[crystal][Constants.META_XTAL_FILES].get(Constants.META_XTAL_CIF)
        if ligand_cif:
            ligands = ligand_cif[Constants.META_LIGANDS]
            if ligands:
                for lig_name in ligands:
                    data = ligands[lig_name]
                    cmpd_code = data.get(Constants.META_CMPD_CODE)
                    cmpd_smiles = data.get(Constants.META_SMILES)
                    mol_data.append((cmpd_code, cmpd_smiles))
    info('read ligand data for', len(cmpd_moldata), 'ligands')

    mmcifgen_refine_tags = None
    mmcifgen_refine_values = None
    mmcifgen_diffrn_tags = None
    mmcifgen_diffrn_values = None
    mmcifgen_refine_item = find_loop_item(mmcifgen_block, '_refine')
    mmcifgen_diffrn_item = find_loop_item(mmcifgen_block, '_diffrn')
    if mmcifgen_refine_item:
        mmcifgen_refine_tags = list(mmcifgen_refine_item.loop.tags)
        mmcifgen_refine_values = list(mmcifgen_refine_item.loop.values)
        mmcifgen_refine_item.erase()
    if mmcifgen_diffrn_item:
        mmcifgen_diffrn_tags = list(mmcifgen_diffrn_item.loop.tags)
        mmcifgen_diffrn_values = list(mmcifgen_diffrn_item.loop.values)
        mmcifgen_diffrn_item.erase()

    (mmcif_gen_entity_tags, mmcif_gen_entity_values) = read_mmcifgen_entity_data_and_erase(mmcifgen_block)

    info("reading soakdb file:", soakdb_file_p)
    df = dbreader.read_pdb_depo(soakdb_file_p)
    info("read {} rows".format(len(df)))
    for index, row in df.iterrows():
        xtal_name = row[Constants.SOAKDB_XTAL_NAME]
        seq_dict = variants.get(xtal_name, default_seq)
        # info('sequences:', seq_dict)
        cmpd_codes = []
        xtals_list.append(xtal_name)
        mmcif = row[Constants.SOAKDB_COL_REFINEMENT_MMCIF_MODEL_LATEST]
        cmpd_code = row[Constants.SOAKDB_COL_COMPOUND_CODE]
        tokens = cmpd_code.split(';')
        for token in tokens:
            cmpd_codes.append(token.strip())
        if mmcif is None or mmcif == 'None':
            refinement_type = Constants.SOAKDB_VALUE_REFMAC
            refinement_prog = 'refmac'
        else:
            refinement_type = Constants.SOAKDB_VALUE_BUSTER
            refinement_prog = 'buster'

        info(index, xtal_name, refinement_type)

        data_processing_logfile = row[Constants.SOAKDB_COL_DATA_PROCESSING_PATH_TO_LOGFILE]
        aimless_ver = None
        data_processing_logfile_p = None
        if data_processing_logfile and data_processing_logfile != 'None':
            data_processing_logfile_p = Path(data_processing_logfile)
            aimless = base_dir / utils.make_path_relative(data_processing_logfile_p).parent / 'aimless.log'
            if aimless.is_file():
                aimless_ver = scrape_aimless_version(aimless)
                info('AIMLESS+VER', aimless_ver)

        xtal_in_path = base_dir / input_path / utils.Constants.DEFAULT_MODEL_BUILDING_DIR / xtal_name
        xtal_out_path = output_dir / xtal_name
        if xtal_out_path.is_dir():
            # delete it
            shutil.rmtree(xtal_out_path)
        # create the output dir
        xtal_out_path.mkdir(parents=True)

        meta_data = crystals.get(xtal_name)
        if meta_data is None:
            warn("Crystal {} not found in metadata. This is strange.".format(xtal_name))
        else:
            mtz_latest = row.get(Constants.SOAKDB_COL_MTZ_LATEST)
            mtz_free = row.get(Constants.SOAKDB_COL_MTZ_FREE)
            pdb = row.get(Constants.SOAKDB_COL_PDB)
            crystallographic_files = meta_data.get(Constants.META_XTAL_FILES)
            ligand_binding_events = crystallographic_files.get(Constants.META_BINDING_EVENT)
            ccp4_files = []
            ccp4_sources = []
            if ligand_binding_events:
                for ligand_binding_event in ligand_binding_events:
                    ccp4_filename = ligand_binding_event.get(Constants.META_FILE)
                    source_filename = ligand_binding_event.get(Constants.META_SOURCE_FILE)
                    if not ccp4_filename:
                        info("no event map file for " + xtal_name)
                    else:
                        ccp4_files.append(str(collator_output_dir.parent / ligand_binding_event[Constants.META_FILE]))
                        ccp4_sources.append('/' + str(Path(source_filename).relative_to(base_dir)))

            mtz_free_path = None
            if mtz_free:
                # fix the path of the mtz_free file as it is probably just the filename, not the path
                mtz_free_path = Path(mtz_free)
                if not mtz_free_path.is_absolute():
                    parent = Path(pdb)
                    while True:
                        parent = parent.parent
                        if parent == Path('/'):
                            info("Couldn't find dir ")
                            break
                        if parent.name == xtal_name:
                            mtz_free_path = parent / mtz_free
                            break
                if not (base_dir / utils.make_path_relative(mtz_free_path)).is_file():
                    warn("couldn't find mtz_free file " + mtz_free)

            info(
                "located files:\n  mtz_latest: {}\n  mtz_free: {}\n  mmcif: {}\n  pdb: {}\n  ccp4: {}".format(
                    mtz_latest, str(mtz_free_path), mmcif, pdb, ccp4_files
                )
            )

            if refinement_type == Constants.SOAKDB_VALUE_BUSTER:
                structure_cif_doc = read_buster_structure(base_dir / utils.make_path_relative(Path(mmcif)), seq_dict)
            else:
                structure_cif_doc = read_refmac_structure(base_dir / utils.make_path_relative(Path(pdb)), seq_dict)

            structure_cif_block0 = structure_cif_doc[0]
            if debug:
                structure_cif_doc.write_file(str(xtal_out_path / 'original.cif'))

            # find and delete the _exptl loop as mmcifgen handles this
            delete_pair_item(structure_cif_doc, '_exptl')

            merge_entity(xtal_name, structure_cif_block0, mmcif_gen_entity_tags, mmcif_gen_entity_values)

            # do some cleaning
            delete_pair_item(structure_cif_doc, '_struct')
            delete_pair_item(structure_cif_doc, '_struct_keywords')
            delete_pair_item(structure_cif_doc, '_pdbx_database_status')
            delete_loop_item(structure_cif_doc, '_exptl')
            delete_loop_item(structure_cif_doc, '_reflns')

            refine_item = find_loop_item(structure_cif_block0, '_refine')
            if refine_item:
                # probably a loop if generated from refmac pdb file
                model_refine_tags = list(refine_item.loop.tags)
                model_refine_values = list(refine_item.loop.values)
                refine_item.erase()
                append_loop_items(
                    structure_cif_block0,
                    model_refine_tags,
                    model_refine_values,
                    mmcifgen_refine_tags,
                    mmcifgen_refine_values,
                )
            else:
                # might be pairs if generated by buster
                d = delete_pairs(structure_cif_block0, '_refine')
                if d:
                    model_refine_tags, model_refine_values = dict_to_tags_values(d)
                    append_loop_items(
                        structure_cif_block0,
                        model_refine_tags,
                        model_refine_values,
                        mmcifgen_refine_tags,
                        mmcifgen_refine_values,
                    )

            # add in the common metadata (generated by mmcif-gen)
            for item in mmcifgen_block:
                to_add = True
                if item.loop is not None:
                    if '_struct.title' in item.loop.tags:
                        # special case of the structure title loop that needs expanding
                        to_add = False
                        if len(item.loop.values) > 1:
                            template = item.loop.values[1]
                            template = template.replace('$CompoundCode', cmpd_code)
                            template = template.replace('$CrystalName', xtal_name)
                            cmpd_codes = cmpd_codes_dict.get(xtal_name)
                            if cmpd_codes:
                                for i, cmpd_code in enumerate(cmpd_codes):
                                    template = template.replace('$ExternalCode' + str(i + 2), cmpd_code)
                            # erase any non-substituted ExternalCodes
                            for i in range(1, 9):
                                template = template.replace('$ExternalCode' + str(i), '')
                            new_loop = structure_cif_block0.init_loop('', item.loop.tags)
                            new_loop.add_row([item.loop.values[0], template])
                if to_add:
                    # otherwise add the item as it is
                    added_item = structure_cif_block0.add_item(item)

            data_processing_log_file = None
            data_processing_prog = row.get(Constants.SOAKDB_COL_DATA_PROCESSING_PROGRAM)
            if data_processing_prog == 'dials':
                data_processing_prog = 'xia2-dials'
            if data_processing_prog and data_processing_prog != 'None':
                data_processing_prog = data_processing_prog.lower()
                info('data processing was done with ' + data_processing_prog)
                if data_processing_logfile_p:
                    stats_cif = (
                        base_dir / utils.make_path_relative(Path(data_processing_logfile).parent) / 'xia2.mmcif.bz2'
                    )
                    if not stats_cif.is_file():
                        p = None
                        # info(str(stats_cif) + ' mmcif file containing data processing stats not found')
                        base_data_processing_logfile_p = base_dir / utils.make_path_relative(data_processing_logfile_p)
                        if base_data_processing_logfile_p.is_file():
                            p = base_data_processing_logfile_p
                            info('processing stats from', data_processing_logfile_p.name)
                            data_processing_log_file = p
                        elif data_processing_prog == 'autoproc' and data_processing_logfile.endswith('.log'):
                            # look for the aimless.log file
                            if aimless.is_file():
                                p = aimless
                                info('processing stats from aimless.log')
                                data_processing_log_file = aimless
                            else:
                                warn('no aimless.log file found')
                        else:
                            warn('could not find logfile', data_processing_logfile)

                        if p:
                            stats_doc = scrape_processing_stats.handle_file(p, data_processing_prog, None, None)
                            if stats_doc is None:
                                warn('stats could not be scraped from log file')
                            else:
                                stats_block0 = stats_doc[0]
                                for item in stats_block0:
                                    structure_cif_block0.add_item(item)

                    else:
                        info('processing stats from xia2.mmcif.bz2')
                        data_processing_log_file = stats_cif
                        stats_item_reflns_shell = None
                        stats_pair_reflns = OrderedDict()
                        stats_pair_diffrn = OrderedDict()
                        with bz2.open(stats_cif) as stats:
                            txt = stats.read()
                            doc = cif.read_string(txt)

                            # grab the reflns and diffrn data from second block
                            for item in doc[1]:
                                if item.loop is None:
                                    p = item.pair
                                    if p[0].startswith('_reflns.'):
                                        stats_pair_reflns[p[0]] = p[1]
                                    if p[0].startswith('_diffrn.'):
                                        stats_pair_diffrn[p[0]] = p[1]
                                else:
                                    for tag in item.loop.tags:
                                        if tag.startswith('_reflns_shell.'):
                                            stats_item_reflns_shell = item

                        if stats_pair_reflns:
                            # completely weird gemmi quirk results in pairs being added to their old locations when
                            # added individually, but not when added in one go
                            d = {}
                            for k, v in stats_pair_reflns.items():
                                d[k[8:]] = v
                            structure_cif_block0.set_pairs('_reflns.', d)
                        if stats_item_reflns_shell:
                            structure_cif_block0.add_item(stats_item_reflns_shell)

            # handle the collection_info.cif file
            collection_info_p = generate_collection_info_path(xtal_in_path, xtal_name)
            collection_info_diffrn_item = None
            if collection_info_p is not None:
                doc = cif.read(str(collection_info_p))
                info('read ' + str(collection_info_p))
                for item in doc[0]:
                    if item.loop:
                        if item.loop.tags[0].startswith('_diffrn.'):
                            # don't add this one as it needs special handling
                            collection_info_diffrn_item = item
                            continue
                    structure_cif_block0.add_item(item)
            combine_diffrn_loops(
                collection_info_diffrn_item, mmcifgen_diffrn_tags, mmcifgen_diffrn_values, structure_cif_block0
            )

            # add in the _software section
            add_software_loop(software_templates, structure_cif_block0, refinement_prog, data_processing_prog)

            # add_sequences(xtal_name, seq_dict, structure_cif_block0)

            # move coordinates to the end
            coordinates_item = find_loop_item(structure_cif_block0, '_atom_site')
            item_count = 0
            coordinates_pos = None
            for item in structure_cif_block0:
                if item == coordinates_item:
                    coordinates_pos = item_count
                item_count += 1
            if coordinates_pos is not None:
                structure_cif_block0.move_item(coordinates_pos, item_count - 1)

            structure_cif_doc.write_file(str(xtal_out_path / (xtal_name + '_struc.cif')))

            # include the ligand CIF file
            cif_file = row[Constants.SOAKDB_COL_CIF]
            if cif_file:
                p = base_dir / utils.make_path_relative(Path(cif_file))
                if p.is_file():
                    p2 = shutil.copy2(p, xtal_out_path / (xtal_name + '_lig.cif'), follow_symlinks=True)
                    info('copied ligand CIF', p)
                else:
                    warn('ligand CIF file not found:', p)

            if debug:
                if data_processing_log_file:
                    p2 = shutil.copy2(data_processing_log_file, xtal_out_path, follow_symlinks=True)
                    info('copied log file', data_processing_log_file)

            merge_sf.run(
                str(base_dir / utils.make_path_relative(Path(mtz_latest))),
                str(base_dir / utils.make_path_relative(mtz_free_path)),
                ccp4_files,
                ccp4_sources,
                str(xtal_out_path / (xtal_name + '_sf.cif')),
                mtz_latest,
                str(mtz_free_path),
                output_individual=debug,
            )

    # write out the ligands to a file
    with open(output_dir / 'ligands.tab', 'wt') as tsv:
        for xtal in xtals_list:
            data = cmpd_moldata.get(xtal)
            if not data:
                warn('compound data not found for', xtal)
            else:
                values = []
                for t in data:
                    values.append(xtal)
                    values.append(t[0])
                    values.append(t[1])
                    mol = Chem.MolFromSmiles(t[1])
                    inchis = Chem.MolToInchi(mol)
                    inchik = Chem.InchiToInchiKey(inchis)
                    values.append(inchis)
                    values.append(inchik)
                tsv.write('\t'.join(values) + '\n')


def add_software_loop(templates_dict, block, refinement_prog, data_processing_prog):
    # TODO - get the tags from the templates rather than hardcoding them
    loop = block.init_loop(
        '_software.',
        [
            'pdbx_ordinal',
            'name',
            'classification',
            'type',
            'version',
            'date',
            'location',
            'description',
            'pdbx_reference_DOI',
        ],
    )

    is_error = False
    refinement_t = templates_dict.get(refinement_prog.lower())
    if not refinement_t:
        error('no software template found for', refinement_prog.lower())
        is_error = True

    data_processing_t = templates_dict.get(data_processing_prog.lower())
    if not data_processing_t:
        error('no software template found for', data_processing_prog.lower())
        is_error = True

    refinement_item = refinement_t[0].find_loop_item('_software.pdbx_ordinal')
    if not refinement_item or not refinement_item.loop:
        error('software template found for', refinement_prog.lower(), 'is not valid. No _software loop found.')
        is_error = True

    data_processing_item = data_processing_t[0].find_loop_item('_software.pdbx_ordinal')
    if not data_processing_item or not data_processing_item.loop:
        error('software template found for', data_processing_prog.lower(), 'is not valid. No _software loop found.')
        is_error = True

    if is_error:
        exit(1)

    refinement_tags = refinement_item.loop.tags
    refinement_values = refinement_item.loop.values
    data_processing_tags = data_processing_item.loop.tags
    data_processing_values = data_processing_item.loop.values

    values = []
    i = 1
    # print(data_processing_tags)
    for j, value in enumerate(data_processing_values):
        if j % len(data_processing_tags) == 0:
            # print(i, j, values)
            if values:
                loop.add_row(values)
                i += 1
            values.clear()
            values.append(str(i))
        else:
            values.append(value)

    loop.add_row(values)
    i += 1
    values.clear()

    # print(refinement_tags)
    for j, value in enumerate(refinement_values):
        if j % len(refinement_tags) == 0:
            # print(i, j, values)
            if values:
                loop.add_row(values)
                i += 1
            values.clear()
            values.append(str(i))
        else:
            values.append(value)
        i += 1
    loop.add_row(values)


def scrape_aimless_version(aimless_p):
    with open(aimless_p, "rt") as f:
        i = 0
        for line in f:
            i += 1
            if i > 10:
                info('could not read aimless version from', aimless_p)
                return None
            # the line we are after looks like this:
            ### CCP4 8.0.016: AIMLESS           version 0.7.13 : 20/07/23##
            match = re.search(r'.*AIMLESS\s+version\s+([\d\.]+)\s+', line)
            if match:
                aimless_ver = match.group(1)
                LOG.info('found aimless version', aimless_ver)
                return aimless_ver
    return None


def dict_to_tags_values(d: dict):
    """
    Convert a dict to lists of tags and values suitable for a loop
    :param d:
    :return:
    """
    tags = []
    values = []
    for t, v in d.items():
        tags.append(t)
        values.append(v)
    return tags, values


def read_mmcifgen_entity_data_and_erase(mmcifgen_block):
    mmcifgen_entity_item = find_loop_item(mmcifgen_block, '_entity')

    mmcif_gen_entity_tags = None
    mmcif_gen_entity_values = None

    if mmcifgen_entity_item.loop:
        mmcif_gen_entity_tags = mmcifgen_entity_item.loop.tags
        mmcif_gen_entity_values = mmcifgen_entity_item.loop.values
    if mmcifgen_entity_item:
        mmcifgen_entity_item.erase()

    return mmcif_gen_entity_tags, mmcif_gen_entity_values


def merge_entity(xtal_name, model_block, mmcif_gen_entity_tags, mmcif_gen_entity_values):
    model_entity_item = find_loop_item(model_block, '_entity')

    model_entity_tags = None
    model_entity_values = None

    if model_entity_item and model_entity_item.loop:
        model_entity_tags = model_entity_item.loop.tags
        model_entity_values = model_entity_item.loop.values
        model_entity_item.erase()
    else:
        warn('_entity loop not present in model CIF for', xtal_name)

    all_entity_tags = []
    if mmcif_gen_entity_tags:
        all_entity_tags.extend(mmcif_gen_entity_tags)
    if model_entity_tags:
        for t in model_entity_tags:
            if t not in all_entity_tags:
                all_entity_tags.append(t)

    loop = model_block.init_loop('', all_entity_tags)
    model_entity_rows = collect_entity_values(
        model_entity_tags, model_entity_values, exclude_pairs=[('_entity.type', 'polymer')]
    )
    mmcifgen_entity_rows = collect_entity_values(mmcif_gen_entity_tags, mmcif_gen_entity_values)
    append_entity_values(all_entity_tags, mmcifgen_entity_rows, loop)
    append_entity_values(all_entity_tags, model_entity_rows, loop)


def append_entity_values(tags, rows, loop):
    for d in rows:
        data = []
        for tag in tags:
            if tag in d:
                data.append(d[tag])
            else:
                data.append('?')
        loop.add_row(data)


def collect_entity_values(tags, values, exclude_pairs=[]):
    rows = []
    if values:
        d = None
        for i, value in enumerate(values):
            if i % len(tags) == 0:
                if d:
                    add_row_if_not_excluded(d, rows, exclude_pairs)
                d = {}
            tag = tags[i % len(tags)]
            d[tag] = value
        add_row_if_not_excluded(d, rows, exclude_pairs)
    return rows


def add_row_if_not_excluded(d, rows, exclude_pairs):
    is_excluded = False
    for exclude in exclude_pairs:
        if d.get(exclude[0]) == exclude[1]:
            is_excluded = True
    if not is_excluded:
        rows.append(d)


def delete_pair_item(doc, prefix):
    for block in doc:
        for item in block:
            if item.pair is not None:
                if item.pair[0].startswith(prefix + '.'):
                    item.erase()


def delete_loop_item(doc, prefix):
    for block in doc:
        item = find_loop_item(block, prefix)
        if item:
            # print('erasing', item)
            item.erase()


def find_loop_item(block, prefix):
    for item in block:
        if item.loop is not None:
            for tag in item.loop.tags:
                if tag.startswith(prefix + '.'):
                    return item
    return None


def delete_pairs(block, prefix):
    """
    Delete the pairs with this prefix and generate a dict of tags and values
    :param block:
    :param prefix:
    :return:
    """
    d = {}
    for item in block:
        if item.pair is not None:
            if item.pair[0].startswith(prefix + '.'):
                d[item.pair[0]] = item.pair[1]
                item.erase()
    return d


def append_loop_items(block_to_add_to: cif.Block, tags1, values1, tags2, values2):
    """Simplistic approach that merges two 1-row loops, add the data from item2 into item1

    :param block_to_add_to: the block the loop needs to be added to
    :param tags1: the first tags to append
    :param values1: the first values to append
    :param tags1: the second tags to append
    :param values1: the second values to append
    :return:
    """

    # its seems that you can't just append to the existing loop's tags and value, you have to create a new loop
    combined_tags = []
    combined_values = []
    combined_tags.extend(tags1)
    combined_tags.extend(tags2)
    combined_values.extend(values1)
    combined_values.extend(values2)

    loop = block_to_add_to.init_loop('', combined_tags)

    loop.add_row(combined_values)


def append_data_to_values(tags1, values1, tags2, values2, ordinal_col):
    """
    Append the data is tags/values2 to that in tags/values1, ensuring the order is correct
    If no value is found then ? is written
    :param tags1:
    :param values1:
    :param tags2:
    :param values2:
    :return: A list of lists of the values
    """

    if ordinal_col:
        ordinal_idx = tags1.index(ordinal_col)
    else:
        ordinal_idx = None

    combined_values = []
    combined_values.append(values1)
    n = len(tags2)
    dicts = []
    d = None
    ordinal_count = len(combined_values)
    for i, val in enumerate(values2):
        ordinal_count += 1
        if i % n == 0:
            d = {}
            dicts.append(d)
        if i % n == ordinal_idx:
            d[tags2[i % n]] = str(ordinal_count)
        else:
            d[tags2[i % n]] = val

    for d in dicts:
        l = []
        combined_values.append(l)
        for tag1 in tags1:
            if tag1 in d:
                l.append(d[tag1])
            else:
                l.append('?')

    return tags1, combined_values


def combine_diffrn_loops(collection_info_item: cif.Item, tags_list: list, values_list: list, into: cif.Block):
    tags = []
    if collection_info_item:
        tags.extend(collection_info_item.loop.tags)
    if tags_list:
        tags.extend(tags_list)
    values = []
    if collection_info_item:
        values.extend(collection_info_item.loop.values)
    if values_list:
        values.extend(values_list)
    if tags:
        new_loop = into.init_loop('', tags)
        new_loop.add_row(values)


def prune_loop(loop, retain):
    for tag in loop.tags:
        if not tag.startswith(retain):
            loop.remove_column(tag)


def read_buster_structure(mmcif_file, seq_dict):
    struc = gemmi.read_structure(str(mmcif_file))
    adjust_chains_and_entities(struc, seq_dict)
    return struc.make_mmcif_document()


def read_refmac_structure(pdb_file, seq_dict):
    struc = gemmi.read_pdb(str(pdb_file), ignore_ter=True)
    struc.setup_entities()
    adjust_chains_and_entities(struc, seq_dict)
    return struc.make_mmcif_document()


def adjust_chains_and_entities(struc, seq_dict):
    """

    :param struc: the gemmi Structure
    :param seq_dict: dict keyed by chain name, values a tuple of (entity_name, entity_seq)
    :return:
    """

    entities_dict = {}
    for ety in struc.entities:
        for chain_name in seq_dict:
            t = seq_dict[chain_name]
            if t[0] == ety.name:
                entities_dict[ety.name] = ety
                info('setting sequence of entity', ety.name, 'to', abbreviate_sequence(t[1]))
                ety.full_sequence = gemmi.expand_one_letter_sequence(t[1], gemmi.ResidueKind.AA)
                break

    for chain_name in seq_dict:
        t = seq_dict[chain_name]
        ety_name = t[0]
        for chain in struc[0]:  # always only a single model
            if chain_name == chain.name:
                for subchain in chain.subchains():
                    if subchain.check_polymer_type() == gemmi.PolymerType.PeptideL:
                        subchain_id = subchain.subchain_id()
                        # info('chain info', ety_name, chain.name, subchain_id, subchain.check_polymer_type())
                        ety = entities_dict.get(ety_name)
                        if not ety:
                            error('entity/chain mismatch. Looking for', ety_name, 'among', entities_dict.keys())
                            exit(1)
                        # info('entity', ety.name, 'has subchains', ety.subchains)
                        if not subchain_id in ety.subchains:
                            info('adding subchain', subchain_id, 'to entity', ety.name)
                            subchains = ety.subchains
                            subchains.append(subchain_id)  # this also removes it from its current entity
                            ety.subchains = subchains

    etys_to_remove = []  # indexes of entities that have to be removed
    for i, ety in enumerate(struc.entities):
        if len(ety.subchains) == 0:
            etys_to_remove.append(i)

    if etys_to_remove:
        etys_to_remove.sort(reverse=True)
        for ety_idx in etys_to_remove:
            info('deleting entity', ety_idx, struc.entities[ety_idx].name)
            del struc.entities[ety_idx]

    struc.setup_entities()
    struc.add_entity_ids(overwrite=True)
    struc.add_entity_types(overwrite=True)
    struc.assign_subchains(force=True)


def read_phasing_software(xtal_dir):
    dimple_log_file_path = Path(xtal_dir) / 'dimple/dimple/dimple.log'
    dimple_ver = None
    if dimple_log_file_path.is_file():
        with open(dimple_log_file_path, 'rt') as dimple:
            looking = False
            for line in dimple:
                line = line.strip()
                if looking:
                    match = re.search(r'^version:\s+([\d\.]+)$', line)
                    if match:
                        LOG.info('found dimple version', match.group(1))
                        dimple_ver = match.group(1)
                        break
                    match = re.search(r'^\[\w+\]$', line)
                    if match:
                        LOG.info('dimple version not found')
                        break
                else:
                    if line == '[workflow]':
                        looking = True

    else:
        LOG.warn('dimple log file not found:', str(dimple_log_file_path))

    return dimple_ver


def create_pairs(data: dict, prefix: str, block: cif.Block):
    block.set_pairs(prefix, data)


def create_loop(data: dict, prefix: str, block: cif.Block):
    loop = block.init_loop(prefix, list(data.keys()))
    size = len(next(iter(data.values())))
    for i in range(size):
        values = []
        for k, v in data.items():
            values.append(str(v[i]))
        loop.add_row(values)


def read_cmpd_codes(filename):
    d = {}
    if filename:
        with open(filename, 'rt') as cmpd_codes:
            line = cmpd_codes.readline()  # header line
            for line in cmpd_codes:
                tokens = line.split(',')
                trimmed = [t.strip() for t in tokens]
                d[trimmed[0]] = trimmed[1:]
    info('read compound codes for', len(d), 'crystals')
    return d


def run(collator_dir, metadata_csv, output_dir, compound_codes_csv=None, debug=False):
    info('run on ' + str(datetime.datetime.now()))
    info('using RDKit version ' + rdBase.rdkitVersion)
    # info('using InCHI version ' + Chem.GetInchiVersion())

    cmpd_codes_dict = read_cmpd_codes(compound_codes_csv)

    collator_path = Path(collator_dir)
    config = utils.read_config_file(collator_path / 'config.yaml')
    meta_collator = utils.read_config_file(collator_path / 'meta_collator.yaml')
    base_dir = config.get(Constants.CONFIG_BASE_DIR)
    base_dir_p = Path(base_dir)
    info('using base dir of', base_dir)

    inputs = utils.find_property(config, Constants.CONFIG_INPUTS)
    info("found {} inputs".format(len(inputs)))
    for i, input in enumerate(inputs):
        if input[Constants.CONFIG_TYPE] == Constants.CONFIG_TYPE_MODEL_BUILDING:
            input_path = input[Constants.CONFIG_DIR]
            soakdb_prop = utils.find_soakdb_file(input)
            print('PROP', soakdb_prop)
            soakdb_file_p = base_dir_p / input_path / soakdb_prop

            print('SOAKDB P', str(soakdb_file_p))

            # run mmcif-gen
            input_name = 'input' + str(i + 1)
            run_mmcifgen(input_name, str(soakdb_file_p), metadata_csv, output_dir)

            output_p = Path(output_dir)

            meta_doc = cif.read(str(output_p / (input_name + '_model.cif')))
            meta_mmcifgen = meta_doc[0]

            info("running for input " + input_path)
            process_input(
                base_dir_p,
                Path(input_path),
                collator_path,
                soakdb_file_p,
                input,
                meta_collator,
                meta_mmcifgen,
                Path(output_dir),
                cmpd_codes_dict=cmpd_codes_dict,
                debug=debug,
            )


def abbreviate_sequence(seq):
    if len(seq) > 25:
        return seq[:10] + '...' + seq[-10:] + ' (' + str(len(seq)) + ' residues)'
    else:
        return seq


def main():
    # Example:
    #   python -m pdbdepo.pdb_deposition -w path_to_collator_output -o pdb_depo -m metadata.csv

    global LOG

    parser = argparse.ArgumentParser(description="pdb deposition")

    parser.add_argument("-w", "--collator-dir", required=True, help="collator's output dir")
    parser.add_argument("-m", "--metadata-csv", required=True, help="CSV file with common metadata (for mmcif-gen)")
    parser.add_argument("-c", "--compound-codes-csv", help="CSV file with compound codes for the title")

    parser.add_argument("-o", "--output-dir", required=True, help="Output directory")
    parser.add_argument("-d", "--debug", action="store_true", help="Include source files in output")
    parser.add_argument("-l", "--log-file", help="File to write logs to")
    parser.add_argument("--log-level", type=int, default=0, help="Logging level (0=INFO, 1=WARN, 2=ERROR)")

    args = parser.parse_args()
    if args.log_file:
        logfile_p = Path(args.log_file)
    else:
        logfile_p = Path(args.output_dir) / 'pdb_depo.log'
    if not logfile_p.parent.is_dir():
        logfile_p.parent.mkdir(parents=True)

    LOG = utils.Logger(logfile=str(logfile_p), level=args.log_level)

    utils.LOG = LOG
    scrape_processing_stats.LOG = LOG
    merge_sf.LOG = LOG
    LOG.info("pdb_deposition: ", args)

    run(
        args.collator_dir,
        args.metadata_csv,
        args.output_dir,
        compound_codes_csv=args.compound_codes_csv,
        debug=args.debug,
    )

    LOG.report()
    LOG.close()


if __name__ == "__main__":
    main()
