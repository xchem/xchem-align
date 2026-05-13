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

"""
Generate a draft assemblies.yaml by inspecting the reference PDB files listed
in config.yaml.

The output is a best-effort first draft. It is correct for monomer targets and
for dimers/oligomers when --chains-per-assembly is supplied. Always review the
file before running the collator.

Usage:
    python -m xchemalign.generate_assemblies -d <working_dir>
    python -m xchemalign.generate_assemblies -d <working_dir> --chains-per-assembly 2
"""

import argparse
import sys
from pathlib import Path

import gemmi

from xchemalign import utils
from xchemalign.utils import Constants

# Names for common assembly sizes
_ASM_SIZE_NAMES = {1: 'monomer', 2: 'dimer', 3: 'trimer', 4: 'tetramer'}


def _sanitise_space_group(sg: str) -> str:
    """'P 43 21 2' → 'P43212'"""
    return sg.replace(' ', '')


_MIN_POLYMER_LENGTH = 20


def _get_protein_chains(st: gemmi.Structure) -> list[str]:
    """Return auth chain names whose polymer span is at least _MIN_POLYMER_LENGTH residues."""
    st.setup_entities()
    return [c.name for c in st[0] if c.get_polymer().length() >= _MIN_POLYMER_LENGTH]


def _count_distinct_polymer_entities(st: gemmi.Structure) -> int:
    return sum(1 for e in st.entities if e.entity_type == gemmi.EntityType.Polymer)


def _assembly_type_name(chains_per: int) -> str:
    return _ASM_SIZE_NAMES.get(chains_per, f'{chains_per}mer')


def _read_cryst1_fast(pdb_path: Path) -> tuple[str, gemmi.UnitCell] | None:
    """
    Read only the CRYST1 record from a PDB file without parsing the whole structure.
    Returns (space_group_hm, UnitCell) or None if the record is absent or unreadable.
    """
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith('CRYST1'):
                    a = float(line[6:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                    alpha = float(line[33:40])
                    beta = float(line[40:47])
                    gamma = float(line[47:54])
                    sg = line[55:66].strip()
                    return sg, gemmi.UnitCell(a, b, c, alpha, beta, gamma)
                elif line.startswith(('ATOM', 'HETATM', 'MODEL')):
                    break  # CRYST1 always precedes coordinates
    except Exception:
        pass
    return None


def _cells_match(c1: gemmi.UnitCell, c2: gemmi.UnitCell, length_tol: float, angle_tol: float) -> bool:
    return (
        abs(c1.a - c2.a) <= length_tol
        and abs(c1.b - c2.b) <= length_tol
        and abs(c1.c - c2.c) <= length_tol
        and abs(c1.alpha - c2.alpha) <= angle_tol
        and abs(c1.beta - c2.beta) <= angle_tol
        and abs(c1.gamma - c2.gamma) <= angle_tol
    )


def _latest_pdb_for_crystal(xtal_dir: Path) -> Path | None:
    """Return the PDB from the highest-numbered Refine_NNNN subdir, preferring non-split files."""
    for refine_dir in reversed(sorted(xtal_dir.glob('Refine_*'))):
        pdbs = list(refine_dir.glob('*.pdb'))
        if pdbs:
            non_split = [p for p in pdbs if '.split.' not in p.name]
            return non_split[0] if non_split else pdbs[0]
    return None


def _collect_pdbs_for_scan(base_path: Path, inputs_cfg: list) -> list[tuple[Path, str]]:
    """
    Collect (pdb_path, source_type) pairs from all input directories.

    manual         → every *.pdb directly in the input dir
    model_building → one PDB per crystal subdir (latest Refine_NNNN)
    """
    results: list[tuple[Path, str]] = []
    for inp in inputs_cfg:
        inp_dir = inp.get(Constants.CONFIG_DIR)
        if not inp_dir:
            continue
        inp_path = base_path / Path(inp_dir)
        inp_type = inp.get(Constants.CONFIG_TYPE)

        if not inp_path.is_dir():
            continue

        if inp_type == Constants.CONFIG_TYPE_MANUAL:
            for pdb in sorted(inp_path.glob('*.pdb')):
                results.append((pdb, 'manual'))

        elif inp_type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            mb_root = inp_path / Constants.DEFAULT_MODEL_BUILDING_DIR
            if mb_root.is_dir():
                for xtal_dir in sorted(mb_root.iterdir()):
                    if xtal_dir.is_dir():
                        pdb = _latest_pdb_for_crystal(xtal_dir)
                        if pdb:
                            results.append((pdb, 'model_building'))
    return results


def _cluster_by_cell(
    entries: list[tuple[Path, str, str, gemmi.UnitCell]],
    length_tol: float,
    angle_tol: float,
    preferred_names: set[str],
) -> list[dict]:
    """
    Greedy clustering of PDB entries by space group + unit cell similarity.

    Each cluster dict contains:
        sg         – space group HM string (sanitised)
        cell       – UnitCell of the seed/representative
        rep        – Path to the representative PDB
        rep_type   – 'manual' | 'model_building'
        count      – number of structures in the cluster
    """
    clusters: list[dict] = []

    for pdb_path, src_type, sg, cell in entries:
        matched = next(
            (cl for cl in clusters if cl['sg'] == sg and _cells_match(cl['cell'], cell, length_tol, angle_tol)),
            None,
        )
        if matched is None:
            clusters.append(
                {
                    'sg': sg,
                    'cell': cell,
                    'rep': pdb_path,
                    'rep_type': src_type,
                    'count': 1,
                }
            )
        else:
            matched['count'] += 1
            # Upgrade representative: prefer ref_datasets names, then manual over model_building
            is_preferred = pdb_path.stem in preferred_names
            cur_preferred = matched['rep'].stem in preferred_names
            if is_preferred and not cur_preferred:
                matched['rep'] = pdb_path
                matched['rep_type'] = src_type
            elif not cur_preferred and src_type == 'manual' and matched['rep_type'] == 'model_building':
                matched['rep'] = pdb_path
                matched['rep_type'] = src_type

    return clusters


def _find_ref_pdb(ref_name: str, base_path: Path, inputs_cfg: list) -> Path | None:
    """
    Search all input directories for a PDB file whose stem matches ref_name.

    For manual inputs the canonical location is <input_dir>/<ref_name>.pdb.
    For model_building inputs we try the model_building subdirectory tree.
    Falls back to a recursive search under each input dir if the canonical
    location is not found.
    """
    candidates = []

    for inp in inputs_cfg:
        inp_dir = inp.get(Constants.CONFIG_DIR)
        if not inp_dir:
            continue
        inp_path = base_path / Path(inp_dir)
        inp_type = inp.get(Constants.CONFIG_TYPE)

        if inp_type == Constants.CONFIG_TYPE_MANUAL:
            p = inp_path / f'{ref_name}.pdb'
            if p.is_file():
                return p

        elif inp_type == Constants.CONFIG_TYPE_MODEL_BUILDING:
            mb_root = inp_path / Constants.DEFAULT_MODEL_BUILDING_DIR / ref_name
            if mb_root.is_dir():
                for f in sorted(mb_root.rglob('*.pdb')):
                    if '.split.' not in f.name:
                        return f
                for f in sorted(mb_root.rglob('*.pdb')):
                    return f

        if inp_path.is_dir():
            for f in sorted(inp_path.rglob(f'{ref_name}.pdb')):
                candidates.append(f)

    return candidates[0] if candidates else None


def _group_chains(chains: list[str], per: int) -> list[list[str]]:
    """Partition chains into groups of size `per`."""
    return [chains[i : i + per] for i in range(0, len(chains), per)]


def _locate_config(working_dir: Path) -> tuple[Path, Path]:
    """
    Return (config_path, output_dir) for the given working directory.
    Looks first in upload-current/, then in the directory itself.
    """
    upload_dir = working_dir / 'upload-current'
    if upload_dir.exists():
        return upload_dir / 'config.yaml', upload_dir
    return working_dir / 'config.yaml', working_dir


def _build_ref_info_from_refs(
    ref_datasets: list[str],
    base_path: Path,
    inputs_cfg: list,
    logger: utils.Logger,
) -> dict[str, dict]:
    """
    ref_datasets mode: build ref_info from the ref_datasets list in config.yaml.
    Returns ref_info dict keyed by dataset name.
    """
    ref_info: dict[str, dict] = {}

    for ref_name in ref_datasets:
        pdb_path = _find_ref_pdb(ref_name, base_path, inputs_cfg)
        if pdb_path is None:
            logger.warn(f'PDB not found for reference {ref_name!r} — skipping')
            continue
        try:
            st = gemmi.read_structure(str(pdb_path))
        except Exception as exc:
            logger.warn(f'Could not read {pdb_path}: {exc} — skipping {ref_name!r}')
            continue

        protein_chains = _get_protein_chains(st)
        if not protein_chains:
            logger.warn(f'No protein chains in {pdb_path} — skipping {ref_name!r}')
            continue

        sg = _sanitise_space_group(st.spacegroup_hm)
        ref_info[ref_name] = {
            'space_group': sg,
            'protein_chains': protein_chains,
            'distinct_entities': _count_distinct_polymer_entities(st),
            'pdb_path': pdb_path,
        }
        logger.info(ref_name)
        logger.info('  space group:    ', sg)
        logger.info('  protein chains: ', protein_chains)

    return ref_info


def _build_ref_info_from_scan(
    base_path: Path,
    inputs_cfg: list,
    length_tol: float,
    angle_tol: float,
    preferred_names: set[str],
    logger: utils.Logger,
) -> dict[str, dict]:
    """
    Scan mode: discover crystal forms by reading CRYST1 from all available PDB
    files and clustering by space group + unit cell similarity.
    Returns ref_info dict keyed by crystalform name (sg or sg_N).
    """
    all_pdbs = _collect_pdbs_for_scan(base_path, inputs_cfg)
    logger.info(f'Scanning {len(all_pdbs)} PDB file(s) across all input directories...')

    entries: list[tuple[Path, str, str, gemmi.UnitCell]] = []
    skipped = 0
    for pdb_path, src_type in all_pdbs:
        result = _read_cryst1_fast(pdb_path)
        if result is None:
            skipped += 1
            continue
        sg_raw, cell = result
        entries.append((pdb_path, src_type, _sanitise_space_group(sg_raw), cell))

    if skipped:
        logger.info(f'{skipped} PDB file(s) had no readable CRYST1 record and were skipped')

    if not entries:
        logger.error('No PDB files with CRYST1 records found')
        return {}

    clusters = _cluster_by_cell(entries, length_tol, angle_tol, preferred_names)

    sg_counts: dict[str, int] = {}
    for cl in clusters:
        sg_counts[cl['sg']] = sg_counts.get(cl['sg'], 0) + 1

    logger.info(f'Found {len(clusters)} crystal form(s) from {len(entries)} structure(s)')
    logger.info(f'  length tolerance: {length_tol} Å, angle tolerance: {angle_tol}°')

    ref_info: dict[str, dict] = {}
    sg_seen: dict[str, int] = {}

    for cluster in clusters:
        sg = cluster['sg']
        sg_seen[sg] = sg_seen.get(sg, 0) + 1
        cf_name = sg if sg_counts[sg] == 1 else f'{sg}_{sg_seen[sg]}'
        rep = cluster['rep']
        cell = cluster['cell']

        logger.info(f'{cf_name}: {cluster["count"]} structure(s)')
        logger.info(
            f'  cell: a={cell.a:.2f} b={cell.b:.2f} c={cell.c:.2f}'
            f'  α={cell.alpha:.2f} β={cell.beta:.2f} γ={cell.gamma:.2f}'
        )
        logger.info(f'  representative: {rep.stem}  ({cluster["rep_type"]})')

        try:
            st = gemmi.read_structure(str(rep))
        except Exception as exc:
            logger.warn(f'Could not read representative {rep}: {exc} — skipping {cf_name}')
            continue

        protein_chains = _get_protein_chains(st)
        if not protein_chains:
            logger.warn(f'No protein chains in {rep} — skipping {cf_name}')
            continue

        ref_info[cf_name] = {
            'space_group': sg,
            'protein_chains': protein_chains,
            'distinct_entities': _count_distinct_polymer_entities(st),
            'pdb_path': rep,
            'rep_name': rep.stem,
            'count': cluster['count'],
            'cell': cell,
        }

    return ref_info


def _emit_yaml(
    ref_info: dict[str, dict],
    ref_datasets: list[str],
    chains_per_assembly: int,
    assembly_name_override: str | None,
    scan_mode: bool,
    logger: utils.Logger,
) -> tuple[str, int]:
    """
    Build the YAML string from ref_info.
    Returns (yaml_content, updated_chains_per_assembly).
    """
    if not scan_mode:
        primary_key = next((r for r in ref_datasets if r in ref_info), next(iter(ref_info)))
    else:
        primary_key = next(iter(ref_info))

    primary_info = ref_info[primary_key]
    primary_chains = primary_info['protein_chains']

    if chains_per_assembly > len(primary_chains):
        logger.warn(
            f'--chains-per-assembly {chains_per_assembly} exceeds the number of protein chains '
            f'({len(primary_chains)}) in the primary reference — falling back to 1'
        )
        chains_per_assembly = 1

    assembly_chains = primary_chains[:chains_per_assembly]
    asm_name = assembly_name_override or _assembly_type_name(chains_per_assembly)
    biomol_str = ','.join(assembly_chains)
    primary_ref_label = primary_info.get('rep_name', primary_key) if scan_mode else primary_key

    # In ref_datasets mode: use the space group as the crystalform name when unique
    use_sg_as_cf_name = False
    if not scan_mode:
        sg_counts: dict[str, int] = {}
        for info in ref_info.values():
            sg = info['space_group']
            sg_counts[sg] = sg_counts.get(sg, 0) + 1
        use_sg_as_cf_name = all(v == 1 for v in sg_counts.values())
        if not use_sg_as_cf_name:
            logger.info(
                'Multiple references share the same space group; using reference names as '
                'crystalform names. Rename them in the output file if you prefer different names.'
            )

    lines: list[str] = [
        '# Draft assemblies.yaml generated by xchemalign.generate_assemblies',
        '# Review carefully before running the collator, especially:',
        '#   - assembly name and type (monomer / dimer / etc.)',
        '#   - biomol / chains mapping',
        '#   - crystalform names (rename from space-group strings if desired)',
        '#   - chain groupings for multi-copy crystal forms',
        '',
        'assemblies:',
        f'    {asm_name}:',
        f'        reference: {primary_ref_label}',
        f'        biomol: {biomol_str}',
        f'        chains: {biomol_str}',
        '',
        'crystalforms:',
    ]

    for cf_key, info in ref_info.items():
        protein_chains = info['protein_chains']
        groups = _group_chains(protein_chains, chains_per_assembly)

        remainder = len(protein_chains) % chains_per_assembly
        if remainder != 0:
            omitted = protein_chains[-remainder:]
            logger.warn(
                f'{cf_key!r} has {len(protein_chains)} protein chain(s), not a multiple of '
                f'--chains-per-assembly {chains_per_assembly}. Omitted: {omitted}.'
            )
            groups = [g for g in groups if len(g) == chains_per_assembly]

        if scan_mode:
            cf_name = cf_key
            ref_label = info['rep_name']
        else:
            cf_name = info['space_group'] if use_sg_as_cf_name else cf_key
            ref_label = cf_key

        parts = [f'space group: {info["space_group"]}']
        if scan_mode:
            parts.append(f'{info["count"]} structure(s)')
            cell = info['cell']
            parts.append(
                f'cell: {cell.a:.1f} {cell.b:.1f} {cell.c:.1f} / ' f'{cell.alpha:.1f} {cell.beta:.1f} {cell.gamma:.1f}'
            )
            parts.append(f'rep: {info["rep_name"]}')
        else:
            parts.append(f'ref: {cf_key}')
        comment = '  # ' + ', '.join(parts)

        lines.append(f'    {cf_name}:{comment}')
        lines.append(f'        reference: {ref_label}')
        lines.append('        assemblies:')
        for i, group in enumerate(groups, start=1):
            lines.append(f'            {i}:')
            lines.append(f'                assembly: {asm_name}')
            lines.append(f'                chains: {",".join(group)}')

    return '\n'.join(lines) + '\n', chains_per_assembly


def generate(
    working_dir: Path,
    chains_per_assembly: int = 1,
    assembly_name_override: str | None = None,
    force: bool = False,
    output_file: Path | None = None,
    scan: bool = False,
    length_tol: float = 5.0,
    angle_tol: float = 2.0,
    dry_run: bool = False,
    logger: utils.Logger | None = None,
) -> bool:
    if logger is None:
        logger = utils.Logger()

    config_file, out_dir = _locate_config(working_dir)

    if not config_file.is_file():
        logger.error(f'config.yaml not found at {config_file}')
        return False

    out_file = output_file if output_file is not None else out_dir / Constants.ASSEMBLIES_FILENAME
    if not dry_run and out_file.is_file() and not force:
        logger.error(f'{out_file} already exists. Use --force to overwrite.')
        return False

    config = utils.read_config_file(config_file)
    base_path = utils.find_path(config, Constants.CONFIG_BASE_DIR)
    if base_path is None:
        logger.error('base_dir not defined in config.yaml')
        return False

    ref_datasets: list[str] = utils.find_property(config, Constants.CONFIG_REF_DATASETS, default=[])
    inputs_cfg: list = utils.find_property(config, Constants.CONFIG_INPUTS, default=[])

    logger.info('Reading config from:', config_file)
    logger.info('base_dir:', base_path)

    if scan:
        logger.info('Mode: scan (discovering crystal forms from all PDB files)')
        ref_info = _build_ref_info_from_scan(
            base_path,
            inputs_cfg,
            length_tol,
            angle_tol,
            preferred_names=set(ref_datasets),
            logger=logger,
        )
    else:
        if not ref_datasets:
            logger.error('No ref_datasets defined in config.yaml. Use --scan to discover crystal forms automatically.')
            return False
        logger.info(f'Mode: ref_datasets ({len(ref_datasets)} reference(s))')
        ref_info = _build_ref_info_from_refs(ref_datasets, base_path, inputs_cfg, logger)

    if not ref_info:
        logger.error('No crystal forms could be determined.')
        return False

    content, _ = _emit_yaml(
        ref_info,
        ref_datasets,
        chains_per_assembly,
        assembly_name_override,
        scan_mode=scan,
        logger=logger,
    )

    if dry_run:
        logger.info('Dry run — output below (not written to disk)')
        print(content)
    else:
        with open(out_file, 'w') as fh:
            fh.write(content)
        logger.info('Wrote:', out_file)

    return True


def main():
    parser = argparse.ArgumentParser(description='Generate a draft assemblies.yaml from reference PDB files')
    parser.add_argument('-d', '--dir', default='.', help='Working directory (default: current dir)')
    parser.add_argument(
        '--chains-per-assembly',
        type=int,
        default=1,
        metavar='N',
        help='Number of protein chains that form one biological assembly unit (default: 1)',
    )
    parser.add_argument(
        '--assembly-name',
        metavar='NAME',
        help='Override the auto-generated assembly name (e.g. monomer, dimer)',
    )
    parser.add_argument(
        '--scan',
        action='store_true',
        help='Discover crystal forms by scanning all PDB files rather than using ref_datasets',
    )
    parser.add_argument(
        '--length-tolerance',
        type=float,
        default=5.0,
        metavar='Å',
        help='Unit cell length tolerance in Ångströms for crystal form clustering (default: 5.0)',
    )
    parser.add_argument(
        '--angle-tolerance',
        type=float,
        default=2.0,
        metavar='°',
        help='Unit cell angle tolerance in degrees for crystal form clustering (default: 2.0)',
    )
    parser.add_argument(
        '-o',
        '--output',
        metavar='FILE',
        help='Write output to this file instead of the default assemblies.yaml location',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print the generated YAML to stdout without writing any file',
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite an existing assemblies.yaml',
    )
    parser.add_argument(
        '--log-level',
        type=int,
        default=0,
        help='Logging level: 0=INFO, 1=WARN, 2=ERROR (default: 0)',
    )
    args = parser.parse_args()

    logger = utils.Logger(level=args.log_level)
    utils.LOG = logger

    working_dir = Path(args.dir).resolve()
    ok = generate(
        working_dir=working_dir,
        chains_per_assembly=args.chains_per_assembly,
        assembly_name_override=args.assembly_name,
        force=args.force,
        output_file=Path(args.output) if args.output else None,
        scan=args.scan,
        length_tol=args.length_tolerance,
        angle_tol=args.angle_tolerance,
        dry_run=args.dry_run,
        logger=logger,
    )
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
