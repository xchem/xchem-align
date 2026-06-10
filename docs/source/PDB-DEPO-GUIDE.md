# PDB Deposition User Guide

The `pdb_deposition` tool generates CIF files suitable for deposition to the Protein Data Bank via PDBe.
It processes all crystals that have a `RefinementOutcome` status of `5 - Deposition ready` in SoakDB,
producing per-crystal structure CIF and structure-factor CIF files enriched with common metadata.

## Environment

You need a standard XCA environment in which to run this.
See the [User guide](USER-GUIDE.md) for details.

## Sequence and entity information

Protein sequence and entity information must be defined as described in the
[Sequence information section](USER-GUIDE.md#sequence-information) of the User Guide.

## Common metadata (`-m` / `--metadata-csv`)

You must supply a `deposition-metadata.csv` file containing information that applies to all crystals:
authors, institutions, grant details, crystal growth conditions, experimental method, etc.
Use an existing file as a starting point.

The file uses an mmcif-gen template format where each row defines one CIF field.
Columns are: category, item, headers, required, info, and one or more value columns.

### Token substitution

Several fields in the template support placeholder tokens that are expanded per-crystal at run time.
Tokens may appear in any field value, including `_struct.title` and `_struct_keywords.text`.

| Token | Replaced with |
|---|---|
| `$CompoundCode` | The compound code from SoakDB for this crystal |
| `$CrystalName` | The crystal name from SoakDB |
| `$ExternalCode2`…`$ExternalCode9` | Values from the compound-codes CSV (see below) |
| `$PoseID` | Comma-separated pose codes from the Fragalysis download CSV (see below) |

Any `$ExternalCode` token that cannot be substituted (no entry in the file, or column absent) is
removed from the output rather than left as a literal token.

### Example title template

```
_struct,title,,,,Crystal structure of Lysozyme in complex with $CompoundCode ($CrystalName) (Compound ID $ExternalCode3),,,,,,,,,,
,,,,,,,,,,,,,,,
```

### Example keywords template

```
_struct_keywords,text,,,,Diamond Light Source, I04-1, ASAP Discovery Consortium, $PoseID,,,,,,,,,,
```

## Compound codes CSV (`-c` / `--compound-codes-csv`)

An optional CSV file that supplies additional external identifiers for the `$ExternalCode` tokens.

- The **first column** must contain the crystal name.
- A **header row** must be present.
- The integer `n` in `$ExternalCoden` selects the **n-th column** (1-indexed), so `$ExternalCode3`
  takes the value from the third column.

Example:

```
CrystalName,CompoundCode,OpenBindId
Zika_NS5A-x0351,Z274555794,OB-00000142
Zika_NS5A-x0501,Z52214433,OB-00000150
Zika_NS5A-x0354,Z275165822,OB-00000163
```

## Fragalysis CSV (`-f` / `--fragalysis-csv`)

An optional CSV file downloaded from Fragalysis that provides pose codes for the `$PoseID` token.
This is the `metadata.csv` from a Fragalysis dataset download.

- **Column 1** — pose code (e.g. `A71EV2A-x0836a`)
- **Column 3** — crystal / experiment code that matches the SoakDB crystal name (e.g. `A71EV2A-x0836`)

When a crystal has more than one pose, all codes are written comma-separated:

```
A71EV2A-x0836a, A71EV2A-x0836b
```

If a crystal has no entry in this file, `$PoseID` is replaced with an empty string.

## Running pdb_deposition

```
python -m pdbdepo.pdb_deposition \
    -w <path-to-collator-output>/upload-current/upload_1 \
    -m <path-to>/deposition-metadata.csv \
    -c <path-to>/compound_codes.csv \
    -f <path-to>/fragalysis-download/metadata.csv
```

### All arguments

| Argument | Required | Description |
|---|---|---|
| `-w` / `--collator-dir` | Yes | Path to the collator output directory (the `upload_N` dir containing `config.yaml` and `meta_collator.yaml`) |
| `-m` / `--metadata-csv` | Yes | Path to the common metadata CSV used by mmcif-gen |
| `-c` / `--compound-codes-csv` | No | Path to the compound codes CSV supplying `$ExternalCode` values |
| `-f` / `--fragalysis-csv` | No | Path to the Fragalysis download `metadata.csv` supplying `$PoseID` values |
| `-d` / `--debug` | No | Copy source input files into the output directory alongside the generated CIF files |
| `--log-level` | No | Logging verbosity: `0` = INFO (default), `1` = WARN, `2` = ERROR |

## Output

All output is written to a `pdb_depo_files/` directory inside the collator output dir:

1. `pdb_depo.log` — log file for the run
2. `ligands.tab` — SMILES and InChI strings for all processed ligands
3. One subdirectory per crystal, each containing:
   - `<crystal>_struc.cif` — the 3D model CIF, ready for PDBe deposition
   - `<crystal>_sf.cif` — merged structure factors (from MTZ latest, MTZ free, and any PanDDA event maps)

When `-d` / `--debug` is specified, the intermediate source files used to build those outputs are
also copied into each crystal's subdirectory.
