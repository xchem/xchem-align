# PDB deposition User Guide

The `pdb_deposition` tool generates data suitable for deposition to the Protein Databank through PDBe.
Entries that have a `RefinementOutcome` status in soakDB of `5 - Deposition ready` are processed, generating CIF files
suitable for deposition.

## Environment

You need a standard XCA environment in which to run this.
See the [User guide](USER-GUIDE.md) for details.

## Sequence and Entity information

These must be defined, as descriibed in the [Sequence information section](USER-GUIDE.md#sequence-information) of the
User Guide.

## Common metadata

You need to create a `deposition-metadata.csv` file containing common information that applies to all crystals that
are processed. This includes information such as authors, institutions etc.
Use an existing file as a starting point for this.

## MMCIF title field

The `deposition-metadata.csv` file has a `_struct,title` line that allows the `_struct` property to be defined.
The value is a template that is expanded with values know to the `pdb_depo` process. A example would be:
```
_struct,title,,,,Crystal structure of Lysozyme in complex with $CompoundCode ($CrystalName) (Compound ID $ExternalCode3),,,,,,,,,,
,,,,,,,,,,,,,,,
```

The 6th item is the template. The `$CompoundCode` and `$CrystalName` tokens are replace with `CompoundCode` and `CrystalName`
values from SoakDB for the specific crystal.

The `$ExternalCoden` token (when n is an integer between 1 and 9 inclusive) are values that can be defined in a CSV file
that looks like this:

```
CrystalName,CompoundCode,OpenBindId
Zika_NS5A-x0351,Z274555794,OB-00000142
Zika_NS5A-x0501,Z52214433,OB-00000150
Zika_NS5A-x0354,Z275165822,OB-00000163
```

- The first column MUST contain the crystal name.
- A header line MUST be present describing the columns.
- The integer value in the token is used to specify which column to use for the value e.g. in the example above,
`$ExternalCode3` refers to the 3rd column, the OpenBindId.

If any `$ExternalCoden` with n between 1 and 9 are not substituted (e.g. there was no record in the file for that crystal
or the column was not present) then the token is just removed from the title.

## Running pdb_deposition

Run as follows:

```
python -m pdbdepo.pdb_deposition -w <path-to-xca-collator-output>/upload-current/upload_1 -o <output-dir> -m <path-to>/deposition-metadata.csv -c <path-to>/compound_codes.csv
```
It is recommended to use a directory named pdb-depo for the -o option.

## Output

The output directory with contain:

1. `ligands.tab` containing the SMILES and InCHI for thee ligands
2. `pdb_dep.log` the log file (the location can be changed using the `-l` or `--log-file` option)
3. A directory for each crystal containing `*_struc.cif` (the 3D model), `*_sf.cif` (the merged structure factors from)
   the MTZ latest, MTZ free and any PanDDAs event maps.

If you specify the `-d` or `--debug` option the files that are used to generate those files are also copied to this directory.
