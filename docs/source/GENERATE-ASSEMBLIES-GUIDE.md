# Generating assemblies.yaml automatically

The `generate_assemblies` tool produces a draft `assemblies.yaml` by inspecting the
crystallographic data directly, removing the need to write the file by hand from scratch.
The output is always a **first draft** — review it before running the collator.

---

## Background

The `assemblies.yaml` file (described in the
[User Guide](USER-GUIDE.md#23-the-assemblies-yaml)) requires knowledge of:

* which space groups and crystal forms are present in the data
* which protein chains form the biological assembly
* how many copies of that assembly sit in the asymmetric unit of each crystal form

Filling this in manually is error-prone, especially for targets with many crystal forms.
`generate_assemblies` automates as much of this as possible from the PDB files themselves.

---

## Quick start

Run from your working directory (the one containing the `upload-current` symlink):

```commandline
python -m xchemalign.generate_assemblies
```

Or specify a different directory:

```commandline
python -m xchemalign.generate_assemblies -d /path/to/working/dir
```

The file is written to `upload-current/assemblies.yaml` by default.

---

## Two operating modes

### Mode 1: ref_datasets (default)

Uses the `ref_datasets` list already in `config.yaml` as the set of reference structures.
One crystal form is created per entry.

```commandline
python -m xchemalign.generate_assemblies
```

This is the fastest mode and is appropriate when `ref_datasets` has already been populated.

### Mode 2: scan

Scans **all** PDB files across every input directory, reads only their `CRYST1` record
(space group and unit cell dimensions), and clusters them into crystal forms automatically.
No prior knowledge of which structures are references is needed.

```commandline
python -m xchemalign.generate_assemblies --scan
```

`ref_datasets` entries are still used as preferred representatives when they fall inside a
discovered cluster (see [Representative selection](#representative-selection) below).

---

## Specifying the biological assembly

The tool cannot determine the biological assembly size from the PDB file alone —
two protein chains in the asymmetric unit could be two independent monomers or one dimer.
Use `--chains-per-assembly` to supply this information:

| Assembly type | Flag |
|---|---|
| Monomer (default) | *(omit the flag)* |
| Homodimer / heterodimer | `--chains-per-assembly 2` |
| Trimer | `--chains-per-assembly 3` |
| Tetramer | `--chains-per-assembly 4` |

The assembly name written into the file (`monomer`, `dimer`, `trimer`, `tetramer`) is
derived from this count. Override it with `--assembly-name`:

```commandline
python -m xchemalign.generate_assemblies --chains-per-assembly 2 --assembly-name dimer
```

---

## Scan mode: crystal form clustering

In `--scan` mode the tool groups structures by **space group** and **unit cell similarity**.
Two structures are placed in the same crystal form if they share the same space group and
their unit cell parameters are within the specified tolerances:

| Parameter | Flag | Default |
|---|---|---|
| Length tolerance (Å) | `--length-tolerance` | `5.0` |
| Angle tolerance (°) | `--angle-tolerance` | `2.0` |

Tighten the tolerances to split forms that differ subtly; widen them to absorb outlier
unit cells into the dominant form.

Example — tighter length tolerance:

```commandline
python -m xchemalign.generate_assemblies --scan --length-tolerance 2.0
```

### Crystal form naming

* If all discovered forms have a **unique space group**, the space group string is used as
  the crystalform name (e.g. `P43212`, `P1211`).
* If **multiple forms share a space group** (common when a target crystallises in several
  different packings of the same symmetry), forms are named with a counter suffix:
  `C121_1`, `C121_2`, etc.

These names can be freely renamed in the output file.

### Representative selection

For each cluster, the representative structure (used as the `reference` in the YAML) is
chosen by the following priority:

1. A structure whose dataset name appears in `ref_datasets` — these are the curated
   references the user has already identified.
2. A structure from a `manual` input directory — manually placed reference files are
   preferred over model-building outputs.
3. The first structure encountered in iteration order (alphabetical within each input
   directory, in the order directories appear in `config.yaml`).

The representative is used both as the `reference` for the crystal form and to determine
the protein chains and assembly type. If the automatically chosen representative is a
poor-quality structure, either add it to `ref_datasets` to override the selection or
edit the output file directly.

---

## CLI reference

```
python -m xchemalign.generate_assemblies [options]

  -d, --dir DIR              Working directory containing config.yaml (default: .)
  --chains-per-assembly N    Protein chains per biological assembly unit (default: 1)
  --assembly-name NAME       Override the auto-generated assembly name
  --scan                     Discover crystal forms by scanning all PDB files
  --length-tolerance Å       Unit cell length tolerance for clustering (default: 5.0)
  --angle-tolerance °        Unit cell angle tolerance for clustering (default: 2.0)
  --dry-run                  Print the generated YAML to stdout without writing a file
  -o, --output FILE          Write output to FILE instead of the default location
  --force                    Overwrite an existing assemblies.yaml
```

---

## Output location

By default the file is written to `<working_dir>/upload-current/assemblies.yaml`.
If `upload-current` does not exist (the working directory has not been initialised),
the file is written directly to `<working_dir>/assemblies.yaml`.

Use `-o` to override:

```commandline
python -m xchemalign.generate_assemblies --dry-run          # preview only
python -m xchemalign.generate_assemblies -o /tmp/draft.yaml  # write elsewhere for testing
```

---

## What is and is not auto-detected

| Information | Auto-detected? | Notes |
|---|---|---|
| Space group | Yes | Read from `CRYST1` |
| Unit cell dimensions | Yes | Read from `CRYST1` |
| Crystal form grouping | Yes (scan mode) | Clustering with tolerances |
| Protein chain names | Yes | Via `gemmi.setup_entities()` |
| Number of ASU copies | Yes | Chain count ÷ chains-per-assembly |
| Biological assembly size | **No** | Requires `--chains-per-assembly` |
| Assembly name (monomer/dimer/…) | Derived | From chain count; override with `--assembly-name` |
| Symmetry-operator chains | **No** | Must be added manually |
| Crystal form names | Derived | Space group string; rename as desired |

---

## Typical workflows

### Simple monomer target, refs already in config.yaml

```commandline
python -m xchemalign.generate_assemblies --dry-run
# review output, then write:
python -m xchemalign.generate_assemblies --force
```

### Dimer target, starting from scratch

```commandline
python -m xchemalign.generate_assemblies --scan --chains-per-assembly 2 --assembly-name dimer --dry-run
# review, then write:
python -m xchemalign.generate_assemblies --scan --chains-per-assembly 2 --assembly-name dimer --force
```

### Many crystal forms, exploring clustering thresholds

```commandline
# see how many forms are found with default tolerances
python -m xchemalign.generate_assemblies --scan --dry-run

# tighten if too few forms are being distinguished
python -m xchemalign.generate_assemblies --scan --length-tolerance 2.0 --dry-run

# widen if singleton outlier clusters are appearing
python -m xchemalign.generate_assemblies --scan --length-tolerance 8.0 --dry-run
```

---

## After generation

The output file contains comments flagging decisions that need human review.
Before running the collator:

1. **Check the assembly definition** — confirm `biomol` and `chains` are correct for
   your biological assembly.
2. **Rename crystalforms** if the auto-generated names (space group strings or `C121_1`
   style) are not meaningful for your target.
3. **Remove spurious crystal forms** — singleton clusters in scan mode are often outlier
   structures rather than genuine distinct forms.
4. **Validate** — the collator validates `assemblies.yaml` against its schema on startup;
   errors will be reported clearly.
