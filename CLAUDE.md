# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project does

XChem Align processes raw XChem crystallography data and produces output suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/). The pipeline runs in three sequential steps: **collate → align → upload**.

## Commands

**Install** (Python 3.10 or 3.11 required):
```bash
pip install .
pip install .[dev]   # includes black, pylint, mypy, pre-commit
```

**Run the pipeline** (working directory must contain a `upload-current` symlink to the current version dir):
```bash
python -m xchemalign.collator -d <working_dir>
python -m xchemalign.aligner -d <working_dir>
python -m xchemalign.copier -h
python -m xchemalign.uploader -h
```

**Tests:**
```bash
pytest                              # all tests
pytest tests/test_integration.py   # full collate+align integration test
pytest tests/test_utils.py -k test_name  # single test
```
Integration tests write output to `test-data/outputs/` and depend on session-scoped fixtures that run in order (see `pytest-order`). The env var `XCA_GIT_REPO` should point to the repo root.

**Linting / formatting:**
```bash
black src tests                     # format (line-length 119)
pylint src/xchemalign               # lint (LNA package excluded by default)
mypy src tests                      # type check
pre-commit run --all-files          # runs black + pylint
```

## Architecture

### Package layout

```
src/
  xchemalign/               # main pipeline package
  ligand_neighbourhood_alignment/  # LNA: alignment algorithm (integrated sub-package)
    dep/                    # deprecated predecessor code, kept for reference
  pdbdepo/                  # PDB deposition tooling (separate concern)
```

### Pipeline steps

| Module | Step | Reads | Writes |
|---|---|---|---|
| `collator.py` | 1 – Collate | SoakDB sqlite, raw model-building dirs, `config.yaml`, `assemblies.yaml` | versioned `upload_N/` dirs with `meta_collator.yaml` |
| `aligner.py` | 2 – Align | collator output + assemblies | aligned structures, `meta_aligner.yaml` |
| `copier.py` | 1 alt | manual input dirs (no SoakDB) | same versioned upload dirs |
| `uploader.py` | 3 – Upload | aligner output | Fragalysis API or S3 |

### Working directory layout

```
<working_dir>/
  upload-current -> upload-vN/   # symlink, required at startup
  upload-vN/
    upload_1/                    # first version's data
    upload_2/                    # incremental updates
    config.yaml                  # user-supplied pipeline config
    assemblies.yaml              # crystallographic assembly definitions
```

### Key design patterns

- **`Constants` class** (`src/xchemalign/utils.py`): all config YAML keys, metadata field names, and directory name conventions are defined here. Always use `Constants.*` rather than raw strings.
- **`dt.py` (LNA)**: Pydantic v2 models for all crystallographic entities — `DatasetID`, `LigandID`, `XtalForm`, `Assembly`, `CanonicalSite`, etc. These are the primary data-transfer objects between pipeline stages.
- **`decoder/decoder.py`**: validates `config.yaml` and `assemblies.yaml` against JSON schemas before the pipeline runs. Validation errors are accumulated and reported rather than thrown immediately.
- **`dbreader.py`**: reads XChem SoakDB (SQLite), returns pandas DataFrames. Used only by `collator.py` and `copier.py` for the `model_building` input type; not used for `manual` input type.
- **Two input modes** in config: `type: model_building` (reads SoakDB) vs `type: manual` (explicit file lists).

### The LNA sub-package

`ligand_neighbourhood_alignment` (LNA) was originally a separate project and retains its own style (loguru for logging, pydantic models, `networkx` graphs). The `dep/` subdirectory contains the old implementation and is not used by the current pipeline. When modifying alignment logic, work in the top-level LNA modules (`alignment_core.py`, `update.py`, `structure_alignment.py`, etc.), not in `dep/`.

### gemmi workaround

There is an active workaround in `pdb_deposition.py` and `rename_chains_2a.py` for a gemmi upstream bug (tracked in project memory). See `memory/project_gemmi_workaround.md` for details on what to revert when the bug is fixed.

## Config file format

The user-supplied `config.yaml` (validated by `decoder/decoder.py`) specifies `target_name`, `base_dir`, `output_dir`, and a list of `inputs`. See `config-example.yaml` and the JSON schemas in `src/xchemalign/decoder/` for the full schema. The `assemblies.yaml` file describes biological assemblies and crystal forms; see `assemblies.yaml` in `test-data/` for examples.
