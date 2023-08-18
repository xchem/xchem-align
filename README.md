# XChem Align

![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/xchem/xchem-align?include_prereleases)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![test](https://github.com/xchem/xchem-align/actions/workflows/test.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/test.yaml)
[![mypy](https://github.com/xchem/xchem-align/actions/workflows/mypy.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/mypy.yaml)

Tools to generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

## Prerequisites

* **Python 3.10** or later

## Getting started (to contribute)

Project dependencies are defined in the `pyproject.toml` file.

You will need to use Python 3.10 or later (a requirement of the `pyproject.toml` file).
If you prefer to use [conda] you can create a Python 3.10 environment using the
`environment.yaml` in this project, otherwise, if you have Python 3.10 or later,
you can create an environment using the built-in `venv` module: -

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip

Make sure you create the venv using Python 3.10 (or later).

From your clean virtual environment you can now install the run-time and development
dependencies like this: -

    pip install .[dev]

The project also relies on CI that is run in GitHub using the actions defined
in the files you'll find in the `.github/workflows` directory.

We also require the use of the Git [pre-commit] framework.
To get started review the pre-commit utility and then install
the pre-commit hooks with the command: -

    pre-commit install

Now the project's rules will run on every commit, and you can check the
current health of your clone with: -

    pre-commit run --all-files

## Getting started (to use)

To run the XChem Align tools you can use a development environment
as described above or create a suitable user (run-time) environment, that does
not install the packages used for development: -

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip

    pip install .

## The tools

The following tools are being created, or planned. They typically are run in the order described.

To run these tools you need access to the Diamond file system, or to data from there that has been copied using
the `copier` tool.

### 1. Copier

This copies the necessary data from the Diamond file system to create an independent set of files that
can be worked with locally. The data for a single crystal is HUGE and it is not realistic to copy the
complete set of files.
If you can run directly against the Diamond files system you do not necessarily need to use this tool.

This tool reads the SoakDB file (found at `processing/database/soakDBDataFile.sqlite`), reads the files that
are needed from `mainTable` and copies those files to the output directory. From there they can be tared up and
copied to your local system.

Usage:

```commandline
python -m xchemalign.copier -b / -i dls/labxchem/data/2020/lb18145-153 \
  -p processing/analysis/pandda_2_2020_04_25/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_04_30/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_05_10/analyses/pandda_inspect_events.csv \
  -o xce-outputs/lb18145-158
```

That example copies the required data for the visit found at `2020/lb18145-158` to the output directory
`xce-outputs/lb18145-158`. The full path names are retained, and typically the required files can be located
anywhere under the `processing` directory in unpredictable locations. Only the
SoakDB database knows where the necessary files are actually located.

The files copied are:
* the soakdb sqlite database file
* the PDB and MTZ files for each crystal
* the CIF file for the ligand in the crystal
* the panddas event map data in `pandda_inspect_events.csv` files
* the *best* panddas event map (`.ccp4` file) identified from the `.csv` file

### 2. Collator

This prepares the input data that is needed and puts it in a standard location. It also creates a
`metadata_collator.yaml` file listing those files and other necessary data. As such, this provides a consistent
staring point for the remaining steps. If your data is not coming from Diamond then you will need to
provide your own mechanism to generate this consistent input, and then you can utilise the following
steps.

Collator provides a mechanism for generating well-defined *releases* of data for your target. This is
done by creating directories in the output directory named `upload_1`, `upload_2` etc. To start a new release
just create the next `upload_*` directory and that will be used. The data previously generated in earlier
releases is not modified, but can be used in the subsequent steps.

Thus, to create your first release you must create your output directory and in it a directory named `upload_1`.

You must create a `config.yaml` file containing the configuration data. It can look like this:

```yaml
target_name: Mpro
base_dir: data/inputs
output_dir: data/outputs/Mpro
copy_dir: data/inputs/copied
inputs:
- dir: dls/labxchem/data/2020/lb18145-153
  soakbd: processing/database/soakDBDataFile.sqlite
  panddas_event_files:
  - processing/analysis/pandda_2_2020_04_25/analyses/pandda_inspect_events.csv
  - processing/analysis/pandda_2_2020_04_30/analyses/pandda_inspect_events.csv
  - processing/analysis/pandda_2_2020_05_10/analyses/pandda_inspect_events.csv

overrides:
  deprecations:
    Mpro-x0104: Cat did a wee on the crystal
    Mpro-x0165: Bad karma
```

`base_dir` is required and is used to specify a path prefix should your input directories not reside at their
specified absolute file path (e.g. if you have used the *copier* tool). If running directly against the files at Diamond
then the base path should be `/`.

The `output_dir` must exist and contain the required `upload_*` directories (initially just a `upload_1` directory).

The `inputs` property is a list of *visits* (Diamond nomenclature) containing data for a number of crystals.
There can be one or more inputs. `inputs` is specified **relative** to `base_dir`. The properties of `inputs` are the
data required to be processed, specified as paths **relative** to the `dir` path which is the base dir of the data
for that visit.

Then you can run the collator tool. You can use the `--validate` option to just perform validation.

```commandline
 python -m xchemalign.collator -c config.yaml --validate
```

Messages are reported at INFO, WARN and ERROR levels. You should look at WARN messages and assess whether they
need to be assessed, but you can proceed even if you see those. However, if you see ERROR messages then these
must be fixed before proceeding.

To run the full collator process run like this:

```commandline
 python -m xchemalign.collator -c config.yaml --log-file collator.log
```

We have specified a file to log the messages to (in addition to seeing them on the console).
The output directory will now contain data like this:

```
└── upload_1
    ├── config.yaml
    ├── crystallographic
    │   ├── Mpro-x0072
    │   │   ├── Mpro-x0072.cif
    │   │   ├── Mpro-x0072.mtz
    │   │   └── Mpro-x0072.pdb
    │   ├── Mpro-x2910
    │   │   ├── Mpro-x2910.cif
    │   │   ├── Mpro-x2910-event_1_1-BDC_0.36_map.ccp4
    │   │   ├── Mpro-x2910.mtz
    │   │   ├── Mpro-x2910.pdb
    │   ├── ...
    └── metadata.yaml
```

Note that many crystals have been omitted for brevity. Typically, there will be 10's or even 100's of these.

The files for each crystal are present in a separate directory, and the directories and files have standardised
names. The `config.yaml` file is copied to the directory and a `metadata_collator.yaml` file is created. This looks like
this:

```yaml
run_on: '2023-03-16 09:38:09.076971'
input_dirs:
- data/inputs/dls/science/groups/i04-1/conor_dev/xchemalign/mpro_pdb
- data/inputs/dls/labxchem/data/2020/lb27963-4
- data/inputs/dls/labxchem/data/2020/lb27995-1
- data/inputs/dls/labxchem/data/2021/lb29612-3
output_dir: data/outputs/Mpro
crystals:
  ...
  Mpro-J0013:
    type: model_building
    last_updated: '2022-12-09 09:56:00'
    crystallographic_files:
      xtal_pdb:
        file: upload_1/crystallographic_files/Mpro-J0013/Mpro-J0013.pdb
        sha256: 61e5d6d1e26ec35eb466df466ce93e74ac6e7770229858a83a8301db7291ea63
      xtal_mtz:
        file: upload_1/crystallographic_files/Mpro-J0013/Mpro-J0013.mtz
        sha256: d907e97b73d969b1cc701ab1f4ace3485552d61f54ca4b02d6b690cc7a95896c
      ligand_cif:
        file: upload_1/crystallographic_files/Mpro-J0013/Mpro-J0013.cif
        sha256: 109bcf8bfb4664d43b9e64ea685603194b79fb792c88c32b6c488ed57f84ef3c
        smiles: CC(C)N(Cc1ccccc1)C(=O)Cn1nnc2ccccc21
      panddas_event_files:
      - file: upload_1/crystallographic_files/Mpro-J0013/1_A_501.ccp4
        sha256: 4727c319d23df78ad85722eb94ce522fce20c238128a9c5ad74a25066dfdaca1
        model: 1
        chain: A
        res: 501
        index: 1
        bdc: 0.29
      - file: upload_1/crystallographic_files/Mpro-J0013/1_C_1.ccp4
        sha256: 7ad23dbac264096939ce814ff427614152b1a9c6c40174f6cd6901074a0b8fac
        model: 1
        chain: C
        res: 1
        index: 2
        bdc: 0.28
    status: new
  ...
  Mpro_Nterm-x0029:
    type: model_building
    last_updated: '2021-03-29 23:31:00'
    crystallographic_files:
      xtal_pdb:
        file: upload_1/crystallographic_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029.pdb
        sha256: 26bbb4e4602798c01eeaaeea197e135a464f9d80f734344c49f671162ba0f6c2
      xtal_mtz:
        file: upload_1/crystallographic_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029.mtz
        sha256: b21e7cb020ebcde067620ce5057de9c6405a24a00836d168c4ba867f411b8014
      ligand_cif:
        file: upload_1/crystallographic_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029.cif
        sha256: e529b25d0a9076b4d94d621b22539f5c1752d0b201f1dfa2c00ae7f8907c5979
        smiles: O=C1C[C@@H](Oc2cc(Cl)cc(OCCNC(=O)c3cc(=O)[nH]c4ccccc34)c2)N1
      panddas_event_files:
      - file: upload_1/crystallographic_files/Mpro_Nterm-x0029/1_A_501.ccp4
        sha256: 5bac7692c34f4ed3b38f1f78c4d9121d255aaba380cc53efc8e4120b7a6d9b2b
        model: 1
        chain: A
        res: 501
        index: 1
        bdc: 0.57
      - file: upload_1/crystallographic_files/Mpro_Nterm-x0029/1_B_401.ccp4
        sha256: f9aa4b331c2b4fd3713ac613f1d27128c4122df7a2ece6ecff77e0bcf321c812
        model: 1
        chain: B
        res: 401
        index: 3
        bdc: 0.55
    status: new
  ...
```

Again, data for many crystals have been omitted for brevity.
The output lists the files and a sha256 digest for each one to allow to determine if it has changed between versions.

### 3. Aligner

This tool is in preparation and will generate alignments of the individual protein chains based on sites that
can be determined automatically or manually. Its output will be aligned PDB files and are placed in the current
version directory (`upload_n`) in the `aligned` directory.

The expectation is that you will run *aligner* multiple times, tweaking the configuration (e.g. the site definitions)
until you are happy with the results. Then you would run the following tools to generate the release data for
Fragalysis. When you have new data to process you would create a new release (create a new `upload_*` directory)
and start again. Data in your previous `upload_*` directories would not be modified, but will still be used in
generating data (e.g. alignments) for your new release.

In addition to generating the aligned structures, this tool also extracts various components of the aligned molecules
to standard file formats that can be used for purposes outside Fragalysis (these wiles will also be downloadable from
Fragalysis). In summary these are:
- the aligned PDB file
- the aligned crystallographic artefacts
- the aligned Panddas event maps
- the aligned MTZ file
- the aligned PDB file without the ligand
- the aligned PDB file without the ligand or solvent (e.g. just the protein)
- only the solvent molecules
- the ligand in MDL molfile format
- the ligand in PDB format
- the ligand SMILES

The metadata generated by the `collator` tool is appended to with information on the aligned structures and written to
a `metadata_aligner.yaml` file. An example (using the same crystals as above) is:

```yaml
...
Mpro-J0013:
  type: model_building
  last_updated: '2022-12-09 09:56:00'
  crystallographic_files:
    ... (as above)
  aligned_files:
    A:
      501:
        0: {structure: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0.pdb, artefacts: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_artefacts.pdb,
            event_map: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_event.ccp4, x_map: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0.ccp4,
            pdb_apo: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_apo.pdb, pdb_apo_solv: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_apo-solv.pdb,
            pdb_apo_desolv: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_apo-desolv.pdb,
            ligand_mol: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_ligand.mol, ligand_pdb: aligned_files/Mpro-J0013/Mpro-J0013_A_501_0_ligand.pdb,
            ligand_smiles: CC(C)N(Cc1ccccc1)C(=O)Cn1nnc2ccccc21}
...
  Mpro_Nterm-x0029:
    type: model_building
    last_updated: '2021-03-29 23:31:00'
    crystallographic_files:
      ... (as above)
    aligned_files:
      A:
        501:
          0: {structure: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0.pdb,
            artefacts: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_artefacts.pdb,
            event_map: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_event.ccp4,
            x_map: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0.ccp4, pdb_apo: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_apo.pdb,
            pdb_apo_solv: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_apo-solv.pdb,
            pdb_apo_desolv: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_apo-desolv.pdb,
            ligand_mol: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_ligand.mol,
            ligand_pdb: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_A_501_0_ligand.pdb,
            ligand_smiles: 'O=C1C[C@H](Oc2cc(Cl)cc(OCCNC(=O)c3cc(=O)[nH]c4ccccc34)c2)N1'}
      B:
        401:
          0: {structure: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0.pdb,
            artefacts: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_artefacts.pdb,
            event_map: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_event.ccp4,
            x_map: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0.ccp4, pdb_apo: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_apo.pdb,
              pdb_apo_solv: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_apo-solv.pdb,
            pdb_apo_desolv: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_apo-desolv.pdb,
            ligand_mol: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_ligand.mol,
            ligand_pdb: aligned_files/Mpro_Nterm-x0029/Mpro_Nterm-x0029_B_401_0_ligand.pdb,
            ligand_smiles: 'O=C1C[C@H](Oc2cc(Cl)cc(OCCNC(=O)c3cc(=O)[nH]c4ccccc34)c2)N1'}
  ...
```

In that data you may notice that the crystal `Mpro_Nterm-x0029` where the asymmetric unit comprises a dimer where the
monomers are not quite identical results in two sites that get aligned (one for the active site in each chain of the
dimer).

The `metadata_aligner.yaml` file also contains information about the macromolecular assemblies, the crystalforms that
were observed, the crystalform sites and the canonical sites, and the conformer sites (TODO - document these).

### 4. Releaser

This tool is in preparation and will prepare a zip file for the release containing just the necessary files (e.g.
those that are new or updated) ready for loading into Fragalysis.

---

[conda]: https://docs.conda.io/en/latest/
[pre-commit]: https://pre-commit.com
