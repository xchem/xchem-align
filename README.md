# XChem Align

![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/xchem/xchem-align?include_prereleases)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![lint](https://github.com/xchem/xchem-align/actions/workflows/lint.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/lint.yaml)
[![test](https://github.com/xchem/xchem-align/actions/workflows/test.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/test.yaml)

Tools to generate data suitable for loading into Fragalysis.

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

## Available tools

The following tools are being created, or planned. They typically are run in the order described.

To run them you should create a conda environment using the `environment.yaml` file.

```commandline
conda env create -f environment.yaml
conda activate xchem-align
```

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
`metadata.yaml` file listing those files and other necessary data. As such, this provides a consistent
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
    │   ├── Mpro-x2910.cif
    │   ├── Mpro-x2910-event_1_1-BDC_0.36_map.ccp4
    │   ├── Mpro-x2910.mtz
    │   └── Mpro-x2910.pdb
    │   ├── ...
    └── metadata.yaml
```

Note that many crystals have been omitted for brevity. Typically, there will be 10's or even 100's of these.

The files for each crystal are present in a separate directory, and the directories and files have standardised
names. The `config.yaml` file is copied to the directory and a `metadata.yaml` file is created. This looks like
this:

```yaml
run_on: '2023-03-16 09:38:09.076971'
input_dirs:
- /dls/labxchem/data/2020/lb18145-153
- /dls/labxchem/data/2020/lb18145-158
output_dir: data/outputs/Mpro
crystals:
  Mpro-x0072:
    last_updated: 2021-08-12 10:39
    crystallographic_files:
      xtal_pdb:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0072/Mpro-x0072.pdb
        sha256: 18c2141bd45a211d37fbdb3169ce61c9cde28e190f693034986ee568b5f5fe7a
      xtal_mtz:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0072/Mpro-x0072.mtz
        sha256: 911b0db787fd38a1121e089691d9807c341d83b0c2185cacd87b875a4ee69f35
      ligand_cif:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0072/Mpro-x0072.cif
        sha256: 3f292dfd7781772fc7e6ce4be265121a1f9a4a2592fbcaa95da380ee71c51a2e
  Mpro-x0104:
    last_updated: 2021-09-21 23:20
    crystallographic_files:
      xtal_pdb:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0104/Mpro-x0104.pdb
        sha256: f011bf57137ffc33a299478b016100e218e7e611cc7f0e77aea126445d725d4a
      xtal_mtz:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0104/Mpro-x0104.mtz
        sha256: 28ebb1264f6d3e7a0a7d1c6a06da0970b64ed5df1c1d294728a880e545e2375e
      ligand_cif:
        file: data/outputs/Mpro/upload_1/crystallographic/Mpro-x0104/Mpro-x0104.cif
        sha256: d82f49029708aeaa5acb246030fedde86f5bc596f4cb38bd1cd8fc46203b8524
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

### 4. Extractor

This tool is in preparation and will extract out various representation of those aligned data, such as apo forms
of the protein, and the ligands in various format.

### 5. Releaser

This tool is in preparation and will prepare a zip file for the release containing just the necessary files (e.g.
those that are new or updated) ready for loading into Fragalysis.

## Building and contributing

Project dependencies are defined in the `pyproject.toml` file. From a
clean virtual environment you can install the run-time and development
dependencies like this:

    pip install --upgrade pip
    pip install .[dev]

The project also relies on CI that is run in GitLab using the actions defined
in the files you'll find in the `.githib/workflows` directory.

We also require the use of the Git [pre-commit] framework.
To get started review the pre-commit utility and then install
the pre-commit hooks: -

    pre-commit install

Now the project's rules will run on every commit, and you can check the
current health of your clone with: -

    pre-commit run --all-files

---

[pre-commit]: https://pre-commit.com
