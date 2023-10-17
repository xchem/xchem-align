# XChem Align User Guide

A guide for using the tools that generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

## Getting started (to use)

To run the XChem Align tools you can use a development environment
as described above or create a suitable user (run-time) environment, that does
not install the packages used for development: -

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip

    pip install .

Make sure you use Python 3.10. Earlier versions will not work, and later ones may work, but have not been tested.

## The tools

The following tools are being created, or planned. They typically are run in the order described.

To run these tools you need access to the Diamond file system, or to data from there that has been copied using
the `copier` tool.

### 1. Copier

This copies the necessary data from the Diamond file system to create an independent set of files that
can be worked with locally. The data for a single crystal is HUGE and it is not realistic to copy the
complete set of files, so `copier` just copies those that are required.

If you can run directly against the Diamond files system you do not need to use this tool.

This tool reads the SoakDB file (found at `processing/database/soakDBDataFile.sqlite`), reads the files that
are needed from the `mainTable` database table and copies those files to the output directory.
From there they can be tared up and  copied to your local system.

Usage:
```commandline
$ python -m xchemalign.copier --help
usage: copier.py [-h] -b BASE_DIR -i INPUT_DIR [-s SOAKDB_FILE] [-p [PANDDAS_FILES ...]] -o OUTPUT_DIR [-l LOG_FILE] [--log-level LOG_LEVEL]

copier

options:
  -h, --help            show this help message and exit
  -b BASE_DIR, --base-dir BASE_DIR
                        Base directory. If running against the Diamond file system use /
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Input directory (relative to base-dir) e.g. the dir with the data for your visit. e.g. dls/labxchem/data/2020/lb18145-153
  -s SOAKDB_FILE, --soakdb-file SOAKDB_FILE
                        Path to soakdb file relative to input-dir. Default is processing/database/soakDBDataFile.sqlite
  -p [PANDDAS_FILES ...], --panddas-files [PANDDAS_FILES ...]
                        Paths to CSV files with panddas data relative to input-dir
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory
  -l LOG_FILE, --log-file LOG_FILE
                        File to write logs to
  --log-level LOG_LEVEL
                        Logging level (0=INFO, 1=WARN, 2=ERROR)
```

Example:
```commandline
python -m xchemalign.copier -b / -i dls/labxchem/data/2020/lb18145-153 \
  -p processing/analysis/pandda_2_2020_04_25/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_04_30/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_05_10/analyses/pandda_inspect_events.csv \
  -o xce-outputs/lb18145-158
```

That example copies the required data for the visit found at `2020/lb18145-158` to the output directory
`xce-outputs/lb18145-158`. The full paths of the files are retained, and typically the required files can be located
anywhere under the `processing` directory in unpredictable locations.
Only the SoakDB database knows where the necessary files are actually located.

The files copied are:
* the soakdb sqlite database file
* the PDB and MTZ files for each crystal
* the CIF file for the ligand in the crystal
* the panddas event map data in `pandda_inspect_events.csv` files
* the *best* panddas event map (`.ccp4` file) identified from the `.csv` file

### 2. Collator

This prepares the input data that is needed and puts it in a standard location. It also creates a
`meta_collator.yaml` file listing those files and other necessary data. As such, this provides a consistent
staring point for the remaining steps. If your data is not coming from Diamond then you will need to
provide your own mechanism to generate this consistent input, and then you can utilise the following
steps.

IMPORTANT: before you run `collator` you need the following information:

1. The 'visits' that you want to include. e.g. at Diamond data for a visit will be located in a directory like this:
`/dls/labxchem/data/2021/lb29612-3`
2. The location of the soakDB database file. Usually this is at `processing/database/soakDBDataFile.sqlite` within the
visit directory, and this is used as the default.
3. The CSV files containing information about the Panddas event maps. This could be anywhere under `processing/analysis`
in the visit dir, and you, the user, will need to know which ones to use.
4. Any additional PDB files that you want to include (e.g. non-Diamond files). Each PDB file (and an optional .mtz file)
must be located in a directory with a name of the PDB code. See the `test-data/inputs_2/pdb` directory as an example.
5. Which datasets you need to specify as 'reference datasets'. These could be structures without ligands that you want
to use for  comparison.

`Collator` provides a mechanism for generating well-defined *releases* of data for your target. This is
done by creating directories in the output directory named `upload_1`, `upload_2` etc. To start a new release
just create the next `upload_*` directory and that will be used. The data previously generated in earlier
releases is not modified, but can be used in the subsequent releases.

Thus, to create your first release you must create your output directory and in it a directory named `upload_1`.

You must create a `config.yaml` file containing the configuration data. It can look like this:

```yaml
target_name: Mpro
base_dir: data/inputs
output_dir: data/outputs/Mpro
ref_datasets:
- 5rgs
- 8dz9
- 8e1y

inputs:

# 10 PDB files
- dir: dls/science/groups/i04-1/conor_dev/xchemalign/mpro_pdb
  type: manual

# 14 Mpro-J00* files
- dir: dls/labxchem/data/2021/lb29612-3
  type: model_building
  soakbd: processing/database/soakDBDataFile.sqlite
  panddas_event_files:
  - processing/analysis/panddas/pandda_20210716/analyses/pandda_inspect_events.csv
```

`base_dir` is required and is used to specify a path prefix should your input directories not reside at their
specified absolute file path (e.g. if you have used the `copier` tool). If running directly against the files at Diamond
then the base path should be `/`.

The `output_dir` must exist and contain the required `upload_*` directories (initially just a `upload_1` directory).

The `inputs` property is a list of *visits* (Diamond nomenclature) containing data for a number of crystals.
There can be one or more inputs. `inputs` is specified **relative** to `base_dir` (but if you need a different `base_dir`
then it can be specified here where it will override the one defined at the top level). The properties of `inputs` are the
data required to be processed, specified as paths **relative** to the `dir` path which is the base dir of the data
for that visit.

Once you have created the `config.yaml` file you can run the collator tool.

Usage:
```commandline
$ python -m xchemalign.collator --help
usage: collator.py [-h] [-c CONFIG_FILE] [-l LOG_FILE] [--log-level LOG_LEVEL] [-v]

collator

options:
  -h, --help            show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Configuration file
  -l LOG_FILE, --log-file LOG_FILE
                        File to write logs to
  --log-level LOG_LEVEL
                        Logging level
  -v, --validate        Only perform validation
```

We will use some test data that can be found in the `test-data` directory to illustrate how to run the tools.

You can use `collator` with the `--validate` option to just perform validation.

```commandline
$ python -m xchemalign.collator -c test-data/config_1.yaml --validate
INFO: initialising logging at level 0 at 2023-08-30 09:26:16.539489
INFO: collator:  Namespace(config_file='test-data/config_1.yaml', log_file=None, log_level=0, validate=True)
INFO: found 1 inputs
INFO: adding input dls/labxchem/data/2020/lb27995-1
INFO: version is 1
INFO: setting version dir to upload_1
INFO: Using version dir upload_1
INFO: validating paths
INFO: testing test-data/inputs_1/dls/labxchem/data/2020/lb27995-1
INFO: validating data
INFO: opening DB file: test-data/inputs_1/dls/labxchem/data/2020/lb27995-1/processing/database/soakDBDataFile.sqlite
INFO: adding crystal (model_building) Mpro-IBM0045
INFO: adding crystal (model_building) Mpro-IBM0058
INFO: adding crystal (model_building) Mpro-IBM0078
INFO: validator handled 3 rows from database, 3 were valid
```

Messages are reported at INFO, WARN and ERROR levels. You should look at WARN messages and assess whether they
need to be assessed, but you can proceed even if you see those. However, if you see ERROR messages then these
must be fixed before proceeding.

To run the full `collator` process run like this:

```commandline
$ python -m xchemalign.collator -c test-data/config_1.yaml --log-file collator.log
INFO: initialising logging at level 0 at 2023-08-30 10:01:08.723210
INFO: collator:  Namespace(config_file='test-data/config_1.yaml', log_file='collator.log', log_level=0, validate=False)
INFO: found 1 inputs
INFO: adding input dls/labxchem/data/2020/lb27995-1
INFO: version is 1
INFO: setting version dir to upload_1
INFO: Using version dir upload_1
INFO: validating paths
INFO: testing test-data/inputs_1/dls/labxchem/data/2020/lb27995-1
INFO: validating data
INFO: opening DB file: test-data/inputs_1/dls/labxchem/data/2020/lb27995-1/processing/database/soakDBDataFile.sqlite
INFO: adding crystal (model_building) Mpro-IBM0045
INFO: adding crystal (model_building) Mpro-IBM0058
INFO: adding crystal (model_building) Mpro-IBM0078
INFO: validator handled 3 rows from database, 3 were valid
INFO: running collator...
INFO: coping files ...
INFO: using cryst_dir of test-data/outputs/upload_1/crystallographic_files
INFO: creating cryst_dir test-data/outputs/upload_1/crystallographic_files
INFO: Mpro-IBM0045 has 3 files to copy
INFO: Mpro-IBM0058 has 3 files to copy
INFO: Mpro-IBM0078 has 3 files to copy
INFO: munging the history ...
INFO: munging current metadata
INFO: metadata 1 has 3 items
INFO: munging resulted in 3 total xtals, 3 are new or updated
INFO: writing metadata ...
INFO: copying config ...
INFO: run complete
```

We have specified a file to log the messages to (in addition to seeing them on the console).
The output directory will now contain data like this:

```
├── upload_1
  ├── config.yaml
  ├── crystallographic_files
  │   ├── Mpro-IBM0045
  │   │   ├── Mpro-IBM0045_1_A_1101.ccp4
  │   │   ├── Mpro-IBM0045.cif
  │   │   ├── Mpro-IBM0045.mtz
  │   │   └── Mpro-IBM0045.pdb
  │   ├── Mpro-IBM0058
  │   │   ├── Mpro-IBM0058_1_A_1101.ccp4
  │   │   ├── Mpro-IBM0058.cif
  │   │   ├── Mpro-IBM0058.mtz
  │   │   └── Mpro-IBM0058.pdb
  │   └── Mpro-IBM0078
  │       ├── Mpro-IBM0078_1_A_1101.ccp4
  │       ├── Mpro-IBM0078.cif
  │       ├── Mpro-IBM0078.mtz
  │       └── Mpro-IBM0078.pdb
  └── meta_collator.yaml

```

With real data there will be many more crystals than this.

The files for each crystal are present in a separate directory, and the directories and files have standardised
names. The `config.yaml` file is copied to the directory and a `meta_collator.yaml` file is created. This looks like
this:

```yaml
run_on: '2023-08-30 10:07:35.721039'
input_dirs:
- test-data/inputs_1/dls/labxchem/data/2020/lb27995-1
output_dir: test-data/outputs
version_number: 1
version_dir: upload_1
previous_version_dirs: []
crystals:
  Mpro-IBM0045:
    reference: true
    type: model_building
    last_updated: '2022-08-08 16:11:00'
    crystallographic_files:
      xtal_pdb:
        file: upload_1/crystallographic_files/Mpro-IBM0045/Mpro-IBM0045.pdb
        sha256: bf24e507365fc150732a5f5df09e771030630e44543c406a8c0dda4b7473290d
      xtal_mtz:
        file: upload_1/crystallographic_files/Mpro-IBM0045/Mpro-IBM0045.mtz
        sha256: b5e9d5e729b4284e4ff43bfbaccc4d2cf05fbe1744ee1603a4048ae32bb2b459
      ligand_cif:
        file: upload_1/crystallographic_files/Mpro-IBM0045/Mpro-IBM0045.cif
        sha256: 8bae84d77c55878f31ac3795720899da7bb52e14a77b36027da524502a7e6e11
        smiles: CN(Cc1ccc(Cl)c(Cl)c1)c1ccc(S(N)(=O)=O)cn1
      panddas_event_files:
      - file: upload_1/crystallographic_files/Mpro-IBM0045/Mpro-IBM0045_1_A_1101.ccp4
        sha256: abab5237ad7a44cf237df0a4fe9b8c6f30f18b545f7e8e3848424282a130a91f
        model: 1
        chain: A
        res: 1101
        index: 1
        bdc: 0.21
    status: new

  ...
```

For brevity the data for only one crystal, IBM0045, is shown.
The output lists the files and a sha256 digest for each one to allow to determine if it has changed between versions.

### 3. Aligner

This tool generates alignments of the individual protein chains based on sites that
can be determined automatically or manually. Its output will be aligned PDB, map and artefact files and are placed in
the current version directory (`upload_n`) in the `aligned_files` directory.

The expectation is that you will run *aligner* multiple times, tweaking the configuration (e.g. the site definitions)
until you are happy with the results. Then you would run the following tools to generate the release data for
Fragalysis. When you have new data to process you would create a new release (create a new `upload_*` directory)
and start again. Data in your previous `upload_*` directories would not be modified, but will still be used
(e.g. alignments) for your new release.

In addition to generating the aligned structures, this tool also extracts various components of the aligned molecules
to standard file formats that can be used for purposes outside Fragalysis (these files will also be downloadable from
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

To run the aligner you need to produce two other important bits of information, for the biological assemblies and the
crystal forms.

The biological assemblies (probably you only have one) is defined in `assemblies.yaml`, the default location being
in the output dir. For our example data it looks like this:

```yaml
dimer:
    reference: Mpro-IBM0045
    biomol: A,B
    chains: A, A(-x,y,-z)
```

The crystal forms are defined in `xtalforms.yaml`, the default location being in the output dir.
An example can be found [here](test-data/outputs/xtalforms.yaml) and a guide for preparing this can be found
[here](xtalforms_example.md). Again you may only have a single crystal form, or you may have more than one.

With those files created we can now run teh aligner.

Usage:
```commandline
$ python -m xchemalign.aligner -h
usage: aligner.py [-h] -d VERSION_DIR [-m METADATA_FILE] [-x XTALFORMS] [-a ASSEMBLIES] [-l LOG_FILE] [--log-level LOG_LEVEL] [--validate]

aligner

options:
  -h, --help            show this help message and exit
  -d VERSION_DIR, --version-dir VERSION_DIR
                        Path to version dir
  -m METADATA_FILE, --metadata_file METADATA_FILE
                        Metadata YAML file
  -x XTALFORMS, --xtalforms XTALFORMS
                        Crystal forms YAML file
  -a ASSEMBLIES, --assemblies ASSEMBLIES
                        Assemblies YAML file
  -l LOG_FILE, --log-file LOG_FILE
                        File to write logs to
  --log-level LOG_LEVEL
                        Logging level
  --validate            Only perform validation
```

With the `assemblies.yaml` and `xtalforms.yaml` present in their default locations you just need to specify the version
directory.

Example:
```commandline
$ python -m xchemalign.aligner -d test-data/outputs/upload_1 --log-file aligner.log
aligner:  Namespace(version_dir='test-data/outputs/upload_1', metadata_file='meta_collator.yaml', xtalforms=None, assemblies=None, log_file='aligner.log', log_level=0, validate=False)
INFO: initialising logging at level 0 at 2023-08-30 11:36:25.917039
INFO: Running aligner...
INFO: Performing alignments
INFO: First run! Not updating!

... lots of output omitted

INFO: extracting components
INFO: handling Mpro-IBM0045
INFO: extracting components Mpro-IBM0045 A 1101 Mpro-IBM0078+A+1101 test-data/outputs/upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101.pdb
INFO: handling Mpro-IBM0058
INFO: extracting components Mpro-IBM0058 A 1101 Mpro-IBM0078+A+1101 test-data/outputs/upload_1/aligned_files/Mpro-IBM0058/Mpro-IBM0058_A_1101_Mpro-IBM0078+A+1101.pdb
INFO: handling Mpro-IBM0078
INFO: extracting components Mpro-IBM0078 A 1101 Mpro-IBM0078+A+1101 test-data/outputs/upload_1/aligned_files/Mpro-IBM0078/Mpro-IBM0078_A_1101_Mpro-IBM0078+A+1101.pdb
```

The metadata generated by the `collator` tool is appended to with information on the aligned structures and written to
a `meta_aligner.yaml` file. An example crystal is:

```yaml
  Mpro-IBM0045:
    reference: true
    type: model_building
    last_updated: '2022-08-08 16:11:00'
    crystallographic_files:
      ... as before
    status: new
    assigned_xtalform: 1
    aligned_files:
      A:
        '1101':
          Mpro-IBM0078+A+1101: {structure: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101.pdb,
                                artefacts: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_artefacts.pdb,
                                event_map: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_event.ccp4,
                                x_map: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101.ccp4,
                                pdb_apo: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_apo.pdb,
                                pdb_apo_solv: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_apo-solv.pdb,
                                pdb_apo_desolv: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_apo-desolv.pdb,
                                ligand_mol: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_ligand.mol,
                                ligand_pdb: upload_1/aligned_files/Mpro-IBM0045/Mpro-IBM0045_A_1101_Mpro-IBM0078+A+1101_ligand.pdb,
                                ligand_smiles: CN(Cc1ccc(Cl)c(Cl)c1)c1ccc(S(N)(=O)=O)cn1}
```

Where the asymmetric unit comprises multiple monomers multiple sites get aligned (one for the site in each chain of the
dimer).

The `metadata_aligner.yaml` file also contains information about the macromolecular assemblies, the crystal forms that
were observed, the crystal form sites and the canonical sites, and the conformer sites.

#### Creating subsequent versions

When you get new or updated data then you can create a new version ready for upload.
the process is as follows:

1. create the new version dir e.g. `upload_2`. Do not modify anything in the old version dirs.
2. Update your `config.yaml` file. You can either edit the existing one (a copy is present in the old version dir as a
record of what was used for that version) or create a new one. Add the new data that you want to process e.g. add
an additional visit.
3. Update the `assemblies.yaml` and `xtalforms.yaml` files if there are new crystal forms or biological assemblies
observed (again these are copied into the version dir to keep a permanent record of what was used).
4. Run the collator tool again. It will automatically pick up the new version dir if you have named it correctly.
5. Run the Aligner tool again, specifying the new version dir to use.

### 4. Releaser

This tool is in preparation and will prepare a gzipped tar file for the release containing just the necessary files
(e.g. those that are new or updated) ready for loading into Fragalysis.

---

[conda]: https://docs.conda.io/en/latest/
[pre-commit]: https://pre-commit.com
