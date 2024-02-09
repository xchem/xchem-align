# _XChemAlign_ User Guide

_XChemAlign_ is a small suite of tools for preparing PDB models for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

* It formalises sites and packing artefacts across crystal forms and conformations, aligning models, maps and artefacts to common origins for each binding site.
* It handles model updates and multiple repeat experiments (e.g. to resolve stereochemistry).
* It allows for fast release cycles through supporting incremental updates.
* It assists efficient curation of auto-identified features, by running fast and on the minimal set of files in any given iteration.

## Overview

There are a few steps involved.
1. [**Enable**](#1-enabling-the-xchemalign-environment) the XChemAlign environment
2. [**Declare**](#2-declaring-things) a few things about your data in two structured files in `yaml`[^1]
3. [**Collate**](#3-collating-files) your files in a new (specific) directory structure
4. [**Align**](#4-aligning-everything) all binding sites to common origins
5. [**Release**](#5-upload-to-fragalysis) the data to Fragalysis
6. [**Re-release**](#6-creating-subsequent-versions) additional data

If you won't run this at Diamond, you will first have to set up your environment and copy over files. See the [instructions here](#non-diamond-instructions)

### Debugging common errors

* [Reporting version of the code](#reporting-version-of-the-code)
* [Missing PanDDA Event Files Warning When You Have Event Maps](#missing-pandda-event-files-warning-when-you-have-event-maps)
* [Missing PanDDAs Event Files Warning When No Event Maps Have Been Generated](#missing-panddas-event-files-warning-when-no-event-maps-have-been-generated)
* [Multiple Reference Structures](#missing-panddas-event-files-warning-when-no-event-maps-have-been-generated)
* [Adding PDB Structures To The Alignment](#adding-pdb-structures-to-the-alignment)

[^1]: "yet another markup language"

##

## 1. Enabling the XChemAlign environment

Uploading data from Diamond Light Source is as simple as running a few commands.

**EVERY TIME** you log in to run any XChemAlign tools: 

```commandline
source /dls/science/groups/i04-1/software/xchem-align/act
conda activate /dls/science/groups/i04-1/software/xchem-align/env_xchem_align
```

**THE VERY FIRST TIME** you must run this (the scripts above will fail and tell you).  _(It enables your linux account to read version info from the code's Git repository.)_:

```commandline
git config --global --add safe.directory /dls/science/groups/i04-1/software/xchem-align
```



## 2. Declaring things

In order to run XChemAlign you will need to create two files:
1. The config.yaml file
2. The assemblies.yaml file

### 2.1. The Config Yaml

The config yaml defines what data to collect for collation. This includes raw crystalographic data, PanDDA data and ligand information.

```yaml
# DO NOT USE TABS FOR THE WHITESPACE!

target_name: Mpro        # The name of your target.
                         # ??~~If adding to data already on Fragalysis, use that 'target' name~~??

base_dir: /some/path/to/test-data/inputs_1   # All _inputs_ are relative to this directory.
                                             # This is usually '/', certainly at Diamond

output_dir: /some/path/to/test-data/outputs  # Directory where ALL uploads folders go.
                                             # Must be an absolute path, NOT relative to base_dir
                                             
ref_datasets:        # List of datasets with reference conformations; these get aligned to every ligand binding site.
                     # You generally want at least one for each major class of conformation
  - Mpro-IBM0045        # Provide here the crystal ids as they appears in the model_building directory
  - Mpro-IBM0175

inputs:        # The datasources to collate

  - dir: dls/labxchem/data/2020/lb27995-1   # The visit directory; assumes processing/analysis/model_building is present
                                            # Path is _relative_ to 'base_dir'.

       type: model_building             # "model_building" means: XChem data

       soakdb: processing/database/soakDBDataFile.sqlite  # The path to the soakdb database relative to 'dir'.

       exclude: [Mpro-IBM0057, Mpro-IBM0108]   # Datasets to be ignored (e.g. if buggy)

       panddas_event_files:           # Tables written out by pandda_inspect; list here all pandda runs (XCA will figure it out from there)
         - processing/analysis/panddas/analyses/pandda_inspect_events.csv  # relative path, starting from 'dir'.

       panddas_missing_ok: [ Mpro-x0089, Mpro-x0211 ]    # Crystals for which XCA should ignore that event maps are missing.


  - dir: dls/labxchem/data/lb32633/lb32633-6/processing/analysis/additional_pdbs_forXCA
    type: manual       # each downloaded pdb file (cif!) goes in sub-directory.

```


### 2.2. The assemblies YAML

This file specifies both the biological *assemblies* and *crystalforms* relative to some reference PDBs.
YAML has a strict formatting specification. Make sure to use spaces and not tabs for whitespace.
The diagram below illustrates the format of the assemblies.yaml file:

![assemblies-yaml-example](https://github.com/xchem/xchem-align/assets/36866506/5c3ad74e-b1ff-4f44-8adb-3a76fbdc42b3)

An example file can be found [here](test-data/outputs/assemblies.yaml). The `biomol` and `chains` directives specify
the mapping between chains in the PDB file (`chains`) to chains in the assembly (`biomol`).
i.e. in the example above the assembly "dimer-inhibited" is formed of three chains **A,B,C** which correspond to chains
**C,E,A** in the **largecellpdb**.

### 2.3 Example configs

Here are some example configs that you can look at and run to hel get your head round how all this works.

#### 2.3.1 Minimal simple example.

This example illustrates only the minimal required configuration.
You will probably need to use additional configuration features, but this should help you understand the basics.

The example is contained in the file `example-simple.tgz`. Expand this file using:

```commandline
tar xvfz example-simple.tgz
```

Take a look at the two configuration files which are:

example-simple/work/config_1.yaml
example-simple/work/assemblies.yaml

You can run XChemAlign with this data using the instructions in the example-simple/README.txt file.

TODO - create a more complex example.

## 3. Collating files

The first step is to collate your data. This process analyses your crystallographic data, PanDDA events, and ligand files and automatically determines the links between them.

```commandline
mkdir <path to your output_dir>/upload_1
python /dls/science/groups/i04-1/software/xchem-align/scripts/collate.py -c <your upload config file>
```

## 4. Aligning everything

The next step is performing local alignments of your ligand bound models and their associated crystallographic maps.

```commandline
python /dls/science/groups/i04-1/software/xchem-align/scripts/align.py -d <your upload directory> -a <your assemblies file>
```
Note: the -a option is only needed if your assemblies file is not named `assemblies.yaml` and is not in `base_dir`.


## 5. Upload to Fragalysis

An automatic tool for Fragalysis upload has not yet been written.

To generate the gzipped tar file needed to manually upload the data move into your output dir and run this command
(updating it for the specific upload version and target name):

```
tar cvfz <target_name>.tgz upload_1
```

First to log in to Fragalysis and authenticate (for the staging server) and log in with your FedID:
https://fragalysis.xchem.diamond.ac.uk/viewer/react/landing

The gzipped tar file can then be uploaded to Fragalysis via (for the staging Fragalysis server):
[/viewer/upload_tset/ page](https://fragalysis.xchem.diamond.ac.uk/api/upload_target_experiments/)

## 6. Creating subsequent versions

When you have new or updated data you need to create a new version of the upload.
The data in `upload_1` must remain. Create a new directory named `upload_2` in the same place, update your
`config.yaml` to reflect the changes (or alternatively use names like `config_1.yaml` and `config_2.yaml`)
and then re-run *collator* and *aligner*. The commands will look like this:

```commandline
mkdir <path to your output_dir>/upload_2
python /dls/science/groups/i04-1/software/xchem-align/scripts/collate.py -c <your upload config file>
python /dls/science/groups/i04-1/software/xchem-align/scripts/align.py -d <your upload directory> -a <your assemblies file>
```

For a third version it's the same, just create and use a directory named `upload_3`.

Collator will automatically find and use the most recent version of the `upload_?` directory.
These must be named in sequence `upload_1`, `upload_2`, `upload_3` ...

When complete tar gzip the relevant `upload_?` dir and load into Fragalysis as before.
Fragalysis also only accepts uploads in the strict sequence described.

## 3. Debugging Errors

### Reporting version of the code.

If you successfully ran *collator* then the file `meta_collator.yaml` in your upload directory will contain a section like this:
```yaml
xca_git_info:
  origin_url: git@github.com:xchem/xchem-align.git
  branch: master
  sha: 4ba23bc73b89b6696d96baa80442122e975a0797
  tag: null
  dirty: false
```
This uniquely identifies the version of the code. Please report this if there is any doubt about the version being used.
If you can't run *collator* then the same info can be generated using: `python -m xchemalign.repo_info`

### Missing PanDDA Event Files Warning When You Have Event Maps

You may have missing datasets in your upload directory because the corresponding PanDDA event maps have not been found. The solution to this is to find the pandda_inspect_events.csv file corresponding to the PanDDA in which the structure was modelled and add it to the config.yaml.

```yaml
# config.yaml
# DO NOT USE TABS FOR THE WHITESPACE!
target_name: Mpro  # The name of your target. If you already have data on Fragalysis it should be the 'target' name that
                   # it appears under
base_dir: /some/path/to/test-data/inputs_1  # The directory that inputs (not output_dir!) are relative to. For users at
                                            # Diamond this should be set to '/'
...
inputs:  # The datasources to collate
  - dir: dls/labxchem/data/2020/lb27995-1  # The target directory. This will pull data from
                                            # 'dir/processing/analysis/modeL_building'. This is relative to 'base_dir'.
    type: model_building  # This will always be model_building unless you have datasets from the pdb you want to align
                          # which is an advanced topic not covered here.
...
    panddas_event_files:  # The paths to the inspect tables of the PanDDAs used to model the bound state.
    - processing/analysis/panddas/analyses/pandda_inspect_events.csv  # Again these are relative to 'dir'.
    - processing/analysis/panddas_2/analyses/pandda_inspect_events.csv # Structures now come from more than one PanDDA so add additional csv!


```

### Missing PanDDAs Event Files Warning When No Event Maps Have Been Generated

You may have refined structures that do not have PanDDA event maps, for example because you are working with follow up compounds for which the electron density is clear in conventional maps. By default XCA will warn you that this is a problem, however you can override this behaviour by setting the panddas_missing_ok key to the config.

```yaml
# config.yaml
# DO NOT USE TABS FOR THE WHITESPACE!
target_name: Mpro  # The name of your target. If you already have data on Fragalysis it should be the 'target' name that
                   # it appears under
base_dir: /some/path/to/test-data/inputs_1  # The directory that inputs (not output_dir!) are relative to. For users at
                                            # Diamond this should be set to '/'
...
    panddas_event_files:  # The paths to the inspect tables of the PanDDAs used to model the bound state.
    - processing/analysis/panddas/analyses/pandda_inspect_events.csv  # Again these are relative to 'dir'.
panddas_missing_ok: [  # List the dataset names that are in your model building directory that you want to export refined
  Mpro-x0089           # models for but that do not have corresponding PanDDA maps
]

```

### Multiple Reference Structures

You may have multiple reference datasets, for example because there are two major conformations that are present. This can be easily handled by adding multiple reference datasets in the config.

```yaml
# DO NOT USE TABS FOR THE WHITESPACE!
target_name: Mpro  # The name of your target. If you already have data on Fragalysis it should be the 'target' name that
                   # it appears under
base_dir: /some/path/to/test-data/inputs_1  # The directory that inputs (not output_dir!) are relative to. For users at
                                            # Diamond this should be set to '/'
output_dir: /some/path/to/test-data/outputs  # The directory that will contain all your upload folders. This path is
                                             # NOT relative to base_dir.
ref_datasets:  # A set of exemplar datasets that you want aligned to every ligand binding site. If you have multiple
              # major classes of conformations there should be at least one of each class.
  - Mpro-IBM0045  # There are given with the dataset folder name/crystal id as it appears in the
                  # model_building directory
  - Mpro-x0089 # A second reference dataset which will be aligned to every site discovered
...

```


### Adding PDB Structures To The Alignment

If you have structures from the PDB or some other, non-PanDDA, source to add, this can be managed by creating a "manual" minput directory, which is structures like a model_building directory (dataset names for folders which contain structures), and adding it to the config.

```yaml
# DO NOT USE TABS FOR THE WHITESPACE!
target_name: Mpro  # The name of your target. If you already have data on Fragalysis it should be the 'target' name that
                   # it appears under
base_dir: /some/path/to/test-data/inputs_1  # The directory that inputs (not output_dir!) are relative to. For users at
                                            # Diamond this should be set to '/'
output_dir: /some/path/to/test-data/outputs  # The directory that will contain all your upload folders. This path is
                                             # NOT relative to base_dir.
...
inputs:  # The datasources to collate
  - dir: dls/labxchem/data/2020/lb27995-1  # The target directory. This will pull data from
                                            # 'dir/processing/analysis/modeL_building'. This is relative to 'base_dir'.
    type: model_building  # This will always be model_building unless you have datasets from the pdb you want to align
                          # which is an advanced topic not covered here.
    soakdb: processing/database/soakDBDataFile.sqlite  # The path to the soakdb database relative to 'dir'.
    # Datasets that are not to be processed with XChemAlign can be added to a list to exclude
    exclude: [  
      Mpro-IBM0057,
    ]
    panddas_event_files:  # The paths to the inspect tables of the PanDDAs used to model the bound state.
    - processing/analysis/panddas/analyses/pandda_inspect_events.csv  # Again these are relative to 'dir'.
  - dir: path/to/some/dir  # Folder containing directories which contain PDB structures (and possibly corresponding MTZs)
    type: manual

```

# Non-Diamond instructions

1. **Set up** _(only once)_ your runtime environment _(easy)_
2. **Copy** relevant files from Diamond _(if not at Diamond)_

## 1. Setting up runtime environment _(only once)_

To run the XChemAlign tools you need to setup a Python environment.
This is described in more detail in the [Developer guide](DEV-GUIDE.md), but just to run the tools do this:

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip
    pip install .

Make sure you use Python 3.10. Earlier versions will not work, and later ones have not been tested.
Those steps just install what you need to run the tools, not to develop them.
You only need to set up this environment once.

## 2. Copying files from Diamond _(if not at Diamond)_

This copies the necessary data from the Diamond file system to create an independent set of files that
can be worked with locally. The data for a single crystal is HUGE and it is not realistic to copy the
complete set of files, so `copier` just copies those that are required.

If you can run directly against the Diamond files system you do not need to use this tool.

This tool reads the SoakDB file (found at `processing/database/soakDBDataFile.sqlite`), reads the files that
are needed from the `mainTable` database table and copies those files to the output directory.
From there they can be tarred up and copied to your local system.

Copier can be run in 2 modes:
1. SSH to Diamond, and copy the files there, (`copy` mode), then tar of zip up those copied files and copy the archive
   to your working environment.
2. Use SCP to copy the files from Diamond directly your working environment (`scp` mode).

Both modes require you to have SSH login credentials t Diamond (e.g. a FedID account and a SSH key).
`scp` mode is simpler to run, but a lot slower.

The settings can be provided in a `config.yaml` file or as commandline arguments, or as a bit of both.
When running in `copy` mode using commandline arguments is probably easier, when running in `scp` mode using a config
file is probably easier.

A `config.yaml` file will look like this:

```yaml
target_name: Mpro
output_dir: data/outputs/Mpro

scp:
  server: ssh.diamond.ac.uk
  username: gse84885
  key: /home/username/.ssh/id_rsa
  base_dir: /

inputs:
- dir: dls/labxchem/data/2020/lb18145-153
  type: model_building
  soakdb: processing/database/soakDBDataFile.sqlite
  panddas_event_files:
    - processing/analysis/pandda_2_2020_04_25/analyses/pandda_inspect_events.csv
    - processing/analysis/pandda_2_2020_04_30/analyses/pandda_inspect_events.csv
    - processing/analysis/pandda_2_2020_05_10/analyses/pandda_inspect_events.csv
```

Some notes about this file:

1. The same file is used to configure the collator tool, which you will run after copying. As the outputs of `copier`
   are the inputs of `collator` some of the names can be confusing.
2. The `scp` section contains your SSH credentials that are used when using `scp` mode.
3. The default for `scp.server` is `ssh.diamond.ac.uk` so you do not actually need this setting.
4. If you don't specify a SSH key then all your SSH keys (in the ~/ssh directory) are loaded.
5. The default for `scp.base_dir` is `/` so you do not actually need this setting.
6. There can be multiple values for the contents of `inputs`. If so all are copied, unless the `--input-dir` option is
   specified, in which case just that one input is copied.
7. When copying multiple inputs some of the commandline options become ambiguous so can't be used.
8. When not using a config file, only a single input can be copied. Run multiple times to copy multiple visits.
9. Commandline options take preference over settings in the config file.

Usage:
```commandline
$ python -m xchemalign.copier --help
usage: copier.py [-h] [-c CONFIG_FILE] [-b BASE_DIR] [-i INPUT_DIR] [-s SOAKDB_FILE] [-p [PANDDAS_FILES ...]] -m {copy,scp} [--scp-username SCP_USERNAME] [--scp-server SCP_SERVER] [--scp-key SCP_KEY] -o OUTPUT_DIR
                 [-l LOG_FILE] [--log-level LOG_LEVEL]

copier

options:
  -h, --help            show this help message and exit
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        Configuration file
  -b BASE_DIR, --base-dir BASE_DIR
                        Base directory. If running against the Diamond file system use /
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Input directory (relative to base-dir) e.g. the dir with the data for your visit. e.g. dls/labxchem/data/2020/lb18145-153
  -s SOAKDB_FILE, --soakdb-file SOAKDB_FILE
                        Path to soakdb file relative to input-dir. Default is processing/database/soakDBDataFile.sqlite
  -p [PANDDAS_FILES ...], --panddas-files [PANDDAS_FILES ...]
                        Paths to CSV files with panddas data relative to input-dir
  -m {copy,scp}, --mode {copy,scp}
                        Mode of file copying
  --scp-username SCP_USERNAME
                        SCP username
  --scp-server SCP_SERVER
                        SCP server
  --scp-key SCP_KEY     SSH key
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory
  -l LOG_FILE, --log-file LOG_FILE
                        File to write logs to
  --log-level LOG_LEVEL
                        Logging level (0=INFO, 1=WARN, 2=ERROR)
```

Example using `scp` mode with only commandline options:
```commandline
$ python -m xchemalign.copier -b / -i dls/labxchem/data/2020/lb18145-153 \
  -p processing/analysis/pandda_2_2020_04_25/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_04_30/analyses/pandda_inspect_events.csv \
  processing/analysis/pandda_2_2020_05_10/analyses/pandda_inspect_events.csv \
  -o xca-outputs -m scp --scp-user=username --scp-key ~/.ssh/id_rsa
```

That example copies the required data for the visit found at `/dls/labxchem/data/2020/lb18145-153` to the output directory
`xca-outputs`. The full paths of the files are retained, and typically the required files can be located
anywhere under the `processing` directory in unpredictable locations.
Only the SoakDB database knows where the necessary files are actually located.

The files copied are:
* the soakdb sqlite database file
* the PDB and MTZ files for each crystal
* the CIF file for the ligand in the crystal
* the panddas event map data in `pandda_inspect_events.csv` files
* the *best* panddas event map (`.ccp4` file) identified from the `.csv` file

Example using `scp` mode with config file options:
```commandline
$ python -m xchemalign.copier -c config.yaml -o xca-outputs -m scp
```

For using `copy` mode just change `--mode scp` to `--mode copy`.
The key difference is that `scp` mode is run remotely from your working environment whereas to use `copy` mode you must
first ssh to `ssh.diamond.ac.uk`, then run copier (having created a new working Python environment), then zip up the
copied files and copy them back to your working environment.
