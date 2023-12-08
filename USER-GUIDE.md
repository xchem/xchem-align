# _XChemAlign_ User Guide

_XChemAlign_ is a small suite of tools for preparing PDB models for loading into into [Fragalysis](https://fragalysis.diamond.ac.uk/).

* It formalises sites and packing artefacts across crystal forms and conformations, aligning models, maps and artefacts to common origins for each binding site.
* It handles model updates and multiple repeat experiments (e.g. to resolve stereochemistry).
* It allows for fast release cycles through supporting incremental updates.
* It assists efficient curation of auto-identified features, by running fast and on the minimal set of files in any given iteration.

## Overview

There a few steps involved.
1. [**Enable**](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#1-enabling-the-xchemalign-environment) the XChemAlign environment
1. [**Declare**](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#2-declaring-things) a few things about your data in two structured files in `yaml`[^1]
2. [**Collate**](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#3-collating-files) your files in a new (speicific) directory structure
3. [**Align**](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#4-aligning-everything) all binding sites to common origins
4. [**Release**](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#5-releasing-to-fragalysis) the data to Fragalysis
6. **Re-release** additional data by repeating (some or all) of steps 1-6.

If you won't run this at Diamond, you will first have to set up your environment and copy over files. See the [instructions below](https://github.com/mwinokan/xchem-align/blob/master/USER-GUIDE.md#non-diamond-instructions)

[^1]: "yet another markup language"

## 1. Enabling the XChemAlign environment

If you are uploading your data from diamond light source then this is as simple as running the commands:

```commandline
source /dls/science/groups/i04-1/conor_dev/xchem-align/act
conda activate /dls/science/groups/i04-1/conor_dev/xchem-align/env_xchem_align
```

## 2. Declaring things

In order to run XChemAlign you will need to create two files:
1. The config.yaml file
2. The crystalforms.yaml file

### 2.1. The Config Yaml

The config yaml defines what data to collect for collation. This will include raw crystalographic data, PanDDA data and ligand information.

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
inputs:  # The datasources to collate
  - dir: dls/labxchem/data/2020/lb27995-1  # The target directory. This will pull data from
                                            # 'dir/processing/analysis/modeL_building'. This is relative to 'base_dir'.
    type: model_building  # This will always be model_building unless you have datasets from the pdb you want to align
                          # which is an advanced topic not covered here.
    soakdb: processing/database/soakDBDataFile.sqlite  # The path to the soakdb database relative to 'dir'.
    exclude: [  # Datasets that are not to be processed with XChemAlign can be added to a list to exclude
      Mpro-IBM0057,
    ]
    panddas_event_files:  # The paths to the inspect tables of the PanDDAs used to model the bound state.
    - processing/analysis/panddas/analyses/pandda_inspect_events.csv  # Again these are relative to 'dir'.

```

### 2.2. The crystalforms Yamls

This file specifies both the biological *assemblies* and *crystalforms* relative to some reference PDBs. YAML has a strict formatting specification. Make sure to use spaces and not tabs for whitespace. The diagram below illustrates the format of the crystalforms file:

![crystalforms-yaml-example](https://github.com/xchem/xchem-align/assets/36866506/5c3ad74e-b1ff-4f44-8adb-3a76fbdc42b3)

The example file can be found [here](test-data/outputs/xtalforms.yaml). The `biomol` and `chains` directives specify the mapping between chains in the PDB file (`chains`) to chains in the assembly (`biomol`). I.e. in the example above the assembly "dimer-inhibited" is formed of three chains **A,B,C** which correspond to chains **C,E,A** in the **largecellpdb**.

## 3. Collating files

The first step is to collate your data. This process analyses your crystallographic data, PanDDA events, and ligand files and automatically determines the links between them.

```commandline
mkdir <path to your output_dir>/upload_1
python /dls/science/groups/i04-1/conor_dev/xchem-align/scripts/collate.py -c <your upload config file>
```

## 4. Aligning everything

The next step is performing local alignments of your ligand bound models and their associated crystallographic maps.

```commandline
python /dls/science/groups/i04-1/conor_dev/xchem-align/scripts/align.py -d <your upload directory> -x <your xtalforms file> -a <your assemblies file>
```

## 5. Releasing to Fragalysis

## Non-Diamond instructions

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

The next step (Collating) requires files  

> [!NOTE]
> Still reworking these instructions.  FvD 2023-10-17

