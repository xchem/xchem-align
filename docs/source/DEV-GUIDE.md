# XChemAlign Developer Guide

Tools to generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

## Prerequisites

* **Python 3.10** or **Python 3.11**. 3.12 is NOT yet supported.

## Getting started (to contribute)

Project dependencies are defined in the `pyproject.toml` file.

You will need to use Python 3.10 or 3.11 (a requirement of the `pyproject.toml` file).
Python 3.12 cannot currently be used.

If you prefer to use [conda] you can create a Python 3.10 environment using the
`environment.yaml` in this project, otherwise, if you have Python 3.10 or 3.11,
you can create an environment using the built-in `venv` module: -

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip

Make sure you create the venv using Python 3.10 or 3.11 (e.g. change the first command to `python3.11 -m venv venv`
if needed).

From your clean virtual environment you can now install the run-time and development
dependencies like this: -

    pip install -e .[dev]

The project also relies on CI that is run in GitHub using the actions defined
in the files you'll find in the `.github/workflows` directory.

We also require the use of the Git [pre-commit] framework.
To get started review the pre-commit utility and then install
the pre-commit hooks with the command: -

    pre-commit install

Now the project's rules will run on every commit, and you can check the
current health of your clone with: -

    pre-commit run --all-files

## Tools

The main tools are implemented as the following Python modules:

- Copier: xchemalign/copier.py
- Collator: xchemalign/collator.py
- Aligner: xchemalign/aligner.py

## Rollout

There is an environment at Diamond where users run the XChem align code on their data.
This can be found on the Diamond file system at `/dls/science/groups/i04-1/software/xchem-align`.
To roll out a new version of this:

1. Check that the repo is up to date on the `master` branch.
2. Test locally
3. Tag the XCA repo and push the tag: `git tag 1.2.3` and `git push origin 1.2.3` (using the appropriate tag number)
4. ssh to Diamond. Then ssh to wilson. Then run `srun --nodes=1 --ntasks=4 --partition=cs05r --pty bash`  and move into the `/dls/science/groups/i04-1/software/xchem-align` dir
5. `git pull` - update the repo
6. `git tag` - check you have the expected tag
7. `rm -rf env_xchem_align` - remove the old conda environment
8. `conda deactivate` - deactivate the current conda env (if necessary)
9. `conda env create -f environment.yaml -p env_xchem_align` - create the new conda environment

NOTE: the repo MUST be tagged before rolling out to users. Step 3 does this and is assumed to be done by the
developers. If you want to roll out a new environment and the repo is not tagged follow the instructions below

### Rollout if no tag is made yet

1. ssh to Diamond and move into the `/dls/science/groups/i04-1/software/xchem-align` dir
2. `git pull` - update the repo
3. `git tag` - check you have the expected tag
4. Tag the XCA repo and push the tag: `git tag 1.2.3` and `git push origin 1.2.3` (using the appropriate tag number)
5. `rm -rf env_xchem_align` - remove the old conda environment
6. `conda deactivate` - deactivate the current conda env (if necessary)
7. `conda env create -f environment.yaml -p env_xchem_align` - create the new conda environment

### Staging environment rollout

There is a staging environment at Diamond where users run a test version the XChem align code on their data.
This can be found on the Diamond file system at `/dls/science/groups/i04-1/software/xchem-align-staging`.
This code is not tagged and data should only be uploaded to the staging Fragalysis. Production Fragalysis will
refuse to load data from a non-tagged version of XCA.

The procedure for rolling out to the staging environment is identical to for the production site described above,
except that you do this in the `xchem-align-staging` directory. Ensure you checkout the required branch or tag before
building the conda environment.

### Troubleshooting the conda environment

#### `JSONDecodeError` when collecting package metadata

If `conda env create` fails while *Collecting package metadata (repodata.json)* with a traceback ending in
something like:

    conda.gateways.repodata.Response304ContentUnchanged
    ...
    json.decoder.JSONDecodeError: Unterminated string starting at: line 1 column 97806901

then conda's cached copy of the channel metadata is truncated. The server answered `304 Not Modified`
("your cached copy is current"), so conda fell back to reading the local cache and found it incomplete.
Retrying does not help, as the server keeps answering `304` and conda keeps reading the same broken file.

Clear the cache and try again:

    conda clean --index-cache
    conda env create -f environment.yaml -p env_xchem_align

#### Check your channels if it keeps happening

`environment.yaml` asks for `conda-forge`, but a bare channel name is resolved through your `~/.condarc`,
so it may not be fetching from where you expect. Check with:

    conda config --show-sources

and confirm the resolved URLs with:

    conda info | grep -A2 'channel URLs'

`conda-forge` should resolve to `https://conda.anaconda.org/conda-forge`. A `custom_multichannels`,
`custom_channels` or `channel_alias` entry redirecting it to a geographically distant mirror makes the
~100 MB metadata download slow and prone to being truncated, which is what causes the error above. Remove
the offending key, e.g.:

    conda config --remove-key custom_multichannels

leaving a `~/.condarc` of just:

    channel_priority: strict
    channels:
      - conda-forge

Note that the rollout is done inside an `srun` job (step 4 above), so a download interrupted by a job
limit will also leave a truncated cache behind.

Do not add Anaconda's `defaults` channel to work around channel problems. Anaconda's terms require a paid
licence for organisations of Diamond's size, and `conda-forge` alone provides everything `environment.yaml`
needs.


---

[conda]: https://docs.conda.io/en/latest/
[pre-commit]: https://pre-commit.com
