# XChem Align Developer Guide

Tools to generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

## Prerequisites

* **Python 3.10** or later

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

There is an environment at Diamond where users user the XChem align code on their data.
This can be found on the Diamond file system at `/dls/science/groups/i04-1/software/xchem-align`.
To roll out a new version of this:

1. Check that the repos are up to date on the `master` (XCA) and `main` (LNA) branches.
2. Test locally
3. Tag the XCA repo and push the tag: `git tag 1.2.3` and `git push origin 1.2.3` (using the appropriate tag number)
4. ssh to Diamond and move into the `/dls/science/groups/i04-1/software/xchem-align` dir
5. `git pull` - update the repo
6. `git tag` - check you have the expected tag
7. `rm -rf env_xchem_align` - remove the old conda environment
8. `conda deactivate` - deactivate the current conda env (if necessary)
9. `conda env create -f environment.yaml -p env_xchem_align` - create the new conda environment

NOTE: the repo MUST be tagged before rolling out to users. Step 3 does this and is assumed to be done by the
developers. If you want to roll out a new environment and the repo is not tagged (step 3) then the commands in step 3
must instead be performed after step 5 and before step 6.

---

[conda]: https://docs.conda.io/en/latest/
[pre-commit]: https://pre-commit.com
