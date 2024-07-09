# XChem Align Developer Guide

Tools to generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

## Prerequisites

* **Python 3.10** or **Python 3.11**. 3.12 is NOT yet supported.

## Getting started (to contribute)

Project dependencies are defined in the `pyproject.toml` file.

You will need to use Python 3.10 or 3.11 (a requirement of the `pyproject.toml` file).
If you prefer to use [conda] you can create a Python 3.10 environment using the
`environment.yaml` in this project, otherwise, if you have Python 3.10 or 3.11,
you can create an environment using the built-in `venv` module: -

    python -m venv venv
    source venv/bin/activate
    pip install --upgrade pip

Make sure you create the venv using Python 3.10 or 3.11.

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

---

[conda]: https://docs.conda.io/en/latest/
[pre-commit]: https://pre-commit.com
