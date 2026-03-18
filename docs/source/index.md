
# XChemAlign

XChemAlign (XCA) is a suite of command line tools for crystallographic fragment screening that identifies groupings of bound small-molecule fragments and aligns protein structures and ccp4 maps into a common reference frame.

## Getting started

XCA is typically run on fragment screening data that is collected by [XChem](https://https://www.diamond.ac.uk/Instruments/Mx/Fragment-Screening.html) at [Diamond Light Source](https://www.diamond.ac.uk/) and subsequently uploaded to [Fragalysis](https://fragalysis.readthedocs.io) for analysis. If you are new to XCA, we recommend starting with the [**XChemAlign user guide**](USER-GUIDE) to familiarise yourself with how it is run.

## Quick links

- [**fragalysis-app slack channel**](https://xchem-workspace.slack.com/archives/C02RCMA6S0Z)
- [**Fragalysis documentation**](https://fragalysis.readthedocs.io)
- [**DLS Fragalysis Production stack**](https://fragalysis.diamond.ac.uk)
- [**DLS Fragalysis Staging stack**](https://fragalysis.xchem.diamond.ac.uk)

---

# Documentation Pages

## XChemAlign User Guide
```{toctree}
:maxdepth: 1
USER-GUIDE.md
```
## XChemAlign Developer Guide
```{toctree}
:maxdepth: 1
DEV-GUIDE.md
```
## XChemAlign Algorithm Guide
```{toctree}
:maxdepth: 1
ALGORITHM-GUIDE.md
```

![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/xchem/xchem-align?include_prereleases)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![test](https://github.com/xchem/xchem-align/actions/workflows/test.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/test.yaml)
[![mypy](https://github.com/xchem/xchem-align/actions/workflows/mypy.yaml/badge.svg)](https://github.com/xchem/xchem-align/actions/workflows/mypy.yaml)

Tools to generate data suitable for loading into [Fragalysis](https://fragalysis.diamond.ac.uk/).

This supersedes [Fragalysis-API](https://github.com/xchem/fragalysis-api).

Contributors:
* Tim Dudgeon (Informatics Matters)
* Conor Wild (Diamond)
* Alan Christie (Informatics Matters)
* Kalev Takis (Informatics Matters)