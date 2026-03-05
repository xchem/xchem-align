# Ligand Neighbourhood Alignment

# Patch Notes:

## 2025 / 10 / 23

 * **Homolog Alignment**: LNA now supports the alignment of homologogous structures. Alignment is now based on matching atoms using aligned sequences rather than insertion id. This is much more general and flexible.
 * **Artefact Chains**: LNA now outputs alignments of artefact chains that have atoms in the neighbourhood of ligands.
 * **Ligand Altlocs**: LNA now supports ligands with altlocs, and will generate seperate alignments for each altloc of the ligand. In general LNA is now "aware" of altlocs and the relevant altloc will feature in automatically generated names and data structures.
 * **Full Biomol Alignments**: LNA now aligns the full assembly that a ligand is modelled as part of, not just the chains it is bound to.
 * **Paralleized Alignment**: Alignments of structures and maps are now processed in parallel using joblib
 * **Massive Refactoring**: LNA has been heavily refactored, and docstrings have been added to many functions.
 * **Additional Tests**: Additional test cases have been added to handle homolog alignment, artefact chains ligand altlocs.

# Package Structure
  * `alignment_hierarchy.py`: This file contains all the functionality for building the artificial "reference assembly" against which all final alignments are positioned.
  * `map_alignment.py`: This file contains the functionality for aligning xmaps/
  * `io.py`: This file contains the functionality for io not completely handled by external libraries
  * `matching.py`: Functionality for determining if two atoms should be aligned to one-another
  * `fs.py`: A model of the input and output file system structures for LNA
  * `constants.py`: Constants used throughout the program
  * `xtalform_assignment.py`: Heuristic for assigning crystalforms to datasets
  * `ligand_neighbourhoods.py`: Functions for determining the neighbourhoods of ligands
  * `nieghbourhood_graph.py`: Functions for alignment graph manipulation
  * `conformer_sites.py`: Heuristics to determine conformer sites
  * `xtalform_sites.py`: Functions to determine xtalform sites
  * `canonical_sites.py`: Heuristics to determine canonical sites
  * `forced_alignments.py`: Calculate forced alignments of datasets to reference assemblies
  * `structure_alignment.py`: Functionality for aligning structures based on sites
  * `dt.py`: Data structures and types
  * `alignment_landmarks.py`: Functionality for determining landmarks for alignment procedures
