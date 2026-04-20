# _XChemAlign_ Algorithm Guide

## How Do You Know Which of A Set of Ligands Are Suitable For Merging?

The core goal of the XChemAlign algorithm is to help users answer this question by putting them all in the same reference frame.

The design of follow-up compounds from crystallographic fragment screens is in general a very difficult, but it is further complicated by crystallography. In particular, identifying the relative geometries with which the compounds are likely to bind to the biological assembly, and hence what opportunities exist to merge or link them can be non-trivial. Furthermore, understanding where design opportunities may be difficult to explore crystallographically, or when fragment binding may be an artefact of crystallography rather than biologically relevant is also non-trivial.

XChemAlign attempts to address these issues by identifying fragment binding sites, and then aligning structural data locally in order to provide the user with a clear picture of how fragments are likely bind relative to each other _in biology_, without crystallographic complications that make design opportunities difficult to recognise.


## Identification Of Fragment Binding Sites

Superposition of crystallographic structures around a subset of residues is well studied, as is the identification of crystallographic artefacts in structural data (when the biological assembly is known), which means that the primary challenge that XChemAlign deals with is the identification of which subsets are meaningful to perform these alignments around in order to generate those views of the data that best capture the available design opportunities in a biological context.

The primary object that XCA deals with is an **Observation**: an instance of a ligand binding observed in stuctural data. The same ligand binding with multiple modes in the same site, or binding at multiple sites in the same structure, count as multiple observations. 

XChemAlign also identifies three types fragment binding site that these Observations are partitioned among:
1. **Conformer Sites**: The set of binders whose structural context can be aligned, _even indirectly_, with a good RMSD i.e. the binding site is in the same conformation
2. **Canonical Sites**: The union of residues present in overlapping conformer sites i.e. the sequence ligands bind to
3. **Crystalform Sites**: The subset of conformer sites which share a crystalform i.e. the same binding site conformation with the same crystalographic artefacts

As their definitions suggest these three types of site form a heirarchy.

### Conformer Sites: Identifying Mutually Alignable Protein Binding Contexts

The first step in generating the biologically relevant binding site is to identify which compounds bind with a similar structural context.

In crystallography, the structural context is made up of two parts:
1. **The Biological Assembly**: The part of the context believed to resemble the protein in its native biological environment
2. **The Crystallographic Artefacts**: The part of the context which contains atoms believed to only be in the compounds vicinity because of crystal artefacts

Once these two have been deconvoluted, which is relatively straightforward if the biological assembly is known, the biological binding context of a compound can be identified.

In order to generate the conformer site the next step is to identify which compounds have protein contexts which can be aligned to one another well.

The key idea of a conformer site is that this alignment does not have to be direct: two compounds do not need to share any atoms in their local protein context if there exist other compounds whose contexts allow them to be aligned to one another indirectly.

In this way two compounds which bind on the opposite sides of a pocket may still end up in the same conformer site if other binders span the entire pocket, or at least the part of that pocket between these two compounds.


### Canonical Sites: Gathering Binders By Sequence

The key to identifying the biological context of fragments is by linking their binding mode to the sequence that generates the residues they bind to.

Once conformer sites have been calcuated, those with high overlap in sequence space can be identified and merged into the canonical sites. Canonical sites are defined by the union of these residues (a subsequence) and a reference structure which enables that sequence to be embedded into cartesian space for the purpose of generating final alignments. 


### Crystalform Sites: Partitioning Structure By Crystallographic Context

The final step for XCA is to identify groups of crystalographic artefacts. 

XCA is aware of the crystalforms present in the screen, and hence can partition the datasets between them. Then it is straightforward for each conformer site to also be partitioned according to which crystalforms are present among the Observations in it.   

Furthermore, the artefact environemnt of each Crystalform Site is tracked and made available in final alignments, so that it is easy to identify when which ligands might be binding because of crystal artefacts, or which crystalforms new designs might be possible to observe in.