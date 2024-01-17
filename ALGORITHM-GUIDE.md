# _XChemAlign_ Algorithm Guide

## How Do You Know If Two Ligands Are Suitable For Merging?

The goal of the XChemAlign algorithm is to help users answer this question aims to address this question.

The design of follow-up compounds from crystallographic fragment screens is in general a very difficult issue, but it is further complicated by crystallography. In particular, identifying the relative geometries with which the compounds are likely to bind to the biological assembly, and hence what opportunities exist to merge or link them is non-trivial. Furthermore, understanding where design opportunities may be difficult to explore crystallographically, or when fragment binding may be an artefact of crystallography rather than biologically relevant is also non-trivial.

XChemAlign attempts to address these issues by identifying fragment binding sites, and then aligning structural data locally in order to provide the user with a clear picture of how fragments are likely bind relative to each other _in biology_, without crystallographic complications that make design opportunities difficult to recognise.



## Identification Of Fragment Binding Sites

Superposition of crystallographic structures around a subset of residues is well studied, as is the identification of crystallographic artefacts ins tructural data (when the biological assembly is known), which means that the primary challenge that XChemAlign deals with is the identification of which subsets are meaningful to perform these alignments around in order to generate those views of the data that best capture the available design opportunities in biology.

XChemAlign identifies three types fragment binding site:
1. **Conformer Sites**: The set of binders whose structural context can be aligned, _even indirectly_, with a good RMSD
2. **Canonical Sites**: The union of residues present in overlapping conformer sites
3. **Crystalform Sites**: The subset of conformer sites which share a crystalform

As their definitions suggest these three types of site form a heirarchy.

### Conformer Sites: Identifying Mutually Alignable Protein Binding Contexts

The first step in generating the biologically relevant binding site is to identify which compounds bind with a similar structural context.

In crystallography, the structural context is made up of two parts:
1. The Biological Assembly: The part of the context that resembles the protein in its native biological context
2. The Crystallographic Artefacts: The part of the context which contains atoms only in the compounds vicinity because of crystal packing

Once these two have been deconvoluted, which is relatively straightforward if the biological assembly is know, the biological binding context of a compound can be identified.

In order to generate the conformer site the next step is to identify which compounds have protein contexts which can be aligned to one another well. 

The key idea of a conformer site is that this alignment does not have to be direct: two compounds do not need to share any atoms in their local protein context if there exist other compounds whose contexts allow them to be aligned to one another indirectly. 

In this way two compounds which bind on the opposite sides of a pocket may still end up in the same conformer site if compounds which span the entire pocket, or at least the part of that pocket between these two compounds.


### Canonical Sites: Gathering Binders By Sequence

The key to identifying the biological context of fragments is by linking their binding mode to the sequence that generates the residues.

As however 




### Crystalform Sites: Partitioning Structure By Crystallographic Context

