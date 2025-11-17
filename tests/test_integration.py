from pathlib import Path

import pytest
import yaml

from rich import print as rprint
import gemmi

from xchemalign.collator import Collator
from xchemalign.aligner import Aligner
from xchemalign import utils
from xchemalign.utils import Constants

from ligand_neighbourhood_alignment import alignment_heirarchy


def test_collator_upload_1(
    constants,
    test_dir,
    uploads_dir,
    config_1_file,
    upload_1_dir,
    assemblies_file,
):
    c = Collator(test_dir)
    logger = c.logger

    meta = c.validate()

    if meta is None or len(logger.errors) > 0:
        print("There are errors, cannot continue")
        print(logger.errors)
        exit(1)
    else:
        c.run(meta)


@pytest.mark.order(after="test_collator_upload_1")
def test_aligner_upload_1(constants, assemblies_file):
    log = str(Path(constants.TEST_DIR) / "aligner.log")

    a = Aligner(constants.TEST_DIR, log_file=log, log_level=0)
    logger = a.logger
    utils.LOG = logger

    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        print(a.logger.errors)
        exit(1)
    else:
        meta = a.run()

    # Makre sure there are two observations for Mpro-IBM0078 A1101's different altlocs
    assert 'C' in meta[Constants.META_XTALS]["Mpro-IBM0078"][Constants.META_ALIGNED_FILES]['A']['1101']
    assert 'D' in meta[Constants.META_XTALS]["Mpro-IBM0078"][Constants.META_ALIGNED_FILES]['A']['1101']



@pytest.mark.order(after="test_aligner_upload_1")
def test_collator_upload_2(constants, config_2_file, upload_2_dir, uploads_dir, test_dir):
    c = Collator(test_dir)
    logger = c.logger

    meta = c.validate()

    print(meta)

    if meta is None or len(logger.errors) > 0:
        print("There are errors, cannot continue")
        print(logger.errors)
        exit(1)
    else:
        c.run(meta)
    assert (
        # len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {}))
        len(
            meta[Constants.META_XTALS]["Mpro-x0107_fake_P1"][Constants.META_XTAL_FILES].get(
                Constants.META_BINDING_EVENT, {}
            )
        )
        != 0
    )


@pytest.mark.order(after="test_collator_upload_2")
def test_aligner_upload_2(constants):
    log = str(Path(constants.TEST_DIR) / "aligner.log")

    a = Aligner(constants.TEST_DIR, log_file=log, log_level=0)
    logger = a.logger
    utils.LOG = logger

    num_errors, num_warnings = a.validate()

    if num_errors:
        print("There are errors, cannot continue")
        print(a.logger.errors)
        exit(1)
    else:
        a.run()

    # Make sure there are aligned files for Mpro-i0130
    aligned_file_dirs= [x.name for x in (Path(constants.UPLOAD_2_DIR) / "aligned_files").glob("*")]
    if "Mpro-i0130" not in aligned_file_dirs:
        raise Exception(f'Mpro-i0130 should be in {aligned_file_dirs}') 



# @pytest.mark.order(after="test_aligner_upload_2")
# def test_collator_upload_3(
#     constants,
#     test_dir,
#     upload_3_dir,
#     upload_3_data_dir,
#     config_3_file,
# ):
#     c = Collator(str(config_3_file), )
#
#     meta, num_errors, num_warnings = c.validate()
#
#     if meta is None or num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         c.run(meta)
#
#     with open(Path(upload_3_dir) / "meta_collator.yaml", 'r') as f:
#         new_meta = yaml.safe_load(f)
#     assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "crystallographic_files").glob("*")]
#
#
# @pytest.mark.order(after="test_collator_upload_3")
# def test_aligner_upload_3(constants, xtalforms_file, upload_3_dir):
#     a = Aligner(upload_3_dir, constants.METADATA_FILE, xtalforms_file, )
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         a.run()
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "aligned_files").glob("*")]


# what follows is from former LNA package


def test_derive_alignment_heirarchy(constants, lna_assemblies):
    # Run _derive_alignment_heirarchy
    reference_assemblies, chain_priority = alignment_heirarchy._derive_alignment_heirarchy(lna_assemblies, debug=True)

    # Print important info
    # rprint(reference_assemblies)
    ...


def test_chain_to_biochain(constants, lna_assemblies, xtalforms):
    rprint(xtalforms)
    for xtalform_name, xtalform in xtalforms.items():
        for xtalform_assembly_name, xtalform_assembly in xtalform.assemblies.items():
            for chain in xtalform_assembly.chains:
                biochain = alignment_heirarchy._chain_to_biochain(chain, xtalform, lna_assemblies)
                rprint(f"Chain to Biochain: {xtalform_name} : {xtalform_assembly_name} : {chain} -> {biochain}")


def test_structure_to_landmarks(pdb_paths):
    st = gemmi.read_structure(str(pdb_paths["Mpro-IBM0045"]))
    landmarks = alignment_heirarchy.structure_to_landmarks(st)
    rprint(landmarks)


def test_calculate_assembly_transform(pdb_paths, lna_assemblies):
    # Generate assembly structure from reference
    as1_ref = gemmi.read_structure(str(pdb_paths["Mpro-IBM0045"]))
    as1 = alignment_heirarchy._get_assembly_st(lna_assemblies["dimer"], as1_ref)

    as2_ref = gemmi.read_structure(str(pdb_paths["7ql8"]))
    as2 = alignment_heirarchy._get_assembly_st(lna_assemblies["monomer"], as2_ref)

    # Generate assembly landmarks from assembly structure
    as1_lm = alignment_heirarchy.structure_to_landmarks(as1)
    as2_lm = alignment_heirarchy.structure_to_landmarks(as2)

    # Generate alignment hierarchy
    hierarchy, chain_priority = alignment_heirarchy._derive_alignment_heirarchy(lna_assemblies)

    # Determine transform
    transform = alignment_heirarchy._calculate_assembly_transform(
        ref=as1_lm, mov=as2_lm, chain=hierarchy["monomer"][1], debug=True
    )
    rprint(transform)


def test_calculate_assembly_sequence(lna_assemblies):
    # Generate alignment hierarchy
    hierarchy, chain_priority = alignment_heirarchy._derive_alignment_heirarchy(lna_assemblies)

    #
    transform_sequence = alignment_heirarchy._calculate_assembly_sequence(hierarchy, "fake_tetramer")
    rprint(transform_sequence)


def test_calculate_assembly_transform_sequence(pdb_paths, lna_assemblies):
    # Get the assembly structures

    # Get the landmarks for the assembly references
    landmarks = {}
    for assembly_name, assembly in lna_assemblies.items():
        ref_st = gemmi.read_structure(str(pdb_paths[assembly.reference]))
        as_st = alignment_heirarchy._get_assembly_st(assembly, ref_st)
        landmarks[assembly_name] = alignment_heirarchy.structure_to_landmarks(as_st)

    # Get the hierarchy
    hierarchy, chain_priority = alignment_heirarchy._derive_alignment_heirarchy(lna_assemblies)

    # Calculate the full transform
    combined_transform = alignment_heirarchy._calculate_assembly_transform_sequence(
        hierarchy, "fake_tetramer", landmarks, debug=True
    )
    rprint(combined_transform)
    tr = alignment_heirarchy._transform_to_gemmi(combined_transform)
    mov = landmarks["fake_tetramer"][("C", ("GLN", "256"), "CA")]
    ref = landmarks["dimer"][("B", ("GLN", "256"), "CA")]
    aligned_mov = tr.apply(gemmi.Position(mov[0], mov[1], mov[2]))
    rprint(f"Reference position: {ref} <- aligned position: {(aligned_mov.x, aligned_mov.y, aligned_mov.z)}")


def test_get_structure_chain_to_assembly_transform(lna_assemblies, xtalforms, pdb_paths):
    # Get the landmarks for assembly references
    landmarks = {}
    for assembly_name, assembly in lna_assemblies.items():
        ref_st = gemmi.read_structure(str(pdb_paths[assembly.reference]))
        as_st = alignment_heirarchy._get_assembly_st(assembly, ref_st)
        landmarks[assembly_name] = alignment_heirarchy.structure_to_landmarks(as_st)

    # Select a test structure, chain and assembly
    test_dataset = "8e1y"
    st = gemmi.read_structure(str(pdb_paths[test_dataset]))
    xtalform = xtalforms["xtalform3"]
    chain = "B"

    #
    tr = alignment_heirarchy._get_structure_chain_to_assembly_transform(st, chain, xtalform, lna_assemblies, landmarks)
    rprint(tr)


# @pytest.mark.order(after="test_collator_upload_1")
# def test_aligner_upload_1(constants, assemblies_file, upload_1_dir):
#     a = Aligner(upload_1_dir, constants.METADATA_FILE, assemblies_file)
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         print(a.logger.errors)
#         exit(1)
#     else:
#         a.run()
#
#
# @pytest.mark.order(after="test_aligner_upload_1")
# def test_collator_upload_2(
#     constants,
#     test_dir,
#     upload_2_dir,
#     upload_2_data_dir,
#     config_2_file,
# ):
#     c = Collator(
#         str(config_2_file),
#     )
#     logger = c.logger
#
#     meta = c.validate()
#
#     if meta is None or len(logger.errors) > 0:
#         print("There are errors, cannot continue")
#         print(logger.errors)
#         exit(1)
#     else:
#         c.run(meta)
#     assert (
#         len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {}))
#         != 0
#     )
#
#
# @pytest.mark.order(after="test_collator_upload_2")
# def test_aligner_upload_2(constants, assemblies_file, upload_2_dir):
#     a = Aligner(upload_2_dir, constants.METADATA_FILE, assemblies_file)
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         print(a.logger.errors)
#         exit(1)
#     else:
#         a.run()
#     assert "Mpro-i0130" in [x.name for x in (Path(upload_2_dir) / "aligned_files").glob("*")]


# @pytest.mark.order(after="test_aligner_upload_2")
# def test_collator_upload_3(
#     constants,
#     test_dir,
#     upload_3_dir,
#     upload_3_data_dir,
#     config_3_file,
# ):
#     c = Collator(str(config_3_file), )
#
#     meta, num_errors, num_warnings = c.validate()
#
#     if meta is None or num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         c.run(meta)
#
#     with open(Path(upload_3_dir) / "meta_collator.yaml", 'r') as f:
#         new_meta = yaml.safe_load(f)
#     assert len(meta[Constants.META_XTALS]["Mpro-i0130"][Constants.META_XTAL_FILES].get(Constants.META_BINDING_EVENT, {})) != 0
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "crystallographic_files").glob("*")]
#
#
# @pytest.mark.order(after="test_collator_upload_3")
# def test_aligner_upload_3(constants, xtalforms_file, upload_3_dir):
#     a = Aligner(upload_3_dir, constants.METADATA_FILE, xtalforms_file, )
#     num_errors, num_warnings = a.validate()
#
#     if num_errors:
#         print("There are errors, cannot continue")
#         exit(1)
#     else:
#         a.run()
#
#     assert "5rgs" in [x.name for x in (Path(upload_3_dir) / "aligned_files").glob("*")]
