from pathlib import Path

import pytest

from gemmi import cif

from pdbdepo import scrape_processing_stats as sps

# Fixtures are real excerpts grabbed from data/std_test processing output, trimmed to just the
# summary section (or, for the small .table1 file, kept whole) that each parser scans.
DATA_DIR = Path(__file__).parent.parent / "test-data" / "processing_stats"
AUTOPROC_LOG = DATA_DIR / "autoproc_aimless.log"
AUTOPROC_STARANISO_TABLE1 = DATA_DIR / "staraniso_alldata-unique.table1"
XIA_3DII_LOG = DATA_DIR / "xia_3dii.log"
XIA2_MULTIPLEX_HTML = DATA_DIR / "xia2_multiplex.html"
# phenix.merging_statistics log for A71EV2A-x5488, from m2ms/fragalysis-frontend#2370
PHENIX_LOG = DATA_DIR / "phenix_merging_stats.log"


# ---------------------------------------------------------------------------
# handle_autoproc
# ---------------------------------------------------------------------------


def test_handle_autoproc_overall_stats():
    reflns, shell = sps.handle_autoproc(str(AUTOPROC_LOG))

    assert reflns["d_resolution_low"] == "46.38"
    assert reflns["d_resolution_high"] == "1.34"
    assert reflns["number_measured_obs"] == "210486"
    assert reflns["number_obs"] == "28848"
    assert reflns["pdbx_redundancy"] == "7.3"
    assert reflns["pdbx_netI_over_sigmaI"] == "6.5"
    assert reflns["pdbx_Rmerge_I_obs"] == "0.173"
    assert reflns["pdbx_Rrim_I_all"] == "0.186"
    assert reflns["pdbx_Rpim_I_all"] == "0.068"
    assert reflns["pdbx_CC_half"] == "0.997"
    assert reflns["percent_possible_obs"] == "92.5"
    # no chi-squared regex for autoproc (unlike xia_3dii)
    assert "pdbx_chi_squared" not in reflns


def test_handle_autoproc_shell_stats():
    reflns, shell = sps.handle_autoproc(str(AUTOPROC_LOG))

    # first row is outer (high-resolution) shell, second is inner (low-resolution) shell
    assert shell["d_res_low"] == ("1.41", "46.38")
    assert shell["d_res_high"] == ("1.34", "4.24")
    assert shell["number_measured_obs"] == ("27769", "7398")
    assert shell["number_unique_obs"] == ("3860", "1019")
    assert shell["pdbx_CC_half"] == ("0.376", "0.998")


def test_handle_autoproc_missing_file():
    reflns, shell = sps.handle_autoproc(str(DATA_DIR / "does_not_exist.log"))
    assert not reflns
    assert not shell


def test_handle_autoproc_none_file():
    reflns, shell = sps.handle_autoproc(None)
    assert not reflns
    assert not shell


# ---------------------------------------------------------------------------
# handle_autoproc_staraniso
# ---------------------------------------------------------------------------


def test_handle_autoproc_staraniso_overall_stats():
    reflns, shell = sps.handle_autoproc_staraniso(str(AUTOPROC_STARANISO_TABLE1))

    assert reflns["d_resolution_low"] == "36.496"
    assert reflns["d_resolution_high"] == "1.530"
    assert reflns["number_measured_obs"] == "49377"
    assert reflns["number_obs"] == "13779"
    assert reflns["pdbx_redundancy"] == "3.6"
    assert reflns["pdbx_netI_over_sigmaI"] == "6.3"
    assert reflns["pdbx_Rmerge_I_obs"] == "0.086"
    assert reflns["pdbx_Rrim_I_all"] == "0.095"
    assert reflns["pdbx_Rpim_I_all"] == "0.040"
    assert reflns["pdbx_CC_half"] == "0.996"
    # staraniso uses the ellipsoidal completeness figure
    assert reflns["percent_possible_obs"] == "74.4"


def test_handle_autoproc_staraniso_shell_stats():
    reflns, shell = sps.handle_autoproc_staraniso(str(AUTOPROC_STARANISO_TABLE1))

    assert shell["d_res_low"] == ("1.739", "36.496")
    assert shell["d_res_high"] == ("1.530", "4.913")
    assert shell["number_measured_obs"] == ("1312", "4384")
    assert shell["number_unique_obs"] == ("689", "689")


# ---------------------------------------------------------------------------
# handle_xia_3dii
# ---------------------------------------------------------------------------


def test_handle_xia_3dii_overall_stats():
    reflns, shell = sps.handle_xia_3dii(str(XIA_3DII_LOG))

    assert reflns["d_resolution_low"] == "47.74"
    assert reflns["d_resolution_high"] == "1.60"
    assert reflns["number_measured_obs"] == "147405"
    assert reflns["number_obs"] == "21643"
    assert reflns["pdbx_redundancy"] == "6.8"
    assert reflns["pdbx_netI_over_sigmaI"] == "4.6"
    assert reflns["pdbx_Rmerge_I_obs"] == "0.183"
    assert reflns["pdbx_Rrim_I_all"] == "0.199"
    assert reflns["pdbx_Rpim_I_all"] == "0.076"
    assert reflns["pdbx_CC_half"] == "0.984"
    assert reflns["percent_possible_obs"] == "100.0"
    # xia_3dii (unlike autoproc) also captures Mean(Chi^2)
    assert reflns["pdbx_chi_squared"] == "0.70"


def test_handle_xia_3dii_shell_stats():
    reflns, shell = sps.handle_xia_3dii(str(XIA_3DII_LOG))

    assert shell["d_res_low"] == ("1.69", "47.74")
    assert shell["d_res_high"] == ("1.60", "5.07")
    assert shell["pdbx_chi_squared"] == ("0.54", "0.44")


# ---------------------------------------------------------------------------
# handle_xia2_multiplex
# ---------------------------------------------------------------------------


def test_handle_xia2_multiplex_overall_stats():
    reflns, shell = sps.handle_xia2_multiplex(str(XIA2_MULTIPLEX_HTML))

    assert reflns["d_resolution_low"] == "87.11"
    assert reflns["d_resolution_high"] == "1.89"
    assert reflns["number_measured_obs"] == "358302"
    assert reflns["number_obs"] == "30629"
    assert reflns["pdbx_redundancy"] == "11.7"
    assert reflns["pdbx_netI_over_sigmaI"] == "14.6"
    assert reflns["pdbx_Rmerge_I_obs"] == "0.282"
    assert reflns["pdbx_Rrim_I_all"] == "0.295"
    assert reflns["pdbx_Rpim_I_all"] == "0.084"
    assert reflns["pdbx_CC_half"] == "0.998"
    # the '%' suffix must be stripped
    assert reflns["percent_possible_obs"] == "100.00"


def test_handle_xia2_multiplex_shell_stats():
    reflns, shell = sps.handle_xia2_multiplex(str(XIA2_MULTIPLEX_HTML))

    # first row is outer (high resolution) shell, second is inner (low resolution) shell
    assert shell["d_res_low"] == ("1.92", "87.20")
    assert shell["d_res_high"] == ("1.89", "5.13")
    assert shell["number_measured_obs"] == ("12297", "18446")
    assert shell["number_unique_obs"] == ("1503", "1585")
    assert shell["pdbx_CC_half"] == ("0.351", "0.998")
    assert shell["percent_possible_obs"] == ("99.60", "100.00")


def test_handle_xia2_multiplex_missing_file():
    reflns, shell = sps.handle_xia2_multiplex(str(DATA_DIR / "does_not_exist.html"))
    assert reflns == {"entry_id": "UNNAMED", "pdbx_diffrn_id": 1, "pdbx_ordinal": 1}
    assert shell == {"pdbx_diffrn_id": (1, 1), "pdbx_ordinal": (1, 2)}


def test_handle_xia2_multiplex_none_file():
    reflns, shell = sps.handle_xia2_multiplex(None)
    assert reflns == {"entry_id": "UNNAMED", "pdbx_diffrn_id": 1, "pdbx_ordinal": 1}
    assert shell == {"pdbx_diffrn_id": (1, 1), "pdbx_ordinal": (1, 2)}


def test_handle_xia2_multiplex_missing_overall_panel(tmp_path):
    """An HTML file without the expected panel/table shouldn't crash, just warn and return
    the bare reflns/shell dicts."""
    html_file = tmp_path / "no_panel.html"
    html_file.write_text("<html><body><p>no stats here</p></body></html>")

    reflns, shell = sps.handle_xia2_multiplex(str(html_file))

    assert reflns == {"entry_id": "UNNAMED", "pdbx_diffrn_id": 1, "pdbx_ordinal": 1}
    assert shell == {"pdbx_diffrn_id": (1, 1), "pdbx_ordinal": (1, 2)}


# ---------------------------------------------------------------------------
# handle_phenix
# ---------------------------------------------------------------------------


def test_handle_phenix_overall_stats():
    reflns, shell = sps.handle_phenix(str(PHENIX_LOG))

    assert reflns["d_resolution_low"] == "47.03"
    assert reflns["d_resolution_high"] == "1.89"
    assert reflns["number_measured_obs"] == "78692"
    assert reflns["number_obs"] == "11633"
    assert reflns["pdbx_redundancy"] == "6.76"
    assert reflns["pdbx_netI_over_sigmaI"] == "3.8"
    assert reflns["pdbx_Rmerge_I_obs"] == "0.909"
    assert reflns["pdbx_Rrim_I_all"] == "0.985"
    assert reflns["pdbx_Rpim_I_all"] == "0.377"
    assert reflns["pdbx_CC_half"] == "0.903"
    assert reflns["percent_possible_obs"] == "99.98"
    assert "pdbx_chi_squared" not in reflns


def test_handle_phenix_shell_stats():
    reflns, shell = sps.handle_phenix(str(PHENIX_LOG))

    # first row is outer (high-resolution) shell, second is inner (low-resolution) shell
    assert shell["d_res_low"] == ("1.92", "47.04")
    assert shell["d_res_high"] == ("1.89", "5.13")
    assert shell["number_measured_obs"] == ("4039", "3794")
    assert shell["number_unique_obs"] == ("581", "608")
    assert shell["pdbx_Rmerge_I_obs"] == ("2.471", "0.402")
    assert shell["pdbx_CC_half"] == ("0.151", "0.890")


def test_handle_phenix_missing_file():
    reflns, shell = sps.handle_phenix(str(DATA_DIR / "does_not_exist.log"))
    assert not reflns
    assert not shell


def test_handle_phenix_none_file():
    reflns, shell = sps.handle_phenix(None)
    assert not reflns
    assert not shell


def test_handle_phenix_no_table(tmp_path):
    log = tmp_path / "phenix.log"
    log.write_text("Merging statistics\nResolution: 47.03 - 1.89\n")
    reflns, shell = sps.handle_phenix(str(log))
    assert "d_resolution_high" not in reflns
    assert "d_res_high" not in shell


@pytest.mark.parametrize(
    "type_name",
    [
        "autoproc-phenix",
        "autoproc-staraniso-phenix",
        "xia2-3dii-phenix",
        "xia2-dials-phenix",
        "xia2_dials_phenix",
        "autoproc_staraniso_phenix",
    ],
)
def test_is_phenix_type(type_name):
    assert sps.is_phenix_type(type_name)


@pytest.mark.parametrize("type_name", ["autoproc", "xia2-dials", "xia2-multiplex-phenix", None, ""])
def test_is_not_phenix_type(type_name):
    assert not sps.is_phenix_type(type_name)


# ---------------------------------------------------------------------------
# handle_file — dispatch by type, and CIF document assembly
# ---------------------------------------------------------------------------


def test_handle_file_autoproc_creates_reflns_loop():
    doc = sps.handle_file(str(AUTOPROC_LOG), "autoproc", None, None)
    block = doc[0]
    assert block.find_pair("_reflns.d_resolution_high") == ("_reflns.d_resolution_high", "1.34")
    shell_item = block.find_loop_item("_reflns_shell.d_res_high")
    assert list(shell_item.loop.values)[list(shell_item.loop.tags).index("_reflns_shell.d_res_high")] == "1.34"


def test_handle_file_autoproc_staraniso_creates_reflns_loop():
    doc = sps.handle_file(str(AUTOPROC_STARANISO_TABLE1), "autoproc_staraniso", None, None)
    block = doc[0]
    assert block.find_pair("_reflns.d_resolution_high") == ("_reflns.d_resolution_high", "1.530")


def test_handle_file_xia_3dii_creates_reflns_loop():
    doc = sps.handle_file(str(XIA_3DII_LOG), "xia2-3dii", None, None)
    block = doc[0]
    assert block.find_pair("_reflns.d_resolution_high") == ("_reflns.d_resolution_high", "1.60")


def test_handle_file_xia2_multiplex_creates_reflns_loop():
    doc = sps.handle_file(str(XIA2_MULTIPLEX_HTML), "xia2-multiplex", None, None)
    block = doc[0]
    assert block.find_pair("_reflns.d_resolution_high") == ("_reflns.d_resolution_high", "1.89")
    assert block.find_pair("_reflns.percent_possible_obs") == ("_reflns.percent_possible_obs", "100.00")


def test_handle_file_appends_to_existing_doc():
    """handle_file should add its new block to an existing doc rather than replacing it."""
    existing_doc = cif.Document()
    existing_doc.add_new_block("existing")

    doc = sps.handle_file(str(XIA2_MULTIPLEX_HTML), "xia2-multiplex", existing_doc, None)

    assert len(doc) == 2
    assert doc[0].name == "existing"


def test_handle_file_unsupported_type_returns_none():
    assert sps.handle_file(str(AUTOPROC_LOG), "not-a-real-type", None, None) is None


def test_handle_file_writes_output_file(tmp_path):
    out_file = tmp_path / "stats.cif"
    sps.handle_file(str(XIA2_MULTIPLEX_HTML), "xia2-multiplex", None, str(out_file))
    assert out_file.is_file()
    content = out_file.read_text()
    assert "_reflns.d_resolution_high" in content


@pytest.mark.parametrize(
    "type_name",
    ["autoproc-phenix", "autoproc-staraniso-phenix", "xia2-3dii-phenix", "xia2-dials-phenix", "xia2_3dii_phenix"],
)
def test_handle_file_phenix_creates_reflns_loop(type_name):
    doc = sps.handle_file(str(PHENIX_LOG), type_name, None, None)
    block = doc[0]
    assert block.find_pair("_reflns.d_resolution_high") == ("_reflns.d_resolution_high", "1.89")
    assert block.find_loop("_reflns_shell.pdbx_CC_half")[0] == "0.151"
