script = """
if __name__ == "__main__":
    #p = read_pdb("test_output_6/aligned_structures/1/1/Zika_NS3A-A0186_0.pdb")
    #fix_nomenclature_errors(p)
    set_nomenclature_errors_on_read("ignore")
    p = read_pdb("test_output_6/aligned_structures/1/1/Zika_NS3A-A0186_0.pdb")
"""
