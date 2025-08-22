def add_pyscf_to_dict(cls_dict):
    if "pyscf" in cls_dict:
        return

    try:
        from pysisyphus.calculators import PySCF

        cls_dict["pyscf"] = PySCF.PySCF
    except (ModuleNotFoundError, OSError):
        print("PySCF import failed! Did you forget to install it?")
