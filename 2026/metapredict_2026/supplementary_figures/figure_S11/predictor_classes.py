"""Classification of CAID3 predictor flavors into disorder vs binding/linker.

The CAID3 disorder ranking tables (Tables 1-2 of Mehdiabadi et al., 2025) only
include predictors that target intrinsic disorder. Many submitted flavors
instead predict binding sites (protein/RNA/DNA/ion) or linker regions; these are
evaluated on the binding and linker references, not the disorder references.

Every flavor still produces a per-residue score, so they *can* be scored against
the disorder reference, but including them in a disorder benchmark is misleading.
This module lists the binding/linker flavors to exclude so the remaining set
mirrors the CAID3 disorder challenge. The list was built from the CAID3 tables
and the published descriptions of each method, then verified by checking that the
filtered top-10 reproduces the paper's Disorder-PDB table.
"""

from __future__ import annotations

# Flavors that predict binding or linker regions rather than intrinsic disorder.
# These are excluded from the disorder ranking.
BINDING_LINKER_FLAVORS = {
    # AlphaFold binding heads (the rsa / pLDDT heads are disorder proxies, kept)
    "AlphaFold-binding",
    "AlphaFold3-binding",
    # Binding-site predictors
    "AIUPred-2-binding",
    "ANCHOR2",
    "bindEmbed21IDR-idrGeneral",
    "bindEmbed21IDR-idrNuc",
    "bindEmbed21IDR-rawGeneral",
    "bindEmbed21IDR-rawNuc",
    "DeepDISObind-all",
    "DeepDISObind-nucleic",
    "DeepDISObind-protein",
    "DeepDRPBind-DNA",
    "DeepDRPBind-protein",
    "DeepDRPBind-RNA",
    "DISOPRED3-bind",
    "DisoRDPbind-dna",
    "DisoRDPbind-prot",
    "DisoRDPbind-rna",
    "DRPBind-DNA",
    "DRPBind-protein",
    "DRPBind-RNA",
    "EBIND-nucleic-acid",
    "EBIND-protein",
    "ENSHROUD-all",
    "ENSHROUD-nucleic",
    "ENSHROUD-protein",
    "IPA-AF2-Protein-bind",
    "IPA-All-bind",
    "IPA-Nucleotide-bind",
    "IPA-Protein-bind",
    "MoRFchibi",
    "MoRFchibi-light",
    "MoRFchibi-web",
    "NeProc-binding",
    "NeProc2-beta-binding",
    "ProBiPred-nucleic",
    "ProBiPred-protein",
    # DisoFLAG binding/linker heads (only DisoFLAG-IDR targets disorder)
    "DisoFLAG-DB",
    "DisoFLAG-IB",
    "DisoFLAG-LB",
    "DisoFLAG-PB",
    "DisoFLAG-RB",
    "DisoFLAG-DFL",
    # Linker / flexible-linker predictors
    "APOD",
    "CLIP",
    "DFLpred",
    "IPA-AF2-Linker",
    "LINKER-Pred",
    "LINKER-Pred-Lite",
    "LINKER-Pred2",
    "LIPNet",
}


def is_disorder_predictor(name: str) -> bool:
    """Return ``True`` if a flavor targets intrinsic disorder (not binding/linker)."""
    return name not in BINDING_LINKER_FLAVORS
