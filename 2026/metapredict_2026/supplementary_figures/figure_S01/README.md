# Figure S1 — metapredict V1 disorder vs AlphaFold2 pLDDT

metapredict V1 disorder scores and AF2 pLDDT disagree about which regions are disordered, which motivated combining them in V2/V3. Panel A is EIF4EBP1 (Q13541) and panel B is Myozenin-1 (Q9NP98); both are fully disordered in DisProt.

| File | Contents |
| --- | --- |
| `figure_S01A.ipynb`, `figure_S01B.ipynb` | One notebook per panel; each writes `figures/figure_S01<panel>.pdf`. |
| `disorder_vs_plddt.py` | Shared plotting function and the two thresholds (V1 0.42, pLDDT 70). |
| `figtools.py` | Shared manuscript style/colours (copy of `shared/figtools.py`). |

**Running:** open a notebook from this directory and run all cells. Requires metapredict 3.1.x, numpy and matplotlib. V1 scores are predicted live; the pLDDT values are read from `../../data/af2_plddt/` (see its README).

**Differences from the published panel**
- Residues called disordered are shaded as contiguous blocks rather than one grey line per residue; disorder calls use a pale V1-orange wash.
- pLDDT is drawn in black (it is not a metapredict version, and blue is now the V2 colour).
- Protein IDs sit above the top axes instead of being a title.
