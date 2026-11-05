# Figure S15 — metapredict V3 on synthetic and designed proteins

metapredict V3 predictions for 63 non-natural sequences from four studies.

| Panel | Notebook | Contents |
| --- | --- | --- |
| A, B, C | `figure_S15A/B/C.ipynb` | Disorder profiles of the hallucinated folded proteins 0738, 0515 and 0217 (Anishchenko et al. 2021). |
| D | `figure_S15D.ipynb` | Wild-type hnRNPA1 LCD and five designed variants (Pesce et al. 2024). |
| E | `figure_S15E.ipynb` | Fraction predicted disordered for 7 folded and 15 disordered random-sequence proteins (Tretyachenko et al. 2017). |
| F | `figure_S15F.ipynb` | Fraction predicted disordered for 32 GOOSE-designed IDRs (Emenecker, Guadalupe et al. 2023). |

`synthetic_disorder.py` holds the shared plotting and fraction-disordered functions, `figtools.py` is the shared style. The sequences for D–F are in `../../data/synthetic_proteins/` (see its README). The hallucinated sequences for A–C (designs 0738, 0515 and 0217 from Anishchenko *et al.*, *Nature* 600, 547–552 (2021), given as the expressed constructs; 0738 and 0217 carry a 20-residue N-terminal His-tag, so the native initiator Met is residue 21) are written in their notebooks.

**Running:** open a notebook from this directory and run all cells. Requires metapredict 3.1.x, protfasta, numpy and matplotlib. All predictions are made live. Each notebook writes `figures/figure_S15<panel>.pdf`.

**Result with metapredict 3.1.0:** 61 of 63 sequences are called correctly, as in the manuscript. Two disordered Tretyachenko sequences (4954, 9121) are predicted mostly folded, giving MCC = 0.821 for panel E. All 32 GOOSE sequences are predicted fully disordered.

**Differences from the published panel**
- Traces and IDR fills use the metapredict V3 colour. Experimental classes in E and F stay blue (folded) and red (disordered), and F is solid red instead of translucent.
- Sequence labels are drawn by the code instead of added at assembly, and D is one stacked figure instead of six separate PDFs. The His-tag arrows in A and C are still added at assembly.
- E no longer has a title; the MCC is printed instead. The unused AlphaFold3 input-JSON cell is dropped.
- F is drawn vertically and rotated 90° when the figure is assembled, as before.
