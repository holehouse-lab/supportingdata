"""Download the AlphaFold2 models of the Figure S10 examples and render them with PyMOL.

For each example protein (``short_idr_examples.select_examples``) this script

1. downloads the canonical AlphaFold DB model (``AF-<UniProt>-F1``) into
   ``data/alphafold_models/`` (skipped if already there);
2. checks that the model's sequence is identical to the CAID3 sequence and reports
   how closely its pLDDT matches CAID3's AlphaFold-pLDDT values (CAID3 used an
   earlier AlphaFold DB release; only the current one can be downloaded);
3. renders a cartoon coloured by the standard AlphaFold pLDDT bands, with the short
   IDRs metapredict V3 misses drawn as thick, opaque tubes and the rest of the
   protein translucent, to ``data/alphafold_renders/<DisProt ID>.png``.

It writes ``data/alphafold_renders/models.tsv`` with the model used for each example
and the pLDDT agreement. Both ``data/`` directories are at the top of the repository.
Needs network access on the first run and PyMOL (``pymol-open-source``); the notebooks
only read the PNGs.

Run from this directory: ``python render_structures.py``
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from pymol import cmd

from caid3_panels import load_predictions, load_reference
from figtools import SHARED_DATA_DIR, require_metapredict_3_1
from ground_truth_panels import get_uniprot
from short_idr_examples import (
    ALPHAFOLD_BANDS,
    PANEL_LETTERS,
    RENDER_DIR,
    Example,
    plddt_trace,
    render_path,
    select_examples,
)

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
# AlphaFold DB refuses requests without a browser-like User-Agent.
HTTP_HEADERS = {"User-Agent": "Mozilla/5.0"}
HTTP_TIMEOUT_S = 60
MODEL_DIR = SHARED_DATA_DIR / "alphafold_models"
MODELS_TSV = RENDER_DIR / "models.tsv"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Rendering
IMAGE_SIZE_PX = 2000  # ~900 dpi across the 2.2-inch structure column
RAY_ANTIALIAS = 2  # PyMOL ray-tracing supersampling level (default 1)
IMAGE_DPI = 300
CORE_MIN_PLDDT = 50.0  # residues used to orient and frame the model; lower-confidence tails may leave the frame
ZOOM_BUFFER_A = 2.0
REST_TRANSPARENCY = 0.45  # cartoon transparency of everything except the missed short IDRs
MISSED_TUBE_RADIUS_A = 1.1  # loop/tube radius of the missed short IDRs (PyMOL default 0.2)
MISSED_OVAL_WIDTH_A = 0.6  # helix thickness of the missed short IDRs (PyMOL default 0.25)
MISSED_RECT_WIDTH_A = 0.8  # strand thickness of the missed short IDRs (PyMOL default 0.4)
# 180-degree flips that keep the principal axes on the screen axes; the one that
# brings the missed short IDRs closest to the viewer is used.
AXIS_FLIPS = [np.diag(signs) for signs in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]


def fetch_model(accession: str) -> Tuple[Path, int]:
    """Download the canonical AlphaFold DB model of a UniProt accession, if not cached.

    Parameters
    ----------
    accession : str
        UniProt accession.

    Returns
    -------
    tuple of (pathlib.Path, int)
        Path of the PDB file and the AlphaFold DB version.
    """
    request = urllib.request.Request(AFDB_API.format(accession=accession), headers=HTTP_HEADERS)
    with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
        entries = json.load(response)
    canonical = [e for e in entries if e["entryId"] == f"AF-{accession}-F1"]
    assert len(canonical) == 1, f"expected one AF-{accession}-F1 entry, found {len(canonical)}"
    entry = canonical[0]
    path = MODEL_DIR / entry["pdbUrl"].rsplit("/", 1)[1]
    if not path.exists():
        MODEL_DIR.mkdir(exist_ok=True)
        request = urllib.request.Request(entry["pdbUrl"], headers=HTTP_HEADERS)
        with urllib.request.urlopen(request, timeout=HTTP_TIMEOUT_S) as response:
            path.write_bytes(response.read())
        print(f"downloaded {path.name}")
    return path, int(entry["latestVersion"])


def model_residues(path: Path) -> Tuple[str, np.ndarray]:
    """Sequence and per-residue pLDDT (the B-factor of each CA atom) of an AlphaFold model."""
    residues: Dict[int, Tuple[str, float]] = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                residues[int(line[22:26])] = (THREE_TO_ONE[line[17:20]], float(line[60:66]))
    numbers = sorted(residues)
    assert numbers == list(range(1, len(numbers) + 1)), f"{path.name}: residues are not numbered 1..N"
    return "".join(residues[n][0] for n in numbers), np.array([residues[n][1] for n in numbers])


def check_model(example: Example, sequence: str, caid_plddt: np.ndarray, path: Path) -> Dict[str, float]:
    """Require the model sequence to match CAID3 and summarise pLDDT agreement with CAID3.

    Returns
    -------
    dict
        Median and maximum absolute pLDDT difference over the whole protein and over
        the missed short IDRs, and the Pearson correlation.
    """
    model_sequence, model_plddt = model_residues(path)
    assert model_sequence == sequence, f"{example.target}: AlphaFold model sequence differs from CAID3"
    difference = np.abs(model_plddt - caid_plddt)
    in_missed = np.zeros(len(sequence), dtype=bool)
    for idr in example.missed:
        in_missed[idr.start:idr.end] = True
    return {
        "median_abs_diff": float(np.median(difference)),
        "max_abs_diff": float(difference.max()),
        "max_abs_diff_missed": float(difference[in_missed].max()),
        "pearson_r": float(np.corrcoef(model_plddt, caid_plddt)[0, 1]),
    }


def colour_by_plddt(selection: str) -> None:
    """Colour a selection by the standard AlphaFold pLDDT bands (pLDDT is in the B-factor)."""
    for name, lower, upper, _label, _hex in ALPHAFOLD_BANDS:
        cmd.color(name, f"({selection}) and b >= {lower} and b < {upper}")


def model_ca_atoms(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Residue numbers, pLDDT and coordinates (A) of the CA atoms of an AlphaFold model.

    Read from the PDB file rather than through PyMOL, because the PyMOL wheel's
    NumPy-returning calls (``get_coords``) are built against NumPy 1.x.
    """
    numbers, plddt, coords = [], [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                numbers.append(int(line[22:26]))
                plddt.append(float(line[60:66]))
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(numbers), np.array(plddt), np.array(coords)


def view_transform(example: Example, model_path: Path) -> List[float]:
    """Rotation that puts the core on the screen axes with the missed short IDRs toward the viewer.

    After this transform and ``cmd.reset``, the screen x/y axes are the first/second
    principal axes of the core (CA atoms with pLDDT >= ``CORE_MIN_PLDDT``) and +z
    points at the viewer. Of the four rotations that keep this alignment, the one
    that puts the missed short IDRs' CA atoms furthest toward +z is used.

    Returns
    -------
    list of float
        Row-major 4x4 homogeneous matrix, for ``cmd.transform_selection(..., homogenous=1)``.
    """
    numbers, plddt, coords = model_ca_atoms(model_path)
    core = coords[plddt >= CORE_MIN_PLDDT]
    centre = core.mean(axis=0)
    _u, _s, axes = np.linalg.svd(core - centre)
    if np.linalg.det(axes) < 0:
        axes[2] *= -1  # keep a proper rotation
    in_missed = np.zeros(numbers.size, dtype=bool)
    for idr in example.missed:
        in_missed |= (numbers >= idr.start + 1) & (numbers <= idr.end)
    missed_offsets = coords[in_missed] - centre
    best_flip = max(AXIS_FLIPS, key=lambda flip: float(np.mean((missed_offsets @ (flip @ axes).T)[:, 2])))
    rotation = best_flip @ axes
    translation = -rotation @ centre
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation
    return [float(value) for value in matrix.ravel()]


def render(example: Example, model_path: Path, out_png: Path) -> None:
    """Render one AlphaFold model to PNG (see the module docstring for the style)."""
    cmd.reinitialize()
    for name, _lower, _upper, _label, hex_colour in ALPHAFOLD_BANDS:
        cmd.set_color(name, [int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    cmd.load(str(model_path), "afmodel")
    cmd.transform_selection("afmodel", view_transform(example, model_path), homogenous=1)
    missed = "afmodel and resi " + "+".join(f"{idr.start + 1}-{idr.end}" for idr in example.missed)

    cmd.hide("everything")
    cmd.show("cartoon", "afmodel")
    colour_by_plddt("afmodel")
    cmd.set("cartoon_transparency", REST_TRANSPARENCY, "afmodel")
    cmd.create("missed", missed)
    cmd.hide("everything", "missed")
    cmd.show("cartoon", "missed")
    colour_by_plddt("missed")
    cmd.set("cartoon_loop_radius", MISSED_TUBE_RADIUS_A, "missed")
    cmd.set("cartoon_tube_radius", MISSED_TUBE_RADIUS_A, "missed")
    cmd.set("cartoon_oval_width", MISSED_OVAL_WIDTH_A, "missed")
    cmd.set("cartoon_rect_width", MISSED_RECT_WIDTH_A, "missed")

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.15)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("antialias", RAY_ANTIALIAS)
    cmd.reset()
    cmd.zoom(f"(afmodel and b >= {CORE_MIN_PLDDT}) or missed", ZOOM_BUFFER_A, complete=1)
    out_png.parent.mkdir(exist_ok=True)
    cmd.png(str(out_png), width=IMAGE_SIZE_PX, height=IMAGE_SIZE_PX, dpi=IMAGE_DPI, ray=1)
    print(f"wrote {out_png}")


def main() -> None:
    require_metapredict_3_1()
    reference = load_reference()
    predictions = load_predictions(["Metapredict-v3", "AlphaFold-pLDDT"])
    examples = select_examples(reference, predictions["Metapredict-v3"])

    rows: List[str] = ["panel\ttarget\tuniprot\tmodel_file\tafdb_version\tmedian_abs_plddt_diff\t"
                       "max_abs_plddt_diff\tmax_abs_plddt_diff_missed\tpearson_r"]
    for letter, example in zip(PANEL_LETTERS, examples):
        sequence, labels = reference[example.target]
        accession = get_uniprot(example.target)[0]
        model_path, version = fetch_model(accession)
        caid_plddt = plddt_trace(labels, predictions["AlphaFold-pLDDT"][example.target])
        agreement = check_model(example, sequence, caid_plddt, model_path)
        print(f"{letter} {example.target} {accession} v{version}: pLDDT vs CAID3 median |diff| "
              f"{agreement['median_abs_diff']:.2f}, max {agreement['max_abs_diff']:.1f} "
              f"(missed short IDRs {agreement['max_abs_diff_missed']:.1f}), r {agreement['pearson_r']:.3f}")
        render(example, model_path, render_path(example.target))
        rows.append(f"{letter}\t{example.target}\t{accession}\t{model_path.name}\t{version}\t"
                    f"{agreement['median_abs_diff']:.2f}\t{agreement['max_abs_diff']:.1f}\t"
                    f"{agreement['max_abs_diff_missed']:.1f}\t{agreement['pearson_r']:.4f}")
    MODELS_TSV.write_text("\n".join(rows) + "\n")
    print(f"wrote {MODELS_TSV}")


if __name__ == "__main__":
    main()
