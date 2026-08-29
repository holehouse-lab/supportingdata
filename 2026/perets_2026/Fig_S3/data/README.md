# Data for Figure S3

## `dssp/` - secondary structure occupancies (panel A)

Per-replica DSSP occupancies, one directory per simulation:

```
dssp/
├── helical_start/1 ... 10/    DSSP_H.csv, DSSP_E.csv
└── coil_start/1 ... 10/       DSSP_H.csv, DSSP_E.csv
```

Each CSV is a single column of 245 values, one per residue of the DDX4<sup>[1-245]</sup> construct, giving the **fraction of sampled conformations in which that residue was assigned to that DSSP class**:

| file | DSSP class |
|------|------------|
| `DSSP_H.csv` | helical |
| `DSSP_E.csv` | extended (beta strand / sheet) |

A third class, `DSSP_C.csv` (coil), exists in the source data but is not used by this figure; it is the remainder.

The 20 replicas split evenly into two sets by starting conformation - 10 begun from a fully helical structure and 10 from a random coil. This is a sampling control: because the two sets converge on the same per-residue occupancies, the reported secondary structure does not depend on where the simulation started. The violins in panel A show the spread across all 20.

Mean helical occupancy across the whole chain is 0.049 and mean extended occupancy is 0.010, which is the quantitative statement behind "DDX4<sup>[1-245]</sup> is disordered". The largest single feature is around residue 175, which reaches a mean helicity of 0.53.

These values were computed by CAMPARI during the production runs. They are not recomputed by the notebook.

## Trajectories (panels B and C)

Panels B and C are computed directly from the pooled CAMPARI trajectories, which are too large to track in git and therefore live in the shared `data/` directory at the root of `perets_2025`. See [`../../data/README.md`](../../data/README.md) for the full simulation provenance - engine, force field, simulation conditions, replica counts - and for what the self-avoiding random coil reference ensemble is and why it is the right comparison.

Both panels need `campari/ddx4_ntd/` and `campari/ev/`.

## References

Kabsch, W. & Sander, C. *Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features.* **Biopolymers** 22, 2577-2637 (1983). doi:[10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)

Vitalis, A. & Pappu, R. V. *ABSINTH: a new continuum solvation model for simulations of polypeptides in aqueous solutions.* **J. Comput. Chem.** 30, 673-699 (2009). doi:[10.1002/jcc.21005](https://doi.org/10.1002/jcc.21005)
