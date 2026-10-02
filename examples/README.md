# Examples

Two worked examples of `manifold.dcov.test()` on molecular torsion angles.
Each one has a script that downloads its data into `data/`, and an R script
that runs the tests and writes its results to `results/`.

| Example | Data script | Analysis script | Output in `results/` |
|---|---|---|---|
| RNA trinucleotides | `data/build_nucleotides.sh` | `trinucleotides.R` | `trinucleotide_pvalues.txt` (LaTeX table) |
| Tripeptides | `data/build_tripeptides.sh` | `tripeptides.R` | `tripeptide_pvalues.pdf` (plots) |

## 1. Download the data

```bash
cd examples/data
./build_nucleotides.sh    # -> data/nucleotides.npy   (~1 min)
./build_tripeptides.sh    # -> data/tripeptides/      (~580 MB download)
```

- `build_nucleotides.sh` downloads the 0.5 Å trinucleotide conformer libraries
  ([nucleotide-library](https://github.com/sjdv1982/nucleotide-library)) and
  their PDB templates
  ([nucleotide-fragment-templates](https://github.com/sjdv1982/nucleotide-fragment-templates)).
  Then `npy2pdb2angles.py` computes the dihedral angles, and the script saves
  them as `nucleotides.npy`.
- `build_tripeptides.sh` downloads the tripeptide angle archive from
  [moma.laas.fr](https://moma.laas.fr/static/data/tripeptide_angles_data.tar)
  and keeps only its `data_coil/` folder, saved as `tripeptides/`.

## 2. Run the analyses

Run the R scripts from `examples/`, because they use paths relative to it:

```bash
cd examples
Rscript trinucleotides.R
Rscript tripeptides.R
```

The tests use permutations and can take a while. Each script uses 3 cores by
default. To change that, edit `n_cores` at the top of the script.

The p-values are cached in `results/*.Rda`. If those files exist, a script
loads them and only regenerates the table or plots. To rerun the tests, delete
the cached files first:

```bash
rm results/trinucleotide_pvalues.Rda                                    # trinucleotides
rm results/tripeptide_pvalues_H0.Rda results/tripeptide_pvalues_H1.Rda  # tripeptides
```
