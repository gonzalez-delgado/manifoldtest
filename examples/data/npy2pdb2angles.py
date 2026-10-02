"""Compute the dihedral angles of every conformer in a nucleotide library .npy.

The angles are computed directly from the coordinates; the PDB template only
provides the atom names, residue ids and chain of each coordinate row.

Output: angles_<output>.csv with one row per residue per conformer
(conformer-major order), and optionally <output>.pdb with --pdb.
"""

# ----------------------- Imports -----------------------

import argparse
import numpy as np
import pandas as pd

# ----------------------- Parser -----------------------

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter
    )

parser.add_argument('npy', help = 'npy to convert')
parser.add_argument('template', help = 'template for the nucleotide')
parser.add_argument('-o', '--output', required = True, help = 'name of the outputs')
parser.add_argument('--pdb', action = 'store_true', help = 'also write <output>.pdb')

args = parser.parse_args()

# ----------------------- Angle definitions -----------------------

# Each angle is 4 atoms (residue offset, atom name), offset -1/0/+1 being the
# previous/current/next residue of the chain.

ANGLES = {
    # Backbone
    "alpha":   [(-1, "O3'"), (0, "P"),   (0, "O5'"), (0, "C5'")],
    "beta":    [(0, "P"),    (0, "O5'"), (0, "C5'"), (0, "C4'")],
    "gamma":   [(0, "O5'"),  (0, "C5'"), (0, "C4'"), (0, "C3'")],
    "delta":   [(0, "C5'"),  (0, "C4'"), (0, "C3'"), (0, "O3'")],
    "epsilon": [(0, "C4'"),  (0, "C3'"), (0, "O3'"), (1, "P")],
    "zeta":    [(0, "C3'"),  (0, "O3'"), (1, "P"),   (1, "O5'")],

    # Sugar
    "nu0":     [(0, "C4'"),  (0, "O4'"), (0, "C1'"), (0, "C2'")],
    "nu1":     [(0, "O4'"),  (0, "C1'"), (0, "C2'"), (0, "C3'")],
    "nu2":     [(0, "C1'"),  (0, "C2'"), (0, "C3'"), (0, "C4'")],
    "nu3":     [(0, "C2'"),  (0, "C3'"), (0, "C4'"), (0, "O4'")],
    "nu4":     [(0, "C3'"),  (0, "C4'"), (0, "O4'"), (0, "C1'")],

    # Pseudoangles
    "eta":     [(-1, "C4'"), (0, "P"),   (0, "C4'"), (1, "P")],
    "theta":   [(0, "P"),    (0, "C4'"), (1, "P"),   (1, "C4'")],
}

# Chi for Purines/Pyrimidines
CHI = {
    "A": [(0, "O4'"), (0, "C1'"), (0, "N9"), (0, "C4")],
    "G": [(0, "O4'"), (0, "C1'"), (0, "N9"), (0, "C4")],
    "C": [(0, "O4'"), (0, "C1'"), (0, "N1"), (0, "C2")],
    "U": [(0, "O4'"), (0, "C1'"), (0, "N1"), (0, "C2")],
}

# ----------------------- Functions -----------------------

def parse_template(template):
    """
    Return the ATOM/HETATM lines and the residues as
    [(chain, res_id, res_name, {atom name: coordinate row})].
    """
    lines = [l for l in open(template) if l.startswith(("ATOM", "HETATM"))]

    residues = []
    for i, l in enumerate(lines):
        key = (l[21], int(l[22:26]), l[17:20].strip())
        if not residues or residues[-1][:3] != key:
            residues.append(key + ({},))
        residues[-1][3][l[12:16].strip()] = i

    return lines, residues

def dihedral(a, b, c, d):
    """
    Dihedral angles in degrees in [0, 360) for arrays of points (N, 3),
    same convention as Bio.PDB.calc_dihedral.
    """
    ab, cb, db = a - b, c - b, d - c
    u, v = np.cross(ab, cb), np.cross(db, cb)
    w = np.cross(u, v)
    x = (u * v).sum(axis = 1)
    y = (w * cb).sum(axis = 1) / np.linalg.norm(cb, axis = 1)

    return np.degrees(np.arctan2(y, x)) % 360


def get_angle(coords, residues, r, atoms):
    """
    Angle over all conformers for residue index r, or None if an atom is
    missing (chain end or atom absent from the template).
    """
    rows = []
    for offset, name in atoms:
        if not 0 <= r + offset < len(residues):
            return None
        res = residues[r + offset]
        if res[0] != residues[r][0] or name not in res[3]:
            return None
        rows.append(res[3][name])

    return dihedral(*(coords[:, i] for i in rows))


def get_full_structure(coords, residues):
    """
    Compute all angles for all residues and return a DataFrame with one row per
    residue per conformer (conformer-major order).
    """
    n_models, n_res = coords.shape[0], len(residues)
    columns = {"chain": [], "res_id": [], "res_name": []}
    angles = {name: np.full((n_models, n_res), np.nan) for name in [*ANGLES, "chi"]}

    for r, (chain, res_id, res_name, _) in enumerate(residues):
        columns["chain"].append(chain)
        columns["res_id"].append(res_id)
        columns["res_name"].append(res_name)

        definitions = dict(ANGLES)
        if res_name[-1] in CHI:
            definitions["chi"] = CHI[res_name[-1]]

        for name, atoms in definitions.items():
            angle = get_angle(coords, residues, r, atoms)
            if angle is not None:
                angles[name][:, r] = angle

    # DataFrame creation, one row per (model, residue) in model-major order
    df = pd.DataFrame({k: np.tile(v, n_models) for k, v in columns.items()})
    for name, values in angles.items():
        df[name] = values.ravel()

    return df


def write_pdb(coords, lines, filename):
    """
    Write a PDB file with the coordinates of all conformers, one MODEL per conformer.
    """
    with open(filename, "w") as f:
        for m, model in enumerate(coords, start = 1):
            f.write(f"MODEL {m}\n")
            for l, (x, y, z) in zip(lines, model):
                f.write(l[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + l[54:].rstrip("\n") + "\n")
            f.write("ENDMDL\n")



# ----------------------- Output files -----------------------

coords = np.load(args.npy)
coords = coords.reshape(coords.shape[0], -1, 3)
lines, residues = parse_template(args.template)
assert len(lines) == coords.shape[1], (len(lines), coords.shape)

# Round to PDB precision and float32 (as Biopython stores them) so the angles
# match those read back from a PDB file.
coords = coords.round(3).astype(np.float32).astype(float)

get_full_structure(coords, residues).to_csv(f"angles_{args.output}.csv", index = False)

if args.pdb:
    write_pdb(coords, lines, f"{args.output}.pdb")
