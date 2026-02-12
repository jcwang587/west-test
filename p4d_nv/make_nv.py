#!/usr/bin/env python3
"""
Generate QE + WEST input files for NV center in diamond:
- pw.in
- wstat.in

Requirements:
  pip install ase numpy
"""

import numpy as np
from ase import Atoms
from ase.build import bulk


def pairwise_distance_pbc(cell, r1, r2):
    """
    Minimum-image distance under PBC in a general cell.
    """
    # fractional coords
    inv_cell = np.linalg.inv(cell.T)
    s1 = inv_cell @ r1
    s2 = inv_cell @ r2
    ds = s2 - s1
    ds -= np.round(ds)  # wrap to [-0.5,0.5)
    dr = cell.T @ ds
    return np.linalg.norm(dr)


def build_nv_diamond(
    supercell_n=3,  # 3 -> 216 pristine atoms, 4 -> 512 pristine atoms
    a0=3.567,  # diamond lattice constant in angstrom (room-temp-ish)
    charge_state="minus",  # "minus" -> NV-, "neutral" -> NV0
    ecutwfc=50.0,
    pseudo_C="C_ONCV_PBE-1.0.upf",
    pseudo_N="N_ONCV_PBE-1.0.upf",
    outdir="./",
    pseudo_dir="./",
    prefix="nv_diamond",
    n_pdep_eigen=300,
):
    # Build conventional cubic diamond cell (8 atoms), then supercell
    # ASE bulk('C','diamond',cubic=True) gives conventional cubic diamond cell
    base = bulk("C", "diamond", a=a0, cubic=True)
    atoms = base.repeat((supercell_n, supercell_n, supercell_n))

    # Choose a carbon atom near center for substitution (to maximize symmetry distance to images)
    cell = atoms.cell.array
    center = np.mean(cell, axis=0) / 2.0  # approximate geometric center
    positions = atoms.get_positions()

    d_to_center = np.linalg.norm(positions - center, axis=1)
    iN = int(np.argmin(d_to_center))  # index of C -> N
    rN = positions[iN].copy()

    # Find nearest-neighbor carbon to iN for vacancy
    # In diamond, nearest-neighbor distance ~ a*sqrt(3)/4
    nn_candidates = []
    for j, rj in enumerate(positions):
        if j == iN:
            continue
        d = pairwise_distance_pbc(cell, rN, rj)
        nn_candidates.append((d, j))
    nn_candidates.sort(key=lambda x: x[0])

    # pick closest atom as vacancy site
    ivac = nn_candidates[0][1]

    # Apply defect: substitute N and remove neighbor C
    symbols = atoms.get_chemical_symbols()
    symbols[iN] = "N"

    # remove vacancy atom
    keep = [k for k in range(len(atoms)) if k != ivac]
    atoms_def = atoms[keep]
    symbols_def = [symbols[k] for k in keep]
    atoms_def.set_chemical_symbols(symbols_def)

    # Count atoms after defect
    nat = len(atoms_def)
    ntyp = 2

    # Charge
    if charge_state.lower() in ("minus", "nv-", "-1"):
        tot_charge = -1.0
    elif charge_state.lower() in ("neutral", "nv0", "0"):
        tot_charge = 0.0
    else:
        raise ValueError("charge_state must be 'minus' or 'neutral'")

    # K-point grid
    kpoints_block = "K_POINTS gamma"

    # Prepare QE text
    cell_lines = "\n".join(
        f"  {v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}" for v in atoms_def.cell.array
    )

    pos_lines = []
    for s, r in zip(atoms_def.get_chemical_symbols(), atoms_def.get_positions()):
        pos_lines.append(f"{s:2s} {r[0]:16.10f} {r[1]:16.10f} {r[2]:16.10f}")
    pos_block = "\n".join(pos_lines)

    # Note:
    # - assume_isolated is left as 'none' for bulk 3D defect supercell benchmarking.
    # - For charged-defect energetics, finite-size corrections may be needed separately.
    pw_in = f"""&CONTROL
  calculation = 'scf'
  prefix = '{prefix}'
  outdir = '{outdir}'
  pseudo_dir = '{pseudo_dir}'
  verbosity = 'high'
/
&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = {ecutwfc}
  occupations = 'fixed'
  nspin = 1
  tot_charge = {tot_charge}
  input_dft = 'PBE'
  assume_isolated = 'none'
/
&ELECTRONS
  diago_full_acc = .false.
/

ATOMIC_SPECIES
C  12.011  {pseudo_C}
N  14.007  {pseudo_N}

CELL_PARAMETERS angstrom
{cell_lines}

ATOMIC_POSITIONS angstrom
{pos_block}

{kpoints_block}
"""

    wstat_in = f"""input_west:
  qe_prefix: {prefix}
  west_prefix: {prefix}
  outdir: {outdir}

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: {n_pdep_eigen}
"""

    with open("pw.in", "w") as f:
        f.write(pw_in)
    with open("wstat.in", "w") as f:
        f.write(wstat_in)

    # Useful summary
    pristine = 8 * (supercell_n**3)
    print("=== NV diamond input generation complete ===")
    print(f"Supercell: {supercell_n}x{supercell_n}x{supercell_n} conventional")
    print(f"Pristine atoms: {pristine}")
    print(f"Defective atoms (nat): {nat}")
    print(f"Charge state: {'NV-' if tot_charge == -1.0 else 'NV0'}")
    print(f"Chosen C->N index (in pristine indexing): {iN}")
    print(f"Chosen vacancy index (in pristine indexing): {ivac}")
    print("Wrote: pw.in, wstat.in")


if __name__ == "__main__":
    # ---- EDIT THESE ----
    supercell_n = 4  # 3 or 4
    charge_state = "minus"  # "minus" (NV-) or "neutral" (NV0)
    n_pdep_eigen = 400  # start point; converge this
    # --------------------

    build_nv_diamond(
        supercell_n=supercell_n, charge_state=charge_state, n_pdep_eigen=n_pdep_eigen
    )
