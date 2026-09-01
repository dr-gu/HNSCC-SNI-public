#!/usr/bin/env python3
# ============================================================
# 09_fig9_md_analysis.py
# Figure 9 -- molecular dynamics trajectory metrics
#
# Computes the stability metrics reported for the two docked
# complexes (HSPH1 and ST6GALNAC1 with their top candidates):
#   - backbone RMSD (reference: first frame)
#   - protein-ligand contact counts (< 4.5 A)
#   - minimum protein-ligand distance
#   - hydrogen-bond counts
#
# Trajectories were produced with GROMACS (see screening/README.md
# for the simulation protocol). Representative excerpt only.
# ============================================================

import sys

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms, contacts, hydrogenbonds

CONTACT_CUTOFF = 4.5  # Angstrom


def analyse(psf_path: str, trajectory_path: str, ligand_resname: str) -> dict:
    u = mda.Universe(psf_path, trajectory_path)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(f"resname {ligand_resname}")

    # backbone RMSD against the first frame
    rmsd = rms.RMSD(protein, protein, select="backbone").run().results.rmsd

    # protein-ligand contacts and minimum distance along the trajectory
    cg = contacts.Contacts(
        u, select=(f"protein and around {CONTACT_CUTOFF} resname {ligand_resname}",
                   f"resname {ligand_resname}"), radius=CONTACT_CUTOFF).run()
    n_contacts = cg.n_contacts
    min_dist = cg.timeseries.min(axis=1)

    # hydrogen bonds (donor-acceptor distance < 3.0 A, angle > 150 deg)
    hbonds = hydrogenbonds.HydrogenBondAnalysis(
        u, selection1=f"protein and around {CONTACT_CUTOFF} resname {ligand_resname}",
        selection2=f"resname {ligand_resname}",
        distance=3.0, angle=150.0).run()
    n_hbonds = hbonds.count_by_time()

    return {
        "rmsd_mean_nm": float(np.mean(rmsd) / 10.0),   # A -> nm
        "contacts_mean": float(np.mean(n_contacts)),
        "min_dist_max_nm": float(np.max(min_dist) / 10.0),
        "hbonds_mean_per_frame": float(np.mean(n_hbonds)),
        "trajectory_frames": len(rmsd),
    }


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: 09_fig9_md_analysis.py <topology> <trajectory> <ligand_resname>")
        return 1
    metrics = analyse(sys.argv[1], sys.argv[2], sys.argv[3])
    for key, value in metrics.items():
        print(f"{key}: {value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
