# AI-Guided Multi-Stage Screening Pipeline (representative excerpt)

This directory documents the multi-stage computational screening cascade
used to nominate candidate compounds against HSPH1 and ST6GALNAC1
(manuscript Figure 9). The full pipeline integrates four stages; a
representative command sequence for each stage is given below.

## Stage 1 -- drug-target interaction prediction (GraphBAN)

A graph-based bilinear attention network (GraphBAN) jointly encodes drug
molecular graphs and protein sequences to predict binding probability.
Two independently trained models (BioSNAP and KIBA interaction datasets)
were applied in parallel to reduce single-model bias.

```bash
# representative inference call (PyTorch checkpoint)
python predict_graphban.py \
  --checkpoint models/graphban_biosnap.pt \
  --drugs screening_input/compound_graphs.pkl \
  --targets screening_input/target_sequences.fa \
  --output scoring/graphban_biosnap_scores.csv
```

Candidate retention: compounds scoring above the 0.9 percentile for
both targets in both models.

## Stage 2 -- ADMET filtering (ADMET-AI)

ADMET-AI predicts absorption, distribution, metabolism, excretion and
toxicity properties from SMILES.

```bash
admet-ai predict --smiles screening_input/candidates.smi \
  --output scoring/admet_ai_results.csv
```

Candidate retention: compounds passing the ADMET-AI reference criteria
for bioavailability and toxicity.

## Stage 3 -- structure-based docking

Structure-based docking was performed against AlphaFold-predicted
domains: the ATPase/nucleotide-binding domain of HSPH1 and the
catalytic domain of ST6GALNAC1. Binding pockets were transferred from
homologous experimental structures (PDB 3FE1 and 6APL) where bound
ligands are absent in the predicted models.

```bash
# receptor and grid preparation (representative)
prepare_receptor.py -r hspH1_nbd.pdb -o docking/hsph1_receptor.pdbqt \
  --pocket-from-pdb 3FE1
# docking (representative; exhaustiveness 32)
vina --receptor docking/hsph1_receptor.pdbqt \
     --ligand screening_input/candidates.pdbqt \
     --config docking/hsph1_box.conf --exhaustiveness 32 \
     --out docking/hsph1_docked.pdbqt
```

Candidate retention: top-ranked compounds by docking affinity per target
(e.g., ZINC000084931393 for HSPH1, ZINC000013548644 for ST6GALNAC1).

## Stage 4 -- molecular dynamics (GROMACS)

Production simulations of the two protein-ligand complexes (100 ns each)
were performed with GROMACS 2025.4 (CHARMM36m force field, TIP3P water,
physiological NaCl).

```bash
# representative protocol
gmx pdb2gmx -f complex.pdb -o complex_processed.gro -ff charmm36m -water tip3p
gmx editconf -f complex_processed.gro -o box.gro -c -d 1.2 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o neutral.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
gmx grompp -f minim.mdp -c neutral.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
gmx grompp -f md_production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_100ns.tpr
gmx mdrun -deffnm md_100ns -nb gpu
```

Trajectory metrics (RMSD, protein-ligand contacts, minimum distance,
hydrogen bonds) are computed by `../scripts/09_fig9_md_analysis.py`.

> Note: this directory contains a representative excerpt of the
> screening workflow; input libraries, model checkpoints and
> intermediate scoring files are not redistributed.
