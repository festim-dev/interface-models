# When chemical potential continuity fails: kinetic interface models for hydrogen isotope transport

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22112697.svg)](https://doi.org/10.5281/zenodo.22112697)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Manuscript source and all simulation code for the paper by Remi
Delaporte-Mathurin and James Dark (MIT).

## The paper

Macroscopic hydrogen transport codes model material interfaces with local
thermodynamic equilibrium (LTE), imposing continuity of chemical potential as a
per-species algebraic constraint. That condition assumes fast interfacial
equilibration, a single exchange pathway between the two sides, and a carrier
species known in advance on each side. The paper argues that the last two are
the more consequential failures, and that neither survives at a metal/molten
salt interface.

We replace the constraint with reversible reaction channels on the interface
obeying mass action, detailed balance fixing each ratio of rate constants from
the same thermodynamic data that parameterises LTE. LTE is recovered as the
fast-kinetics limit of a *single* channel, so the framework generalises LTE and
does not compete with it. A Damköhler number delimits validity within a
channel, a branching ratio between channels. The framework is implemented in
[FESTIM](https://github.com/festim-dev/FESTIM) and applied to a representative
nickel/FLiBe system, where hydrogen partitions kinetically between molecular
and fluoride carriers and the apparent interfacial law drifts between Sieverts
and Henry with loading and salt redox state.

## Repository layout

| Path | What it is |
|---|---|
| `main.tex` | The whole manuscript (`elsarticle`, single file). |
| `references.bib` | Bibliography. Read the verification notice near the top before editing. |
| `scripts/` | Every simulation, verification and figure script in the paper. |

Every figure and every number quoted in the paper is produced by a script in
`scripts/`, and each script writes its PDF next to itself; `main.tex` sets
`\graphicspath` to those directories, so figures are included by bare filename.
`scripts/README.md` maps each directory and each figure to the section it
backs, and states the Damköhler convention shared by the scripts and the paper.

Broadly:

- `0-schematics/`: illustrative figures, no solve.
- `1-first-order-interface/`: Model 1, first-order exchange on a 1D two-slab
  problem, with its analytical steady state and the LTE-limit verification.
- `2-higher-order-reactions/`: Model 2, recombination into a molecular carrier,
  the isotopologue channel set, and the analytical `(Da, B)` regime map.
- `3-competing-channels/`: Model 3, recombination and fluorination competing on
  one interface, and the apparent-exponent verification.
- `4-hyperion/`: the dimensional Ni/FLiBe application, all parameters pinned in
  `parameters.py`.
- `reference/`: read-only snapshot of the FESTIM interface module the paper
  describes.

## Running the scripts

```bash
cd scripts
conda env create -f environment.yml
conda activate interface-models-env
pip install --no-deps git+https://github.com/festim-dev/FESTIM@interface-flux
```

The interface classes used here (`InterfaceFlux`, `InterfaceReaction`) live on
the `interface-flux` branch of FESTIM and are not in a release yet, hence the
extra install on top of the conda environment. `--no-deps` will not pull in
`io4dolfinx`, which that branch also needs, so install it separately if the
import fails.

Each script is standalone and self-describing: run it directly and it produces
its figure and, where relevant, the CSV of the sweep it performs.

```bash
python 1-first-order-interface/verification_lte_limit.py
```

The committed CSVs are the outputs of those runs, so a rerun can be diffed
against them.

## Building the manuscript

```bash
latexmk -pdf main.tex
```

The class options are `[preprint,11pt]` for drafting and `[5p]` for submission.

Building requires a local LaTeX installation (`elsarticle`, `siunitx` v3); the
manuscript is otherwise built on Overleaf, which this repository is synced to.

## Licence

Released under the MIT licence, see [LICENSE](LICENSE).
