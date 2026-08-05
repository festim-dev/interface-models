# Scripts

Simulation scripts backing the paper. Copied from the companion repo
`festim-dev/interface-models` (local clone at `~/interface-models`), whose
history no longer exists on the GitHub remote, the remote `main` was replaced
by the Overleaf import of this paper.

## Environment

```bash
conda env create -f environment.yml
conda activate interface-models-env
pip install --no-deps git+https://github.com/festim-dev/FESTIM@interface-flux
```

The env `interface-models-env` already exists on this machine with festim
`2.2rc2.dev22+gde0e82d95` (tip of `interface-flux`) and dolfinx `0.10.0`. Note
that `interface-flux` now also needs `io4dolfinx`, which `--no-deps` will not
pull in.

## Contents

| Directory | Paper section | What it does |
|---|---|---|
| `1-first-order-interface/` | Models §2.3, verification §4.1–4.2 | `InterfaceFlux`-style first-order exchange on a 1D two-slab problem. `analytical_solution.py` is the closed-form steady state (interfacial resistance in series with `L/D`); `example_usage.py` runs one case and overlays it; `parametric_study.py` sweeps Da over `logspace(-2, 2)` transient + steady and writes the two PDFs; `verification_lte_limit.py` is the §4.1 verification, an eight-decade steady sweep (`logspace(-2, 6)`) checking the solution against `analytical_solution.py` and measuring the order of convergence to LTE in `1/Da`, writing `verification_lte_limit.csv` and `.pdf` (the latter is Fig. 1 of the paper). |
| `2-higher-order-reactions/` | Models §2.4, §2.6; verification §4.3 | `example_reaction_interface.py`: `2A ⇌ B` across the interface, i.e. recombination into a molecular carrier. `isotopologue_equilibrium.py`: zero-flux test of the H/T channel set, measuring `K_exch = c_HT²/(c_H2 c_T2)`. |
| `reference/` |, | Read-only snapshot of `src/festim/subdomain/interface.py` from `interface-flux` at the commit in `FESTIM_COMMIT.txt`. Reference only; edit FESTIM itself, not this copy. |

Directory names were kebab-cased on import (`1-First order kinetic interface`
→ `1-first-order-interface`, etc.) so paths are usable from scripts and LaTeX.

**Interface trapping is out of scope for this paper.** The `3-interface
trapping/` directory of the companion repo (Be/BeO, W/Cu, raw DOLFINx) was
deliberately not imported; it belongs to the codimension-1 companion paper.
