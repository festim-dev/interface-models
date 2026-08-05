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
| `0-schematics/` | Models §2.2--2.5 | Illustrative figures, no solve involved. `interface_models_schematic.py` draws the LTE closure and Models 1--3 on one two-subdomain geometry. |
| `1-first-order-interface/` | Models §2.3, verification §4.1–4.2 | `InterfaceFlux`-style first-order exchange on a 1D two-slab problem. `analytical_solution.py` is the closed-form steady state (interfacial resistance in series with `L/D`); `example_usage.py` runs one case and overlays it; `parametric_study.py` sweeps Da over `logspace(-2, 2)` transient + steady; `two_slab_setup.py` draws the problem and its series-resistance reading, analytically, for the §4.1 setup figure; `verification_lte_limit.py` is the §4.1 verification, an eight-decade steady sweep (`logspace(-2, 6)`) checking the solution against `analytical_solution.py` and measuring the order of convergence to LTE in `1/Da`. |
| `2-higher-order-reactions/` | Models §2.4, §2.6, §2.7; verification §4.1, §4.3 | `example_reaction_interface.py`: `2A ⇌ B` across the interface, i.e. recombination into a molecular carrier. `analytical_solution_model2.py` is the closed-form steady state of that channel on the two-slab problem (the root of a quadratic in the rate, not a sum of resistances) together with its Sieverts/Henry limit; `verification_model2_lte_limit.py` is the §4.1 metal/liquid verification, the same eight-decade sweep (`logspace(-2, 6)`) as `1-first-order-interface/verification_lte_limit.py` but on a second-order channel, so the sweep is controlled by `2 k_plus c_0 L_m/D_m` and the local Da is reported as a diagnostic. `isotopologue_equilibrium.py`: zero-flux test of the H/T channel set, measuring `K_exch = c_HT²/(c_H2 c_T2)`. `regime_map.py`: the `(Da, B)` regime map of §2.7, analytical, no solve involved. |
| `reference/` |, | Read-only snapshot of `src/festim/subdomain/interface.py` from `interface-flux` at the commit in `FESTIM_COMMIT.txt`. Reference only; edit FESTIM itself, not this copy. |

## Figures

`main.tex` sets `\graphicspath` to the three directories above, so figures are
included by bare filename and stay next to the script that writes them.

| PDF | Paper | Written by |
|---|---|---|
| `interface_models_schematic.pdf` | §2.3, the four closures | `interface_models_schematic.py` |
| `two_slab_setup.pdf` | §4.1, the problem and its resistances | `two_slab_setup.py` |
| `regime_map.pdf` | §2.7, the `(Da, B)` map | `regime_map.py` |
| `parametric_study_damkohler.pdf` | §4.1, the transient | `parametric_study.py` |
| `verification_lte_limit.pdf` | §4.1, Model 1 steady + convergence | `verification_lte_limit.py` |
| `verification_model2_lte_limit.pdf` | §4.1, Model 2 on metal/liquid | `verification_model2_lte_limit.py` |
| `isotopologue_equilibrium.pdf` | §4.3, the zero-flux test | `isotopologue_equilibrium.py` |

`parametric_study_damkohler_steady.pdf` is not used by the paper: its content is
the upper panel of `verification_lte_limit.pdf`, over a wider Da range and from
the same solves as the convergence panel below it.

Directory names were kebab-cased on import (`1-First order kinetic interface`
→ `1-first-order-interface`, etc.) so paths are usable from scripts and LaTeX.

**Interface trapping is out of scope for this paper.** The `3-interface
trapping/` directory of the companion repo (Be/BeO, W/Cu, raw DOLFINx) was
deliberately not imported; it belongs to the codimension-1 companion paper.
