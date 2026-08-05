# CLAUDE.md

Context for working on the *advanced interface models in FESTIM* paper.

## What this repo is

A single-paper LaTeX repo (Overleaf-synced; `origin` =
`github.com/festim-dev/interface-models.git`, only commit so far: "Initial
Overleaf Import").

- `main.tex`: the whole manuscript. Currently a **structured outline**: every
  section body is an `itemize` block of notes, not prose. The job is to convert
  those bullets into prose section by section, keeping the equations and the
  argument order already fixed.
- `references.bib`: bibliography. **Read the plain-text VERIFICATION NOTICE
  block at line ~80 before touching it.** Several entries were assembled from
  web search, not exported from Zotero; placeholders deliberately carry
  `author = {TODO, Replace}` so they fail loudly.
- `scripts/`: all simulation code, imported from the companion repo. See
  `scripts/README.md` for the layout and which paper section each directory
  backs.

Class: `elsarticle`, `[preprint,11pt]` for drafting, `[5p]` for submission.
Target journal is set to **Fusion Engineering and Design** (see the TODO(journal)
comment for the alternatives considered).

There is **no LaTeX toolchain installed** (`latexmk`, `pdflatex`, `tectonic` all
absent). Don't try to compile; check syntax by reading. If a build is needed,
say so rather than installing TeX Live unprompted.

## The argument of the paper (so edits stay on-message)

Macroscopic H-transport codes close material interfaces with local thermodynamic
equilibrium (LTE), continuity of chemical potential, imposed as an algebraic
per-species constraint. Three assumptions hide in that closure:

- (H1) interfacial equilibration is fast vs bulk transport;
- (H2) a **single** exchange pathway connects the two sides;
- (H3) the carrier species on each side is known a priori.

The literature scrutinises (H1). **The paper's claim is that (H2) and (H3) are
the more consequential failures**, and that both are unavoidable at a
metal/molten-salt interface. Replace the algebraic constraint with reversible
mass-action reaction channels on the interface; detailed balance fixes
`k+/k-` from thermodynamics, so LTE is recovered exactly as the fast-kinetics
limit of a *single* channel (Sieverts/Sieverts from Model 1, Sieverts/Henry
from Model 2). Two dimensionless groups delimit validity: a Damköhler number
**within** a channel, a branching ratio `B` **between** channels. Application:
HYPERION (Ni / molten FLiBe / cover gas), single-isotope and H/T.

Recurring rhetorical moves to preserve when writing prose:

- "strict generalisation, not a competing model", because detailed balance ties
  the rate ratio to the same thermodynamics LTE uses.
- "steady state ≠ equilibrium", a permeating interface carries net flux, so
  it is never at equilibrium; large Da is **necessary but not sufficient**.
- with H and T, 2 metal-side species feed 5 salt-side carriers through a
  **bilinear, non-diagonal** map ⇒ per-species LTE is not inaccurate, it is
  **ill-posed**.
- a solubility constant has **no law-independent units** (Sec. `sec:units`);
  the measured pressure exponent is a readout of the branching ratio, not a salt
  property. This is the paper's falsifiable prediction: a redox sweep at fixed
  T should move the exponent continuously between 1 and 2.
- honest caveats are kept, not softened (the Gibbs-minimisation objection in
  `sec:model3`, the parameter-provenance caveat in `sec:limitations`). They
  pre-empt referees; don't delete them when converting bullets to prose.

Interface **trapping** / interfacial inventory is explicitly out of scope, it
is the companion paper (`festim_codim1_inprep`), and that is where the Be/BeO
and W/Cu cases belong.

## Related repos on this machine

| Path | What it is |
|---|---|
| `scripts/` (here) | All simulation code, imported from the companion repo. Work here, not in `~/interface-models`. |
| `/home/remidm/interface-models` | Original companion code repo. **Its history exists only locally**: the GitHub remote `main` was replaced by the Overleaf import of this paper. Kept as the provenance of `scripts/`; superseded by it. |
| `/home/remidm/FESTIM` | FESTIM source. The relevant branch is `interface-flux` (remote `origin` = festim-dev). Currently checked out on another branch, so use `git show origin/interface-flux:<path>` rather than switching. |
| `/home/remidm/hyperion` | HYPERION modelling scripts (dry runs, FLiBe permeability, meshes). Source for Sec. 5 setup numbers. |

`scripts/` layout (details in `scripts/README.md`):

- `1-first-order-interface/`: `example_usage.py` (1D two-slab,
  `InterfaceReaction` with `reactants=[A]`, `products=[A]`),
  `analytical_solution.py` (closed-form interfacial concentrations and flux with
  interfacial resistance in series with `L/D`, this is the analytical backbone
  of `sec:analytical` / `app:analytical`), `parametric_study.py` (Da sweep over
  `logspace(-2, 2)`, transient + steady; produces
  `parametric_study_damkohler.pdf` and `..._steady.pdf`, these are the
  Sec. `sec:lte_limit` figures).
- `2-higher-order-reactions/example_reaction_interface.py`: `2A ⇌ B` across the
  interface, i.e. **Model 2**. Byte-identical to the file of the same name on
  the `interface-flux` branch.
- `reference/interface.py`: read-only snapshot of the branch's
  `src/festim/subdomain/interface.py` at the commit in `FESTIM_COMMIT.txt`.

**Interface trapping is out of scope.** The companion repo's `3-interface
trapping/` (Be/BeO, W/Cu, raw DOLFINx) was deliberately not imported and its
results must not be cited here, they belong to the codimension-1 companion
paper `festim_codim1_inprep`.

Environment: conda env `interface-models-env`
(`/home/remidm/miniconda3/envs/interface-models-env/bin/python`), which already
has festim `2.2rc2.dev22+gde0e82d95` (= tip of `interface-flux`) and dolfinx
`0.10.0` installed. Plotting uses `morethemes` (`mt.set_theme("urban")`) and the
palette `#1a4848, #f7b000, #f46036, #c9f2c7, #aceca1`, top/right spines
removed. Match that style for any new figure.

## The FESTIM API this paper describes

On `interface-flux`, `src/festim/subdomain/interface.py` exports three interface
classes (all re-exported from `festim/__init__.py`), used via
`F.HydrogenTransportProblemDiscontinuous(...)` and its `interfaces` list.

**`Interface`** (pre-existing, = LTE): `penalty` or `nitsche` method. Enforces
`u_0/K_0 = u_1/K_1`, or squares the Sieverts side when the two solubility laws
differ. This is the baseline the kinetic models must reproduce at large Da.

**`InterfaceFlux(id, subdomains, k_plus, k_minus)`** = **Model 1**. Adds
`R = k+·u_0 − k−·u_1` to the residual of subdomain 0 and `−R` to subdomain 1,
per species. Species must exist in both subdomains.

**`InterfaceReaction(id, subdomains, k_plus, k_minus, reactants, products)`** =
**Models 2/3**, the general channel. Key semantics:

- `reaction_term()` = `k+ · Π(reactant concentrations) − k− · Π(product
  concentrations)`; reactants are read on `subdomains[0]`, products on
  `subdomains[1]`.
- **Stoichiometry is encoded by repetition in the list.** `reactants=[A]*2`
  gives `k+·c_A²` *and* adds `R` to A's residual twice, i.e. a factor 2, which
  is exactly Eqs. (`eq:bc_metal`)/(`eq:bc_salt`): flux of atoms out of the metal
  is `2w`, production of the molecular carrier is `w`.
- Reactants and products may live on **different** subdomains and be different
  species (`A` on vol1, `B` on vol2), this is what makes multi-carrier and
  isotopologue channels expressible.
- Multiple channels on one interface means multiple `InterfaceReaction` objects
  with the same `id`, whose residual contributions add. That is Model 3, and it
  currently **crashes upstream**: see the known-limitation note below.
- No interfacial degrees of freedom are introduced; everything is algebraic in
  the trace values. Consistent with the claim in `sec:weak`.

**Mass-action convention (verified in code and numerically).** `Reaction`,
`SurfaceReactionBC` and `InterfaceReaction` all compute the rate as the plain
product of reactant concentrations with **no combinatorial or statistical
prefactor**. `SurfaceReactionBC`'s docstring says so explicitly: "In the special
case where A=B, then the flux of particle entering the surface is 2*K". So the
statistical degeneracy of a mixed isotopologue pair lives in the *value* of
`k_HT+` (= `2 k_HH+` in the mass-independent limit), never in the rate law.
`scripts/2-higher-order-reactions/isotopologue_equilibrium.py` demonstrates
this: equal forward constants give `K_exch = 1`, and `k_HT+ = 2 k_HH+` gives
`K_exch = 4` (the correct classical value, H2:HT:T2 = 1:2:1). All cases
reproduce to 6 significant figures.

The degeneracy **cannot** be encoded structurally. Declaring the mixed channel
twice (`reactants=[H, T]` and `reactants=[T, H]`) applies the whole rate `R` to
each residual twice, doubling forward and reverse equally: it rescales the
channel in time and leaves `K_exch = 1`. Same for any other repetition trick,
because `R` carries both directions. The degeneracy is forward-only (two ways
to pick an H,T pair; one way to dissociate HT), so it must go into a constant:
`k_HT+ = 2 k_HH+`, or equivalently `k_HT- = k_HH-/2`. Worth one sentence of
implementation guidance in `sec:weak` or `sec:params` when Sec. 3 is written;
users will reach for the duplication trick. Any equation
written in the paper must match this convention, since the paper claims to
describe the implementation.

**Several channels on one interface `id` now work.** This used to crash in
`create_formulation` (festim-dev/FESTIM#1222); fixed on `interface-flux` by
deduplicating the integration data, with a system test in
`test/system_tests/test_interface_reactions.py`. Model 3 and the isotopologue
set assemble against stock `HydrogenTransportProblemDiscontinuous`, so no
subclass or workaround is needed.

**Gaps between the manuscript and what the branch actually implements**, check
against the code before claiming a capability in the paper:

- `k_plus`/`k_minus` are plain numbers. There is **no Arrhenius / temperature
  dependence** and no `a_F`-style activity prefactor; `InterfaceReaction.
  get_formulation` ignores its `temperature` and `species` arguments entirely.
  Model 3's `a_F` currently has to be folded into `k_f+` by hand.
- **Detailed balance is not enforced by the code.** The `k+/k−` ratio is the
  user's responsibility. If the paper claims thermodynamic consistency is
  guaranteed, that is a statement about the framework, not about the current
  implementation, phrase accordingly.
- No `E_k` / activation-energy plumbing, no positivity safeguard, no
  interface-residual scaling. The `TODO`s in `sec:numerics` about stiffness and
  positivity are genuinely open, not just unwritten.

## Working rules

- **Never invent a citation.** If a claim needs a reference that is not in
  `references.bib`, insert a `\TODO{}` or `\VERIFY{}` naming what is needed.
  The three annotation macros are `\TODO{}` (red), `\VERIFY{}` (orange),
  `\KEY{}` (blue); all three are deleted before submission.
- **Never silently upgrade a placeholder bib entry.** Entries with
  `author = {TODO, Replace}` (`hyperion_facility`, `hyperion_multidim`,
  `libra_2022`, `mit_thesis_flibe`) must keep failing loudly until re-exported
  from Zotero. Same for the `note` fields recording what is unverified.
- **Do not fabricate numbers.** No parameter table rows, HYPERION dimensions,
  fitted `k_r`/`k_f`, or simulation results unless they come from
  `/home/remidm/interface-models`, `/home/remidm/hyperion`, or a source the user
  supplies. Numerical results quoted in the paper should be traceable to a
  script in the companion repo.
- Notation macros already defined: `\FESTIM \TMAP \MHIMS \HYPERION \Da \kB
  \half \flibe \cm{} \cs{} \Bra`. `\cm{i}` = metal-side (atomic), `\cs{α}` =
  salt-side (molecular/fluoride), `\Bra` = branching ratio. Reuse them; don't
  hand-roll new notation for the same quantity.
- When converting a bullet block to prose, keep every equation and its label;
  labels are cross-referenced throughout (`eq:lte_ss`, `eq:lte_sh`,
  `eq:mass_action`, `eq:model1_flux`, `eq:model2_rate`, `eq:model3_rate`,
  `eq:branching`, `eq:damkohler`, `eq:iso_*`).
- Prose style: British spelling ("realised", "generalisation"), first person
  plural, no bullet lists in the final text except the explicit decision
  procedure in `sec:when_lte`.
- **Never use em dashes** (`---`, `--`, or the character itself). Use a comma,
  a colon, a semicolon, parentheses, or two sentences. This applies to
  everything written here: prose, headings, captions, comments.
- **No AI tells.** Specifically avoid: "it is worth noting", "it is important
  to note", "delve", "leverage", "showcase", "underscore", "crucially",
  "notably", "robust" as filler, "a testament to"; the "not only X but also Y"
  and "it is not X; it is Y" constructions when they are decoration rather than
  the actual argument; three-item lists padded to three for rhythm; a
  summarising sentence at the end of a paragraph that restates the paragraph;
  hedges stacked together ("may potentially suggest"). Write the claim once,
  in the plainest order, and move on. Some emphatic constructions in the
  outline are the author's own voice (e.g. "LTE is not inaccurate; it is
  ill-posed") and should be preserved as written.
- The single most valuable planned figure, per the outline, is the two-axis
  `(Da, B)` regime map in `sec:damkohler`.
