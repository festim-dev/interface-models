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
per-species constraint. Three assumptions hide in that condition:

- (A1) interfacial equilibration is fast vs bulk transport;
- (A2) a **single** exchange pathway connects the two sides;
- (A3) the carrier species on each side is known a priori.

The literature scrutinises (A1). **The paper's claim is that (A2) and (A3) are
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
  `analytical_solution.py` (analytical interfacial concentrations and flux with
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

### Sources and numbers

- **Never invent a citation.** If a claim needs a reference that is not in
  `references.bib`, insert a `\TODO{}` or `\VERIFY{}` naming what is needed.
  The three annotation macros are `\TODO{}` (red), `\VERIFY{}` (orange),
  `\KEY{}` (blue); all three are deleted before submission.
- **Never silently upgrade a placeholder bib entry.** Entries with
  `author = {TODO, Replace}` must keep failing loudly until re-exported from
  Zotero. Same for the `note` fields recording what is unverified.
- **Do not fabricate numbers.** No parameter table rows, HYPERION dimensions,
  fitted `k_r`/`k_f`, or simulation results unless they come from the companion
  repo or a source the user supplies. Every number quoted in the paper should be
  traceable to a script.
- **Sec. 5 is a demonstration, not a validation.** Nothing there is fitted or
  compared with data, and the paper's own argument is that the salt-side
  parameterisation does not currently exist. Do not quote its numbers as
  findings about FLiBe, in the abstract, the introduction, or anywhere else.
  Directional statements ("LTE underestimates the steady flux, by an amount the
  redox state sets") are safe; specific factors presented as measurements are
  not.

### LaTeX conventions

- **One sentence per line.** Never wrap a sentence across lines and never put
  two sentences on one line. Both authors edit in Overleaf simultaneously, so
  this keeps diffs sentence-sized and lets a line in the PDF be found by search
  in the source.
- Notation macros already defined: `\FESTIM \TMAP \MHIMS \HYPERION \Da \kB
  \half \flibe \cm{} \cs{} \Bra`. Reuse them; do not hand-roll new notation for
  the same quantity.
- **Units go through `siunitx`** (v3 macros): `\qty{2.032}{\milli\metre}`,
  `\unit{\metre\per\second}`, `\num{1.31e5}`. Never hand-roll
  `$1.31\times10^{5}$~Pa`. Spell units out in full inside `\qty`/`\unit`
  (`\centi\metre`, never `\cm`): siunitx ships a short `\cm` that collides with
  the concentration macro. Clean powers of ten that are part of a maths
  relation (`$\Da = 10^{2}$`) stay as maths.
- Keep every equation label when rewriting; labels are cross-referenced
  throughout.
- **Figure and table references go at the end of the sentence** as
  `(see Fig.~\ref{...})`, or as the subject of their own sentence
  (`Figure~\ref{...} shows ...`). Never mid-clause. This follows the authors'
  own papers.
- **Never use em dashes** (`---`, `--`, or the character itself). Use a comma,
  a colon, a semicolon, parentheses, or two sentences. Applies to prose,
  headings, captions and comments.

### Prose style

- British spelling ("realised", "generalisation").
- **First person is for our moves and our findings; everything else is
  impersonal.** "We argue that (A2) and (A3) are the more consequential
  failures" keeps its *we*. "Reported FLiBe transport properties span orders of
  magnitude" does not take one. Roughly once per 200 words overall.
- **One term, one meaning.** A term that carries two senses in one paper is an
  error, however clear each sense is locally. See the terminology registry
  below; add to it whenever a new term is introduced.
- **Sentence cap: 25 words** for descriptive prose. Longer is allowed only when
  splitting would break a logical relation that has to be held in one breath.
- **Paragraph cap: six sentences.**
- **Active voice by default.** "Evaluate one Damköhler number per side", not
  "one Damköhler number should be evaluated per side".
- No bullet lists in the final text except the explicit decision procedure in
  `sec:when_lte`.
- **Cut "btw points".** A sentence that gives a further example, or shows that
  an idea also appears in another field, and which no later sentence depends
  on, comes out. Test: delete it, and check whether any later sentence breaks
  or whether the reader loses something needed to use or judge the model. These
  cluster immediately after a claim, as a second supporting citation that adds
  a field rather than a fact. Defensive writing that answers an objection a
  referee will raise is load-bearing and stays.
- **Core derivations belong in the main text.** Appendices carry supporting
  proof for things not central to the argument. If a result is quoted in three
  places and named as the paper's falsifiable content, derive it where the
  reader meets it, however short the algebra.
- **No AI tells.** Avoid: "it is worth noting", "it is important to note",
  "delve", "leverage", "showcase", "underscore", "crucially", "notably",
  "robust" as filler, "a testament to"; the "not only X but also Y" and "it is
  not X; it is Y" constructions when they are decoration rather than the actual
  argument; three-item lists padded to three for rhythm; enumerators
  ("First, ... Second, ... Third, ...") where the preceding sentence has
  already announced the count; a summarising sentence that restates its own
  paragraph; stacked hedges ("may potentially suggest"). Some emphatic
  constructions are the author's own voice (e.g. "LTE is not inaccurate; it is
  ill-posed") and are preserved as written.
- **Three tics to keep counting.** Re-grep before declaring a section done.
  - `rather than`: the authors' own papers use it essentially never. Prefer
    restructuring, "instead of", "and not", or a second sentence, and vary the
    replacement so "and not" does not become the next tic.
  - the trailing gloss `, which is <restatement>`: drop the "which is" and leave
    the bare appositive, or split the sentence.
  - `This is the <noun> that/used in ...` as a cross-reference. Say "We use this
    solution in Sec. 4.1".
- **Banned words: "closure" and "closed form".** For the algebraic LTE
  constraint write "condition", "interface law", "constraint", or "model". For
  the analytical steady states write "analytical solution". An interface is
  "treated" or "described" a given way, it is not "closed" by it.
- **Minimise Latin.** "a priori" and "a fortiori" were removed from the
  introduction. Prefer "in advance", "known beforehand", and drop
  "a fortiori" entirely, since the sentence structure usually carries it.
- **Voice benchmark.** Delaporte-Mathurin et al., arXiv:2603.25751 (PathSim) is
  the closest sample of the first author's own prose: plainer and more
  expository than this draft, mean sentence ~23 words, plain connectives,
  almost no compressed epigrams. Match that register.

## Terminology registry

One term, one meaning. Each row records the sense the paper uses and what was
banned or renamed to get there.

| Term | The one meaning | Notes |
|---|---|---|
| operating point | the temperature, pressure and geometry a case is run at | Was also used for the balance point of competing half-reactions (Sec. 2.5) and for the state a nonlinear rate is linearised about (Secs. 2.7, 4.1, App. C). Both other senses removed: say "the corrosion rate is set by the competition between them", and "linearising about a chosen interfacial loading". |
| Γ | the interface | Never a flux. Fluxes are φ; the in/out pair in `sec:lte` is φ_in, φ_out. |
| λ | the atomic jump distance (Sec. 2.3) | Also the stoichiometric coefficient of the non-hydrogenic constituents in Eq. (3). Resolve before submission: either rename one, or fold the activities into a per-channel factor and drop λ from Eq. (3). |
| mol vs particles | Sec. 2 counts moles, Sec. 5 and the scripts count particles | Second-order constant is m⁴ mol⁻¹ s⁻¹ in Sec. 2 and m⁴ s⁻¹ in Fig. 10 and `parameters.py`. State the switch where it happens. |
| (A1), (A2), (A3) | the three hidden assumptions | Not (H1)–(H3). Check Secs. 5.1 and 6.1. |
| bilinear, non-diagonal, diagonal map | **banned** | Removed 2026-08-17: meaningless without a matrix picture the paper does not draw. Say that HT is made from one H and one T, so its concentration depends on both metal-side concentrations at once, and a condition written one species at a time has nothing to equate it to. |
| Da vs Da⋆ | Da is one-sided, Eq. (24); Da⋆ is two-sided, Eq. (B.2) | Fig. 2 plots Da⋆, so thresholds quoted against that figure say Da⋆. |
| BeF₄²⁻ | the tetrafluoroberyllate anion | Charge is 2−. |