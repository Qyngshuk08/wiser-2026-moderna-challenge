# Presentation Script — Quantum RNA Secondary Structure Prediction
Team Qudit Creons — WISER 2026 Moderna Challenge

Read this while screen-recording the PDF (presentation/presentation.pdf).
Target pace: ~35-45 seconds per slide, ~11-13 minutes total for 18 slides.
Speak naturally — this is written to be SAID, not read verbatim if it feels stiff.

NOTE before recording: Slide 13 (Limitations) still lists "D-Wave's n=100
QPU trigger unexplained" as an open item. That was actually resolved later
in the project — it turned out to be a constraint-density effect from the
branching-topology fix, not an unexplained size threshold. The slide text
is stale. Script below narrates the CORRECT, resolved version verbally.
Consider fixing the slide itself before recording if you have time;
otherwise this script covers you.

---

## [Title Slide]

"This is Quantum RNA Secondary Structure Prediction, submitted by team
Qudit Creons for the WISER 2026 Moderna Challenge. I worked on this
individually. Over the course of this project I built a QUBO formulation
for RNA folding using real thermodynamic energies, ran it on real D-Wave
and IBM quantum hardware, and — maybe more importantly — I want to walk
through two real bugs I found in my own work and how I corrected them,
because that process is a big part of what makes this submission solid."

---

## Slide 1 — The Problem

"RNA secondary structure prediction — figuring out how a strand of RNA
folds back on itself into its minimum free energy shape — is
combinatorially hard. The number of possible folds grows fast with
sequence length. Classical tools like ViennaRNA solve this exactly using
dynamic programming, but that doesn't translate naturally to quantum
hardware.

So the question I set out to answer wasn't just 'can a quantum computer
fold RNA' — it was 'where exactly does a quantum approach succeed, where
does it break down, and why.' That second part — the why — ended up being
most of the actual work.

What I built: a QUBO formulation using real Turner2004 thermodynamic
energies, not made-up reward constants. I solved it on D-Wave, using both
the hybrid solver and direct QPU annealing, and via QAOA on real IBM
hardware. I benchmarked everything against ViennaRNA's real MFE output,
at both small hand-picked scale and a proper statistical scale — 320
random sequences. And every single result was cross-checked against
independent ground truth. That last part caught several real bugs, which
I'll get into."

---

## Slide 2 — Formulation

"Here's the actual math. I represent RNA structure using what I call
'quartet' variables — one binary variable per stacked base pair. So a
variable represents: is pair (i,j) stacked immediately on pair (i+1,
j-1)? That's the dominant energetic contribution in real RNA folding —
adjacent stacked pairs are what actually stabilizes a helix.

Two constraints: each base can pair at most once, obviously, and no
crossing pairs — meaning no pseudoknots, for now. I do come back to that
later, I actually built a pseudoknot-permissive extension, but the base
model excludes them.

The thing I want to highlight here, because it's the actual technical
contribution: the energy weight on each quartet isn't a tunable constant
I picked to make the optimizer behave nicely. It's the real ΔG stacking
energy, pulled directly out of ViennaRNA's own Turner2004 parameter
tables. Prior quantum annealing papers on RNA folding — Fox et al., 
Zaborniak et al. — use heuristic reward constants instead. Mine uses real
thermodynamics. You can see on the right — a 9-nucleotide toy example,
GGGAAACCC — the total I compute, negative 1.20 kilocalories per mole,
matches ViennaRNA's real answer exactly."

---

## Slide 3 — Validation Methodology

"This slide is really about how I approached the whole project, not just
one component. Every layer — brute force, dynamic programming, the QUBO
itself, real hardware — got checked against an independent source of
truth before I trusted it. Not just checked once and moved on — checked
every time something changed.

That discipline is the reason this whole project is defensible. It's
also literally what caught the two bugs I'm about to walk through,
including one that had already made it into a hardware result I'd
published."

---

## Slide 3b — Major Correction: A Foundational Indexing Bug

"This is the most important slide in the whole deck, and I want to be
direct about it rather than bury it.

Months into the project, while I was trying to extend the model to
handle bulge loops — a more complex loop type — I hit something that
didn't make sense. My bulge-aware model, which has strictly MORE options
available than my simpler hairpin-only model, was scoring a structure as
WORSE than the simpler model found for the identical sequence. That's
not possible if both models are correct — a model with more options can
never do worse than a model with fewer.

I didn't dismiss that as noise. I traced it down. The root cause: the
function that looks up real stacking energies from ViennaRNA's tables
was indexing one of the two base pairs in the wrong order — reversed
from what it should have been. This bug had been sitting in four
different files since basically the start of the project.

I confirmed it systematically — not by inspection, by testing all 36
possible canonical base-pair-type combinations against ViennaRNA's own
energy evaluator directly. Only 2 out of 36 matched with the old
indexing. All 36 matched once I reversed it.

Here's the caveat that actually matters for interpreting everything else
in this deck: the STRUCTURES my model predicted barely changed — match
rate statistics moved by about one percentage point in either direction.
But every individual ENERGY VALUE I'd reported anywhere — every DP
result, every QUBO result, every D-Wave and IBM hardware result — was
off by up to about 1.8 kilocalories per mole per stack. That's real. I
fixed it in all four files, revalidated everything from scratch, and I
reran what hardware access I had left against the fix. I could not rerun
everything — hardware access ran out. That's disclosed honestly, not
hidden."

---

## Slide 4 — Finding: Stacking-Only Match Rate Collapses at Scale

"Now into the actual results. First thing I found, and this humbled me a
bit: on a handful of small, hand-picked test sequences, my simplest
model — stacking energy only, no loop penalties — got 4 out of 5 correct
against real ViennaRNA output. That looked good. It was misleading.

Once I ran a real statistical test — 320 randomly generated sequences,
not hand-picked — the true match rate was 10 percent. And it drops to
zero percent for sequences 60 nucleotides or longer, which is exactly
the length range I later used for hardware scaling tests. The energy gap
also grows more negative as sequences get longer — meaning the model
doesn't just occasionally pick a different fold, it systematically
overestimates how stable the structure is, because it's ignoring the
real entropic cost of unpaired loop regions."

---

## Slide 5 — Finding: Real Loop Energies Fix Most of the Gap

"So the obvious next step: add real hairpin loop energy — the
thermodynamic penalty for the loop at the tip of a hairpin. I built a
dynamic program using ViennaRNA's own loop-energy evaluators, ran the
identical 320-sequence test, and match rate went from 10 percent to 53.4
percent. Real, substantial improvement.

Caveat, because I want to be honest about the process, not just the
result: while building this, I found a real bug in my own code — a 'free
multiloop fallback' that was silently zeroing out the hairpin penalty on
every single closure, which meant my early version was systematically
too optimistic. I caught it by cross-checking my computed energy against
ViennaRNA's own structure evaluator, saw they disagreed, and traced it
down rather than assuming my number was right."

---

## Slide 6 — D-Wave: Embedding Density Limits Direct QPU Use

"Moving to actual quantum hardware. On D-Wave, I ran both the hybrid
solver and direct QPU annealing. Direct QPU submission requires
embedding my problem's variables onto D-Wave's physical qubit
connectivity graph — and that's where I found the real bottleneck.

My constraint graph density stays above 47 percent even at 540 variables
— which is way denser than D-Wave's native qubit connectivity, roughly
15 to 20 connections per qubit. What that means in practice: as problem
size grows, wall time explodes — 0.28 seconds up to over 40 seconds — but
if you actually look at the QPU access time specifically, it stays flat,
around 125 to 185 milliseconds the whole time. The bottleneck isn't the
quantum hardware. It's the classical embedding search trying to fit a
dense problem onto a sparse physical graph. And when it does find an
embedding, solution quality on direct QPU was actually 39 percent worse
than the hybrid solver on the same problem — likely from chain breaks in
that dense embedding."

---

## Slide 7 — D-Wave: Hardware Rerun — Including a Confirmed Correction

"This is the second real bug I want to walk through, and it's a
significant one because it affected results I'd already published.

I reran my six standard test sequences against the hairpin-aware model on
real D-Wave hardware. Two of them — the 80 and 100 nucleotide cases — had
energies that turned out to be wrong. Not slightly wrong — the original
published values were negative 6.90 and negative 7.93; the corrected
values are negative 6.20 and negative 7.93 — sorry, negative 6.70. The
root cause was a branching-topology bug — my energy model was allowing a
base pair to be treated as closing an empty hairpin loop, when in some
selected structures that 'empty' loop actually contained another
independent, nested stem. I fixed the constraint logic, reran on real
hardware, and corrected the record. I didn't delete the old numbers —
they're still in the repository, clearly flagged as superseded, so
there's a transparent trail of what was wrong and what I changed."

---

## Slide 8 — QAOA: A Penalty-Landscape Bug — Fixed, But Not Solved

"Switching to IBM and QAOA. First problem I hit: the same penalty weight
I used for classical exact solving — sized to guarantee correctness no
matter what — turned out to be about 11 times larger than any single
favorable energy in these small problems. That's fine for a classical
solver that searches exhaustively. For QAOA, it meant the landscape was
so dominated by that penalty term that the optimizer kept collapsing to
the trivial 'select nothing' answer. Confirmed that empirically — zero
percent success across 8 restarts, and I ruled out 'just needs more
optimizer effort' as the explanation first.

Fixed it with a much tighter, QAOA-specific penalty. That eliminated the
total-collapse failure mode. But — and this is the caveat that matters —
it does not guarantee success. I ran a proper 20-seed statistical sweep,
and the real success rate is only 15 percent on one toy case and 25
percent on another, even at just 4 to 6 qubits. I want to be direct that
an earlier version of this project's documentation overstated this as
'fixed' based on a single successful run. I caught that overstatement
myself and corrected it once I had the real distribution."

---

## Slide 9 — Real IBM Hardware: Confirmed, Then a Noise Wall

"First real hardware run, on ibm_kingston, 4 qubits: exact match to the
correct answer. But only 17 percent of the 2000 shots landed on that
correct answer — compare that to 97 percent confidence on the noiseless
simulator. That gap is real hardware noise, quantified directly, not
estimated.

Second run, 6 qubits on ibm_marrakesh: the best FEASIBLE answer the
hardware found was the trivial empty structure — the correct, non-trivial
answer got completely buried by noise. Two more qubits was enough to
erase it entirely.

I also caught a bug while looking at this: the single most common answer
returned by the hardware was actually infeasible — it violated a
constraint, with an energy that was roughly twice what it should have
been, consistent with two conflicting parts of the structure both firing
at once. I fixed that by ranking all the sampled answers by frequency and
picking the highest-ranked one that's actually feasible, instead of just
trusting whatever came back most often."

---

## Slide 9b — IBM Access Restored — A Correction, Not a Confirmation

"At one point I reported IBM hardware access as permanently exhausted.
That turned out to be wrong — access came back. So I reran both toy
cases against the corrected model.

Both reruns got the wrong answer on real hardware. I want to be clear
that's not a new problem — it's exactly consistent with the 15 and 25
percent success rates I'd already measured. The combined chance of both
runs succeeding was under 4 percent. Failing both was actually the
single most likely outcome. What this really did was correct an
overstated claim in my own earlier documentation that implied the fix
guaranteed success — caught by the same validate-before-trust process
I've been describing this whole time."

---

## Slide 10 — Noise Study: Two Hardware Platforms, Different Signatures

"I also ran a real comparison of noise characteristics across both
hardware platforms, since D-Wave's annealing-based noise and IBM's
gate-based decoherence are fundamentally different physical processes.

On IBM: 3 out of 4 runs correct, across 3 different backends. Interesting
detail — the confidence range for the correct runs, 17 to 24 percent,
actually overlaps with the confidence of the one WRONG run, 15.3 percent.
So confidence alone doesn't reliably tell you whether an answer is
correct.

On D-Wave: all 3 runs found the correct answer, and the chain-break rate
was low — under 7 out of 1000 reads had any break at all. But the
read-level confidence — how often the correct answer specifically came
up — was 8 to 10 times lower than IBM's. Chain breaks alone don't explain
that gap. I want to be honest that with only 3 D-Wave runs, I can't fully
diagnose the mechanism — but the finding itself, that per-read confidence
stays low even with minimal chain breaks, is real and worth reporting."

---

## Slide 11 — Comparing Two Quantum Encodings

"One of the optional advanced tasks: compare my quartet-based encoding
against an alternative — raw base-pair variables, the encoding style used
in earlier literature. Both represent the exact same physical model,
just structured differently.

Real, consistent finding across every sequence length I tested: the raw
pair encoding needs about 3 times more qubits than my quartet encoding —
ratio between 2.94 and 3.67 — but has a somewhat less dense constraint
graph. Real trade-off between qubit efficiency and constraint density,
not a free win either way.

One more thing I caught and verified rather than dismissed: at one
sequence length, both encodings found DIFFERENT structures that were
tied at the exact same optimal energy. Confirmed independently that both
were genuinely valid, equally optimal answers — a real degenerate
optimum in the underlying energy landscape, not a bug in either
encoding."

---

## Slide 12 — Pseudoknots: A Validated Proof-of-Concept

"Pseudoknots — crossing base pairs — are excluded by the base model.
General pseudoknot folding is NP-hard, so I didn't attempt that. What I
did do: relax the no-crossing constraint and prove the solver correctly
exploits it when a crossing structure is genuinely favorable.

I hand-built a test sequence with two independent stems specifically
designed so they'd have to cross to both be used. Verified
programmatically, not by eye, that they really do cross. Without
allowing crossing: best answer uses only one stem, negative 8.70. With
crossing allowed: the solver correctly finds both stems together,
negative 17.40 — which exactly matches what I calculated by hand should
happen if both stems' energies just add together, since they don't share
any bases.

Honest caveat: I did not add real pseudoknot-specific loop-topology
energy penalties — the extra thermodynamic cost that a real crossing
structure incurs beyond simple stacking. My crossing quartets are scored
identically to nested ones. That's a real simplification, stated
directly, not hidden in a footnote."

---

## Slide 13 — Limitations

"Here's the honest state of what's not done. No multiloop support
anywhere — neither the dynamic program nor the QUBO handles multi-branch
loop structures. The QUBO still lacks bulge and internal loop energy —
I proved with a dynamic program that adding it would help, that's the
single biggest remaining gap in the actual QUBO. No real pseudoknot
energy penalties, as I just said. QAOA succeeds only 15 to 25 percent of
the time, measured directly, and I confirmed via classical simulated
annealing that this is a genuine landscape-difficulty problem, not a bug
in my formulation. IBM hardware noise dominates by around 6 qubits on
current hardware without deeper error mitigation.

One correction to what's on this slide as written: it lists D-Wave's
qubit-access trigger pattern as unexplained. I actually resolved that
later — it turned out to be caused by constraint density from the
branching-topology fix, not an unexplained size threshold. The slide
text is a bit stale; the real answer is: I found it, and it's a density
effect, confirmed with real repeated hardware runs."

---

## Slide 14 — Future Work

"Real, prioritized next steps, not a vague wishlist. Highest priority:
port bulge and internal loop energy into the actual QUBO — the dynamic
program already proves it's worth doing. Proper multiloop support next.
On the QAOA side: either problem-specific circuit depth tuning or a
different reformulation entirely, since a blanket depth increase doesn't
generalize — I tested that directly and it helped one test case while
making another one worse. Real pseudoknot energy penalties. Extending the
D-Wave noise study across more qubit counts. And rerunning historical
hardware results against the indexing fix wherever hardware access
becomes available again — that's genuinely the top priority if I get more
hardware time."

---

## [Closing Slide]

"Everything in this presentation came from checking results against
independent ground truth, and being willing to say 'that doesn't look
right' even about my own already-published numbers. I found and
corrected two real bugs during this project, one of which had already
made it into hardware results I'd shared. I think that process — not
just the final numbers — is what this submission actually demonstrates.
Thank you."

---

## Post-recording checklist

- Total runtime should land around 11-14 minutes at a natural pace.
- If you want it shorter for a strict time limit, the safest cuts are:
  trimming slide 10 (noise study) and slide 11 (encoding comparison) down
  to 2-3 sentences each — they're real findings but less central than the
  two bug corrections (3b, 7, 8, 9b), which should stay full-length since
  they're the strongest part of the submission.
- Consider fixing slide 13's stale D-Wave n=100 line in the actual .tex
  file before recording, so the script and slide agree — ask if you want
  this done before you record.
