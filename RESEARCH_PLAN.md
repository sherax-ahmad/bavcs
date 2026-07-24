# BAVCS: Bias-Aware Variant Concordance Score

## The gap (and why it's real, narrow, and unclaimed)

Pathak et al. (bioRxiv 2024.05.20.594987, "Pervasive ancestry bias in
variant effect predictors", Marsh lab, Edinburgh/NTU Singapore, revised
March 2025) rigorously demonstrated that clinical-trained variant effect
predictors (VEPs) — the tools clinicians actually use — show large,
ancestry-correlated distortions in predicted pathogenicity, driven by
uneven ClinVar representation across ancestry groups. Their discussion
explicitly recommends cross-checking clinical-trained VEP calls against
population-free VEPs when there's a discrepancy, but stops short of
building this into a tool. A second, independent paper (bioRxiv
2026.02.14.705914, Feb 2026) reaches a related but distinct conclusion:
ancestry-specific effects mostly disappear once you control for allele
frequency — meaning the *mechanism* of the bias, not just its existence,
is still actively being worked out. This is exactly the kind of live,
unsettled space where a new, well-scoped method has genuine room to
contribute.

**Nobody has published a tool that takes this recommendation and turns
it into a computable, per-variant, per-patient-ancestry score.** That's
the gap BAVCS fills.

## What's novel here (be precise about this in any write-up/application)

- Not novel: the *existence* of ancestry bias in VEPs (established,
  multiple papers).
- Not novel: the *idea* of cross-checking VEP classes (recommended by
  Pathak et al., not built).
- **Novel**: a formula that (a) operationalizes that cross-check
  automatically, (b) weights it by a newly-defined per-gene, per-ancestry
  representation-confidence metric R(g,a) — finer-grained than the
  whole-dataset ascertainment metric the source paper computed — and
  (c) outputs both a bias-adjusted score and an interpretable flag.

## Progress so far (real, not synthetic)

- **Core algorithm implemented and unit-tested**: `bavcs/score.py`
  (`compute_bavcs`, `representation_confidence`).
- **Real data pipeline confirmed working**: `bavcs/myvariant_client.py`
  queries the live myvariant.info API, which mirrors dbNSFP and returns
  real REVEL, CADD, MetaRNN, BayesDel, ESM-1b, AlphaMissense, PrimateAI,
  ClinPred, and VEST4 scores per variant — the same VEP panel used by the
  source paper.
- **Real, non-synthetic proof of concept**: `bavcs/data/braf_real_sample.json`
  holds 9 real BRAF missense variants with real dbNSFP scores, fetched
  live. `scripts/analyze_real_sample.py` computes real population-free
  vs. clinical-trained discordance on this data and produces
  `figures/braf_discordance.png`. Discordance ranged from 0.013 to 0.361
  across these 9 real variants — i.e. real VEP classes really do disagree
  by varying amounts on real variants, which is the load-bearing
  assumption BAVCS depends on.
- **What this does NOT yet show**: an actual ancestry bias effect. The
  representation-confidence weight R(g,a) is not yet computed on real
  data — that requires gnomAD ancestry-split allele counts and ClinVar
  classifications, which need a bulk-data environment with open internet
  access (this sandbox's network is limited to package registries, not
  biological databases). `bavcs/myvariant_client.py` and the data plan
  below are written to be run from your own machine or a cloud VM.

## Data plan (what's needed for full validation)

1. **VEP scores per variant.** `bavcs/myvariant_client.py` already does
   this for a gene panel via the myvariant.info API. For genome-wide
   coverage, a bulk dbNSFP download is more efficient than one API call
   per variant.
2. **Ancestry-labelled variant frequencies.** gnomAD v4.1 publishes
   per-ancestry allele counts, fully public. Gives the "observed variants
   per gene per ancestry" denominator for R(g,a).
3. **ClinVar classifications.** Public, fully downloadable (VCF or XML).
   Gives the "classified variants per gene per ancestry" numerator —
   approximate ancestry the same way the source paper did: matching a
   ClinVar entry's variant to its gnomAD ancestry-group allele frequency
   profile, since ClinVar itself doesn't label ancestry directly.

## Validation plan

1. **Retrospective concordance test** (next priority): the source
   paper released its full dataset at `osf.io/wz2sb`. Running BAVCS
   against their real data and checking whether it flags the same
   MAL/IMX cases they found most ancestry-distorted is the fastest path
   to a credible first result, since ground-truth labels already exist.
2. **MAVE cross-check.** For genes with multiplex assay of variant
   effect (MAVE) data (BRCA1, TP53, others in MaveDB), compare BAVCS's
   bias-adjusted score against experimentally measured functional
   effect — the closest thing to ground truth independent of any VEP's
   training data.
3. **Negative controls.** Confirm BAVCS does NOT flag variants in
   well-studied genes/ancestries as high-risk.


