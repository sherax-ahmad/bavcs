# BAVCS — Bias-Aware Variant Concordance Score

A per-variant, per-ancestry method for flagging when a clinical-trained
variant effect predictor's (VEP) pathogenicity call may be distorted by
ancestry representation bias in its training data, rather than reflecting
genuine biology.

## Why this exists

Pathak et al. (bioRxiv 2024.05.20.594987, Marsh lab, revised March 2025)
showed that clinical-trained VEPs (REVEL, BayesDel, MetaRNN, ClinPred...)
have ancestry-correlated pathogenicity distortions, driven by uneven
ClinVar representation, while population-free VEPs (ESM-1b, CPT-1, GEMME)
are largely immune. Their paper recommends cross-checking the two VEP
classes when they disagree — but stops at the recommendation. **BAVCS
turns that recommendation into a computable per-variant score.**

See [`RESEARCH_PLAN.md`](RESEARCH_PLAN.md) for the full research plan,
what's novel here, and what real-data validation still needs to happen.

## What's real vs. what's a placeholder right now

| Component | Status |
|---|---|
| Core BAVCS algorithm (`bavcs/score.py`) | Implemented, unit-tested |
| Real dbNSFP data pipeline (`bavcs/myvariant_client.py`) | Implemented, confirmed working against the live myvariant.info API |
| Real BRAF sample + discordance analysis (`bavcs/data/`, `scripts/analyze_real_sample.py`) | **Real data**, 9 real BRAF missense variants, real REVEL/CADD/MetaRNN/BayesDel/ESM-1b/AlphaMissense/PrimateAI/ClinPred/VEST4 scores |
| Ancestry representation confidence `R(g,a)` | **Not yet computed on real data** — needs gnomAD ancestry-split allele counts + ClinVar classifications, which requires an environment with open internet access (see below) |
| Validation against the source paper's own dataset | Not started — next priority, see `RESEARCH_PLAN.md` |

## Quick start

```bash
pip install -e ".[dev]"
pytest tests/ -v                        # 6 tests, run against real fetched data
python scripts/analyze_real_sample.py   # regenerates figures/braf_discordance.png
```

## Extending to real, ancestry-weighted BAVCS

This was built inside a sandboxed environment with no outbound network
access to biological databases (only package registries) — so gene-panel
data collection had to stop at a small, real, single-gene proof of
concept. To go further:

```bash
pip install requests
python -m bavcs.myvariant_client BRCA1 TP53 CFTR MLH1 --out data/gene_panel.json
```

Run that from your own machine or a cloud VM (your GitHub Student Pack
DigitalOcean/Azure credits work well here). Then pull gnomAD ancestry-split
allele counts and ClinVar classifications for the same genes to compute
`R(g,a)` for real — see `RESEARCH_PLAN.md` for the exact data sources and
validation plan (including using the source paper's own published dataset
at `osf.io/wz2sb` as a ground-truth check).

## Repository layout

```
bavcs/
  score.py            # core BAVCS algorithm
  real_data.py         # real-data loading, normalization, discordance
  myvariant_client.py  # real API client (run outside this sandbox)
  data/
    braf_real_sample.json  # real dbNSFP data, 9 BRAF variants
tests/
  test_score.py        # unit tests against the real fixture
scripts/
  analyze_real_sample.py  # generates figures/braf_discordance.png
demo.py                 # synthetic walkthrough of the algorithm (for intuition only)
RESEARCH_PLAN.md        # full research plan, novelty case, validation plan
```
