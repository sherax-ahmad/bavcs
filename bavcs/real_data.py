"""
bavcs.real_data
================

Utilities for loading and normalizing REAL dbNSFP-derived VEP scores
(as opposed to the synthetic examples used only for algorithm sanity-checks).

Score direction note
---------------------
Not all VEPs point the same way. Most (REVEL, CADD, MetaRNN, ClinPred,
AlphaMissense, PrimateAI, VEST4, BayesDel) report higher = more damaging.
ESM-1b is a pseudo-log-likelihood: LOWER (more negative) = more damaging.
This module flips ESM-1b's sign before normalization so that, after
processing, higher always means "predicted more damaging" across every
VEP -- required before scores from different VEPs can be combined or
compared.
"""

import json
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from .score import VariantScores, VEP_CLASSES, POPULATION_FREE, CLINICAL_TRAINED

# VEPs whose raw score direction must be flipped so higher = more damaging
INVERTED_DIRECTION_VEPS = {"ESM-1b"}

DATA_DIR = Path(__file__).parent / "data"


def load_real_sample(path: Path = None) -> List[dict]:
    """Load the real, dbNSFP-sourced BRAF variant sample (see data/README)."""
    path = path or (DATA_DIR / "braf_real_sample.json")
    with open(path) as f:
        payload = json.load(f)
    return payload["variants"]


def normalize_real_sample(raw_variants: List[dict]) -> List[VariantScores]:
    """
    Min-max normalize each VEP's raw scores to [0, 1] ACROSS THIS SAMPLE ONLY.

    Important limitation: proper rank-normalization (as used in the source
    paper and as BAVCS is designed for) requires each VEP's score
    distribution across the WHOLE genome (or at least a large variant set),
    not just 9 variants in one gene. Min-max normalization within this small
    real sample is a reasonable stand-in for validating the pipeline
    end-to-end on real data, but the resulting BAVCS numbers here should be
    read as a mechanism check, not a calibrated bias estimate. Real
    calibration needs the dbNSFP-genome-wide rankscore columns (dbNSFP
    ships these directly, e.g. `revel_rankscore`), which requires a bulk
    dbNSFP download -- see scripts/fetch_real_data.py.
    """
    # collect raw values per VEP across the sample
    by_vep: Dict[str, List[float]] = {}
    for v in raw_variants:
        for vep, val in v["scores"].items():
            by_vep.setdefault(vep, []).append(val)

    # compute min/max per VEP, applying direction flip first
    ranges = {}
    for vep, vals in by_vep.items():
        adj = [-x for x in vals] if vep in INVERTED_DIRECTION_VEPS else vals
        ranges[vep] = (min(adj), max(adj))

    normalized = []
    for v in raw_variants:
        norm_scores = {}
        for vep, val in v["scores"].items():
            adj = -val if vep in INVERTED_DIRECTION_VEPS else val
            lo, hi = ranges[vep]
            norm_scores[vep] = (adj - lo) / (hi - lo) if hi > lo else 0.5
        normalized.append(
            VariantScores(
                variant_id=v["variant_id"],
                gene=v["gene"],
                scores=norm_scores,
            )
        )
    return normalized


@dataclass
class DiscordanceResult:
    variant_id: str
    hgvsp: str
    s_pf: float
    s_ct: float
    discordance: float


def compute_real_discordance(raw_variants: List[dict]) -> List[DiscordanceResult]:
    """
    Real, non-synthetic discordance computation: population-free (ESM-1b
    only, in this sample) vs. clinical-trained (REVEL, CADD, MetaRNN,
    BayesDel, ClinPred, VEST4 where available) mean normalized score, per
    real BRAF variant.

    This does NOT yet include the ancestry representation-confidence
    weighting R(g,a) -- that requires gnomAD ancestry-split allele counts
    and ClinVar classifications, which this sandbox cannot bulk-fetch (see
    RESEARCH_PLAN.md). What's below is the real, verifiable half of BAVCS:
    the discordance signal itself, computed on real predictor output.
    """
    normalized = normalize_real_sample(raw_variants)
    hgvsp_by_id = {v["variant_id"]: v["hgvsp"] for v in raw_variants}

    results = []
    for v in normalized:
        s_pf = v.class_mean(POPULATION_FREE)
        s_ct = v.class_mean(CLINICAL_TRAINED)
        if s_pf is None or s_ct is None:
            continue
        results.append(
            DiscordanceResult(
                variant_id=v.variant_id,
                hgvsp=hgvsp_by_id.get(v.variant_id, ""),
                s_pf=s_pf,
                s_ct=s_ct,
                discordance=abs(s_ct - s_pf),
            )
        )
    return results
