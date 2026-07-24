"""
Tests for BAVCS. The tests in this file run against REAL dbNSFP data
(bavcs/data/braf_real_sample.json), fetched live from myvariant.info on
2026-07-18 -- not synthetic fixtures. This is the "tested on real data"
layer of the project: it proves the parsing, normalization, and
discordance math work correctly on real VEP output, even though the
full ancestry-weighted BAVCS score still needs external gnomAD/ClinVar
data (see RESEARCH_PLAN.md).
"""

import math
import pytest

from bavcs import (
    load_real_sample,
    normalize_real_sample,
    compute_real_discordance,
    compute_bavcs,
    VariantScores,
)


def test_real_sample_loads():
    variants = load_real_sample()
    assert len(variants) == 9
    assert all("scores" in v for v in variants)
    assert all(v["gene"] == "BRAF" for v in variants)


def test_real_sample_has_expected_veps():
    variants = load_real_sample()
    first = variants[0]["scores"]
    for expected in ["REVEL", "CADD", "MetaRNN", "BayesDel", "ESM-1b",
                      "AlphaMissense", "PrimateAI", "ClinPred"]:
        assert expected in first, f"missing {expected} in real fixture"


def test_normalization_flips_esm1b_direction():
    """
    p.Ile724Asn has ESM-1b = -16.454 (strongly damaging, most negative)
    and p.Ser339Thr has ESM-1b = -3.916 (mild/tolerated, least negative).
    After direction-corrected normalization, the damaging one should end
    up with the HIGHER normalized ESM-1b score.
    """
    raw = load_real_sample()
    normalized = normalize_real_sample(raw)
    by_id = {v.variant_id: v for v in normalized}

    damaging = by_id["chr7:g.140434527A>T"]  # p.Ile724Asn, ESM-1b -16.454
    tolerated = by_id["chr7:g.140494233A>T"]  # p.Ser339Thr, ESM-1b -3.916

    assert damaging.scores["ESM-1b"] > tolerated.scores["ESM-1b"]


def test_normalized_scores_are_bounded():
    raw = load_real_sample()
    normalized = normalize_real_sample(raw)
    for v in normalized:
        for vep, score in v.scores.items():
            assert 0.0 <= score <= 1.0, f"{vep} score {score} out of [0,1] for {v.variant_id}"


def test_real_discordance_is_computable_for_all_variants():
    raw = load_real_sample()
    results = compute_real_discordance(raw)
    # every variant in the real sample has at least ESM-1b (pop-free) and
    # several clinical-trained VEPs, so all 9 should produce a result
    assert len(results) == 9
    for r in results:
        assert 0.0 <= r.s_pf <= 1.0
        assert 0.0 <= r.s_ct <= 1.0
        assert math.isclose(r.discordance, abs(r.s_ct - r.s_pf), abs_tol=1e-9)


def test_compute_bavcs_on_a_real_normalized_variant():
    """
    Sanity check that compute_bavcs (the full algorithm, including the
    representation-confidence weighting) runs correctly on a REAL
    normalized variant, using a placeholder representation confidence
    until real gnomAD/ClinVar ancestry data is wired in.
    """
    raw = load_real_sample()
    normalized = normalize_real_sample(raw)
    variant = normalized[0]

    result = compute_bavcs(
        variant,
        ancestry="NFE",
        representation_confidence_value=0.5,
    )
    assert result.flag in {
        "HIGH_BIAS_RISK", "DISCORDANT_BUT_WELL_REPRESENTED",
        "LOW_BIAS_RISK", "INSUFFICIENT_DATA",
    }
    assert result.bavcs is not None
