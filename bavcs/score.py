"""
BAVCS: Bias-Aware Variant Concordance Score
=============================================

A per-variant, per-ancestry method for flagging when a clinical-trained
variant effect predictor's (VEP) pathogenicity call may be inflated or
deflated by ancestry representation bias in its training data.

Background
----------
Pathak, Bora, Badonyi, Livesey, Ngeow & Marsh (bioRxiv 2024.05.20.594987,
"Pervasive ancestry bias in variant effect predictors") showed that
clinical-trained VEPs (e.g. BayesDel, MetaRNN, ClinPred) show large,
ancestry-correlated swings in predicted pathogenicity burden, driven by
how thoroughly each ancestry group's variants have been clinically
classified (ClinVar ascertainment bias). Population-free VEPs (e.g.
CPT-1, ESM-1b, GEMME) are largely immune to this because they never see
clinical labels or population allele frequencies.

Their recommendation: when a clinical-trained VEP's call disagrees with
a population-free VEP's call, consider whether the discrepancy is
ancestry-driven. This module operationalizes that recommendation as a
computable per-variant score, using per-gene ancestry representation
confidence (rather than the whole-dataset-level bias they measured) as
the weighting factor.

This is a prototype / reference implementation. It demonstrates the
algorithm on synthetic data — see demo.py. To run on real variants you
need:
  1. Per-VEP scores for the variant (e.g. from dbNSFP, or querying each
     VEP's own database/API).
  2. Gene x ancestry representation data (e.g. computed from gnomAD +
     ClinVar, as in the source paper's Methods).
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import statistics


# ---------------------------------------------------------------------------
# VEP classification (from Livesey & Marsh's scheme, used by the source paper)
# ---------------------------------------------------------------------------

POPULATION_FREE = "population_free"
POPULATION_TUNED = "population_tuned"
CLINICAL_TRAINED = "clinical_trained"

VEP_CLASSES: Dict[str, str] = {
    # population-free: evolutionary / protein-language-model based
    "CPT-1": POPULATION_FREE,
    "ESM-1b": POPULATION_FREE,
    "GEMME": POPULATION_FREE,
    # population-tuned: calibrated with population variants, not clinical labels
    "AlphaMissense": POPULATION_TUNED,
    "popEVE": POPULATION_TUNED,
    "PrimateAI": POPULATION_TUNED,
    # clinical-trained: trained directly on ClinVar-style pathogenic/benign labels
    "REVEL": CLINICAL_TRAINED,
    "BayesDel": CLINICAL_TRAINED,
    "MetaRNN": CLINICAL_TRAINED,
    "ClinPred": CLINICAL_TRAINED,
    "VEST4": CLINICAL_TRAINED,
    "CADD": CLINICAL_TRAINED,
}


@dataclass
class VariantScores:
    """Raw rank-normalised (0-1) scores for one variant, keyed by VEP name."""
    variant_id: str
    gene: str
    scores: Dict[str, float]  # VEP name -> rank-normalised score in [0, 1]

    def class_mean(self, vep_class: str) -> Optional[float]:
        vals = [
            s for name, s in self.scores.items()
            if VEP_CLASSES.get(name) == vep_class
        ]
        return statistics.mean(vals) if vals else None


@dataclass
class BAVCSResult:
    variant_id: str
    gene: str
    ancestry: str
    s_pf: Optional[float]          # population-free class mean
    s_ct: Optional[float]          # clinical-trained class mean
    discordance: Optional[float]   # |S_ct - S_pf|
    representation_confidence: float  # R(g, a) in [0, 1]
    bavcs: Optional[float]         # final bias-adjusted score
    flag: str                      # human-readable risk flag
    rationale: str


# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def representation_confidence(
    gene: str,
    ancestry: str,
    classified_variants_in_gene_ancestry: int,
    observed_variants_in_gene_ancestry: int,
) -> float:
    """
    R(g, a): fraction of a gene's missense variants observed in a given
    ancestry group that have ANY clinical classification in ClinVar.

    This is the per-gene, per-ancestry analogue of the whole-dataset
    "% of variants with ClinVar labels per ancestry group" metric used
    in Figure 1B of the source paper (Pathak et al. 2024/2025) — but
    computed narrowly enough to weight a single variant's confidence,
    not just describe the dataset in aggregate.

    In production this needs:
      - observed_variants_in_gene_ancestry: count of missense variants
        in `gene` seen in ancestry-labelled cohorts (e.g. gnomAD, split
        by the same ancestry groupings as the source paper: AFR, AMR,
        ASJ, EAS, FIN, MID, NFE, SAS, plus MCPS/SG10K groups for
        under-served populations).
      - classified_variants_in_gene_ancestry: subset of those with a
        ClinVar (likely) pathogenic/benign classification.
    """
    if observed_variants_in_gene_ancestry <= 0:
        return 0.0
    return min(
        1.0,
        classified_variants_in_gene_ancestry / observed_variants_in_gene_ancestry,
    )


def compute_bavcs(
    variant: VariantScores,
    ancestry: str,
    representation_confidence_value: float,
    discordance_threshold: float = 0.15,
    low_confidence_threshold: float = 0.10,
) -> BAVCSResult:
    """
    Compute the Bias-Aware Variant Concordance Score for one variant,
    for a patient of a given ancestry.
    """
    s_pf = variant.class_mean(POPULATION_FREE)
    s_ct = variant.class_mean(CLINICAL_TRAINED)

    if s_pf is None or s_ct is None:
        return BAVCSResult(
            variant_id=variant.variant_id,
            gene=variant.gene,
            ancestry=ancestry,
            s_pf=s_pf,
            s_ct=s_ct,
            discordance=None,
            representation_confidence=representation_confidence_value,
            bavcs=None,
            flag="INSUFFICIENT_DATA",
            rationale=(
                "Need at least one population-free and one clinical-trained "
                "VEP score to compute BAVCS."
            ),
        )

    discordance = abs(s_ct - s_pf)
    r = representation_confidence_value

    # Confidence-weighted blend: the clinical-trained score's contribution
    # is discounted in direct proportion to how poorly this gene has been
    # clinically characterised in this patient's ancestry group.
    bavcs = (s_pf + r * s_ct) / (1 + r)

    high_bias_risk = discordance >= discordance_threshold and r <= low_confidence_threshold

    if high_bias_risk:
        flag = "HIGH_BIAS_RISK"
        rationale = (
            f"Clinical-trained VEPs diverge from population-free VEPs by "
            f"{discordance:.2f} (rank scale), and gene {variant.gene} has low "
            f"clinical curation depth in ancestry '{ancestry}' "
            f"(R={r:.2f}). The clinical-trained call may be inflated or "
            f"deflated by ancestry ascertainment bias rather than genuine "
            f"pathogenicity. Recommend weighting the population-free score "
            f"(S_pf={s_pf:.2f}) more heavily than the clinical-trained "
            f"consensus (S_ct={s_ct:.2f})."
        )
    elif discordance >= discordance_threshold:
        flag = "DISCORDANT_BUT_WELL_REPRESENTED"
        rationale = (
            f"Predictor classes disagree (discordance={discordance:.2f}) but "
            f"gene {variant.gene} is reasonably well-characterised in "
            f"ancestry '{ancestry}' (R={r:.2f}), so the discordance is less "
            f"likely to be an ancestry-bias artifact. Still worth manual "
            f"review, but not flagged as a bias-specific concern."
        )
    else:
        flag = "LOW_BIAS_RISK"
        rationale = (
            f"Population-free and clinical-trained VEPs agree closely "
            f"(discordance={discordance:.2f}). No strong signal of "
            f"ancestry-driven distortion for this variant."
        )

    return BAVCSResult(
        variant_id=variant.variant_id,
        gene=variant.gene,
        ancestry=ancestry,
        s_pf=s_pf,
        s_ct=s_ct,
        discordance=discordance,
        representation_confidence=r,
        bavcs=bavcs,
        flag=flag,
        rationale=rationale,
    )
