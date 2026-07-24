"""
BAVCS demo — SYNTHETIC DATA ONLY. For a REAL-DATA run, see
scripts/analyze_real_sample.py instead, which uses actual fetched
dbNSFP scores for 9 real BRAF variants (bavcs/data/braf_real_sample.json).
This synthetic demo exists only to build intuition for how the flags
behave across ancestries; it is not evidence of anything.

This script demonstrates the algorithm's behavior using hand-constructed
example variants whose score patterns are modeled on the qualitative
trends reported in Pathak et al. (bioRxiv 2024.05.20.594987): clinical-
trained VEPs (BayesDel, MetaRNN, ClinPred) showing inflated pathogenicity
predictions for underrepresented ancestries (the paper specifically
flags Malay and Indigenous Mexican groups), while population-free VEPs
(CPT-1, ESM-1b, GEMME) stay stable across ancestries for the same
underlying variant.

This is NOT real variant data. To run on real variants, replace the
VariantScores objects below with scores pulled from dbNSFP (or each
VEP's own precomputed score files) and replace the representation
confidence numbers with real gene x ancestry ClinVar/gnomAD counts.
"""

from bavcs import VariantScores, compute_bavcs


def print_result(r):
    print(f"\nVariant: {r.variant_id}  (gene {r.gene}, ancestry {r.ancestry})")
    print(f"  Population-free mean score : {r.s_pf:.3f}" if r.s_pf is not None else "  Population-free mean score : n/a")
    print(f"  Clinical-trained mean score: {r.s_ct:.3f}" if r.s_ct is not None else "  Clinical-trained mean score: n/a")
    print(f"  Discordance                : {r.discordance:.3f}" if r.discordance is not None else "  Discordance: n/a")
    print(f"  Representation confidence  : {r.representation_confidence:.3f}")
    print(f"  BAVCS (bias-adjusted score): {r.bavcs:.3f}" if r.bavcs is not None else "  BAVCS: n/a")
    print(f"  FLAG: {r.flag}")
    print(f"  Rationale: {r.rationale}")


if __name__ == "__main__":
    # --- Case 1: a variant in a well-studied gene (e.g. BRCA1-like), seen
    # in a European (NFE) patient. Clinical-trained VEPs and population-free
    # VEPs roughly agree; the gene is well characterised for this ancestry.
    v1 = VariantScores(
        variant_id="chr17:g.43094692G>A (synthetic)",
        gene="BRCA1_like",
        scores={
            "CPT-1": 0.78, "ESM-1b": 0.75, "GEMME": 0.80,
            "AlphaMissense": 0.81, "popEVE": 0.77,
            "BayesDel": 0.83, "MetaRNN": 0.85, "ClinPred": 0.82,
        },
    )
    r1 = compute_bavcs(
        v1, ancestry="NFE",
        representation_confidence_value=0.85,  # well-characterised for NFE
    )
    print_result(r1)

    # --- Case 2: the SAME gene, SAME underlying variant class of severity,
    # but the patient is of Malay ancestry, and the gene is far less
    # clinically characterised for this population (low R). Clinical-
    # trained VEPs are pulled toward "damaging" (as the source paper found
    # for MAL specifically), while population-free VEPs stay put.
    v2 = VariantScores(
        variant_id="chr17:g.43094692G>A (synthetic, MAL cohort)",
        gene="BRCA1_like",
        scores={
            "CPT-1": 0.55, "ESM-1b": 0.52, "GEMME": 0.58,
            "AlphaMissense": 0.60, "popEVE": 0.54,
            "BayesDel": 0.88, "MetaRNN": 0.91, "ClinPred": 0.86,
        },
    )
    r2 = compute_bavcs(
        v2, ancestry="MAL",
        representation_confidence_value=0.08,  # poorly characterised for MAL
    )
    print_result(r2)

    # --- Case 3: discordant scores, but in a gene that happens to be
    # reasonably well studied even for an underrepresented ancestry
    # (e.g. a pharmacogenomic gene with dedicated diversity studies).
    v3 = VariantScores(
        variant_id="chr10:g.some_locus (synthetic, AFR cohort)",
        gene="CYP_like",
        scores={
            "CPT-1": 0.40, "ESM-1b": 0.44, "GEMME": 0.38,
            "AlphaMissense": 0.50, "popEVE": 0.46,
            "BayesDel": 0.66, "MetaRNN": 0.70, "ClinPred": 0.64,
        },
    )
    r3 = compute_bavcs(
        v3, ancestry="AFR",
        representation_confidence_value=0.55,  # moderately well studied
    )
    print_result(r3)

    print("\n" + "=" * 70)
    print("Case 1 (NFE, well-represented gene)   -> expect LOW_BIAS_RISK")
    print("Case 2 (MAL, poorly-represented gene) -> expect HIGH_BIAS_RISK")
    print("Case 3 (AFR, moderately represented)  -> expect DISCORDANT_BUT_WELL_REPRESENTED")
