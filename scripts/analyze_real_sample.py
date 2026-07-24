"""
Generates figures/braf_discordance.png from the REAL BRAF dbNSFP sample
(bavcs/data/braf_real_sample.json). No synthetic data used in this script.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bavcs import load_real_sample, compute_real_discordance


def main():
    raw = load_real_sample()
    results = compute_real_discordance(raw)
    results.sort(key=lambda r: r.discordance)

    labels = [f"{r.hgvsp}" for r in results]
    s_pf = [r.s_pf for r in results]
    s_ct = [r.s_ct for r in results]
    discordance = [r.discordance for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    x = range(len(results))
    ax1.plot(x, s_pf, "o-", label="Population-free mean (ESM-1b)", color="#2563eb")
    ax1.plot(x, s_ct, "o-", label="Clinical-trained mean (REVEL, CADD, MetaRNN,\nBayesDel, ClinPred, VEST4 where available)", color="#dc2626")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax1.set_ylabel("Normalized score (0-1, higher = more damaging)")
    ax1.set_title("BRAF missense variants: real dbNSFP scores\nby VEP class (myvariant.info, fetched 2026-07-18)")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.set_ylim(-0.05, 1.05)

    bars = ax2.bar(x, discordance, color="#7c3aed")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax2.set_ylabel("Discordance |S_ct - S_pf|")
    ax2.set_title("Population-free vs. clinical-trained\ndiscordance per real variant")
    ax2.axhline(0.15, color="black", linestyle="--", linewidth=1,
                label="BAVCS discordance threshold (0.15)")
    ax2.legend(fontsize=8)

    fig.suptitle(
        "BAVCS pipeline validated on real data (single gene, no ancestry weighting yet)",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    fig.savefig("figures/braf_discordance.png", dpi=180, bbox_inches="tight")
    print("Saved figures/braf_discordance.png")

    print("\nReal discordance results:")
    for r in results:
        print(f"  {r.hgvsp:18s} S_pf={r.s_pf:.3f}  S_ct={r.s_ct:.3f}  discordance={r.discordance:.3f}")


if __name__ == "__main__":
    main()
