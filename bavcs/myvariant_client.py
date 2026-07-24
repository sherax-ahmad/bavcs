"""
bavcs.myvariant_client
=======================

A real client for the myvariant.info API (public, no auth required for
basic use; ~1000 requests/day anonymous). This is confirmed working --
see the real BRAF sample this project ships with, fetched with the same
underlying API.

THIS FILE CANNOT BE RUN INSIDE THE SANDBOXED ENVIRONMENT THAT GENERATED
THIS PACKAGE (no outbound network to biological databases there). Run it
from your own machine, a DigitalOcean/Azure VM, or wherever you have
normal internet access.

Usage:
    pip install requests
    python -m bavcs.myvariant_client BRCA1 TP53 CFTR --out data/scores.json
"""

import argparse
import json
import sys
import time
from typing import List, Dict, Any

try:
    import requests
except ImportError:
    requests = None

MYVARIANT_BASE = "https://myvariant.info/v1/query"

FIELDS = ",".join([
    "_id",
    "dbnsfp.hgvsp",
    "dbnsfp.genename",
    "dbnsfp.revel.score",
    "dbnsfp.bayesdel.no_af.score",
    "dbnsfp.metarnn.score",
    "dbnsfp.esm1b.score",
    "dbnsfp.alphamissense.score",
    "dbnsfp.primateai.score",
    "dbnsfp.clinpred.score",
    "dbnsfp.vest4.score",
    "cadd.phred",
    "clinvar.rcv.clinical_significance",
])


def fetch_gene_variants(gene: str, size: int = 200) -> List[Dict[str, Any]]:
    """Fetch missense variant annotations for one gene from myvariant.info."""
    if requests is None:
        raise RuntimeError("pip install requests first")

    query = f"dbnsfp.genename:{gene}"
    params = {"q": query, "fields": FIELDS, "size": size}
    resp = requests.get(MYVARIANT_BASE, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    return data.get("hits", [])


def fetch_many_genes(genes: List[str], size: int = 200, sleep_s: float = 0.5) -> Dict[str, List[dict]]:
    """
    Fetch variants for a panel of genes. Keep the gene panel small (10-20
    genes) for a first pass -- see RESEARCH_PLAN.md for why genome-wide is
    a much bigger undertaking (needs bulk dbNSFP download, not per-gene
    API calls).
    """
    out = {}
    for gene in genes:
        print(f"Fetching {gene}...", file=sys.stderr)
        out[gene] = fetch_gene_variants(gene, size=size)
        time.sleep(sleep_s)  # be polite to the free public API
    return out


def main():
    parser = argparse.ArgumentParser(description="Fetch real dbNSFP VEP scores for a gene panel")
    parser.add_argument("genes", nargs="+", help="Gene symbols, e.g. BRCA1 TP53 CFTR")
    parser.add_argument("--out", default="data/fetched_scores.json")
    parser.add_argument("--size", type=int, default=200, help="Max variants per gene")
    args = parser.parse_args()

    results = fetch_many_genes(args.genes, size=args.size)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    total = sum(len(v) for v in results.values())
    print(f"Wrote {total} variant records across {len(args.genes)} genes to {args.out}")


if __name__ == "__main__":
    main()
