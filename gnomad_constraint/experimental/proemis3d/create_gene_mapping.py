#!/usr/bin/env python3
"""
Create a transcript_id -> gene_id mapping TSV from a forward combined TSV.

Reads the forward TSV, collects unique transcript_ids, looks up each transcript's
parent gene ID via the Ensembl REST API, and writes transcript_id\tgene_id.

Usage:
  python -m gnomad_constraint.experimental.proemis3d.create_gene_mapping \\
    --tsv gnomad_constraint/experimental/proemis3d/forward_combined.tsv \\
    --output gnomad_constraint/experimental/proemis3d/gene_mapping.tsv

Requires: requests (pip install requests)
"""
import argparse
import sys
import time
from pathlib import Path

# Ensembl REST API: 15 requests per second recommended
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id/{id}?expand=1"
HEADERS = {"Content-Type": "application/json"}


def _unique_transcript_ids(tsv_path: str) -> list:
    """Read TSV and return unique transcript_id values (first column occurrence only if duplicate headers)."""
    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t", nrows=None, low_memory=False)
    first_tid_idx = next((i for i, c in enumerate(df.columns) if c == "transcript_id"), None)
    if first_tid_idx is None:
        raise ValueError("TSV has no 'transcript_id' column")
    series = df.iloc[:, first_tid_idx]
    series.name = "transcript_id"
    return series.dropna().unique().tolist()


def _fetch_gene_id(transcript_id: str, session) -> str | None:
    """Return parent gene ID for transcript from Ensembl REST API, or None on failure."""
    import requests

    url = ENSEMBL_LOOKUP.format(id=transcript_id)
    try:
        r = session.get(url, headers=HEADERS, timeout=30)
        r.raise_for_status()
        data = r.json()
        return data.get("Parent") or data.get("parent")
    except Exception as e:
        print(f"  Warning: {transcript_id}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(description="Create transcript_id -> gene_id mapping from forward TSV via Ensembl API.")
    parser.add_argument("--tsv", required=True, help="Path to forward_combined.tsv (or any TSV with transcript_id)")
    parser.add_argument("--output", required=True, help="Output TSV path (transcript_id, gene_id)")
    parser.add_argument("--delay", type=float, default=0.07, help="Seconds between API requests (default 0.07 ~ 14/sec)")
    args = parser.parse_args()

    try:
        import requests
    except ImportError:
        print("Install requests: pip install requests", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {args.tsv}...")
    transcript_ids = _unique_transcript_ids(args.tsv)
    print(f"Found {len(transcript_ids)} unique transcript(s).")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    with open(out_path, "w") as f:
        f.write("transcript_id\tgene_id\n")
        for i, tid in enumerate(sorted(transcript_ids)):
            gene_id = _fetch_gene_id(tid, session)
            if gene_id:
                f.write(f"{tid}\t{gene_id}\n")
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(transcript_ids)}")
            time.sleep(args.delay)

    print(f"Wrote {out_path}")
    with open(out_path) as f:
        next(f)  # header
        n_genes = len(set(line.split("\t")[1].strip() for line in f if "\t" in line))
    print(f"Unique genes: {n_genes}")


if __name__ == "__main__":
    main()
