#!/usr/bin/env python3
import os, re, csv, sys
from pathlib import Path

# --- CONFIG: update these two paths for your run ---
DATA_DIR   = Path("/mount/britton/Jose/Jobs/Jason/2025-07-29-novel-isolate-closed-genome-nextflow-run/data")
OUTPUT_CSV = Path("/mount/britton/Jose/Jobs/Jason/2025-07-29-novel-isolate-closed-genome-nextflow-run/workflow/nextflow_readfile.csv")

ILLUMINA_R1_RE = re.compile(r"^(?P<sample>.+?)(?:_S\d+)?_R1(?:_\d+)?\.fastq\.gz$")
ILLUMINA_R2_RE = re.compile(r"^(?P<sample>.+?)(?:_S\d+)?_R2(?:_\d+)?\.fastq\.gz$")
# Heuristic for ONT: allow names like Pseudo_nanopore.fastq.gz or Pseudo_ONT.fastq.gz
ONT_RE         = re.compile(r"^(?P<sample>.+?)(?:_(nanopore|ont|ONT).*)?\.fastq\.gz$")

def abspath(p: Path) -> str:
    return str(p.resolve())

def main():
    if not DATA_DIR.is_dir():
        print(f"ERROR: Data directory not found: {DATA_DIR}", file=sys.stderr)
        sys.exit(1)

    # Scan only top-level; adjust to rglob('*.fastq.gz') if you want recursion
    files = sorted([p for p in DATA_DIR.glob("*.fastq.gz")])

    samples = {}  # sample -> {'r1': Path|None, 'r2': Path|None, 'ont': Path|None}

    # First pass: capture Illumina pairs
    for p in files:
        name = p.name
        m1 = ILLUMINA_R1_RE.match(name)
        m2 = ILLUMINA_R2_RE.match(name)
        if m1:
            sid = m1.group("sample")
            samples.setdefault(sid, {"r1": None, "r2": None, "ont": None})
            # prefer the largest file if duplicates exist (multi-lane)
            if samples[sid]["r1"] is None or p.stat().st_size > samples[sid]["r1"].stat().st_size:
                samples[sid]["r1"] = p
            continue
        if m2:
            sid = m2.group("sample")
            samples.setdefault(sid, {"r1": None, "r2": None, "ont": None})
            if samples[sid]["r2"] is None or p.stat().st_size > samples[sid]["r2"].stat().st_size:
                samples[sid]["r2"] = p
            continue

    # Second pass: assign ONT reads by sample (match prefix before _nanopore/_ont or full stem)
    # Strategy: For each ONT file, infer a candidate sample and attach if that sample exists.
    for p in files:
        name = p.name
        if ILLUMINA_R1_RE.match(name) or ILLUMINA_R2_RE.match(name):
            continue
        if not name.endswith(".fastq.gz"):
            continue
        m = ONT_RE.match(name)
        if not m:
            continue
        candidate = m.group("sample")
        # If exact sample exists, use it
        if candidate in samples:
            if samples[candidate]["ont"] is None or p.stat().st_size > samples[candidate]["ont"].stat().st_size:
                samples[candidate]["ont"] = p
            continue
        # Fallback: if only one sample exists and this ONT file is present, assign it
        if len(samples) == 1:
            only_sid = next(iter(samples))
            if samples[only_sid]["ont"] is None or p.stat().st_size > samples[only_sid]["ont"].stat().st_size:
                samples[only_sid]["ont"] = p

    # Write CSV in the schema Nextflow expects: sample,r1,r2,ont
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(OUTPUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "r1", "r2", "ont"])
        for sid, paths in sorted(samples.items()):
            r1, r2, ont = paths["r1"], paths["r2"], paths["ont"]
            if not (r1 and r2 and ont):
                # Print a clear message but keep going for other samples
                missing = ",".join(k for k, v in paths.items() if v is None)
                print(f"SKIP {sid}: missing {missing}", file=sys.stderr)
                continue
            w.writerow([sid, abspath(r1), abspath(r2), abspath(ont)])
            written += 1

    print(f"Wrote {written} row(s) to {OUTPUT_CSV}")
    if written == 0:
        print("No complete {R1,R2,ONT} sets found — check filenames.", file=sys.stderr)

if __name__ == "__main__":
    main()

