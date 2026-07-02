#!/usr/bin/env python3
import os, re, csv, sys
from pathlib import Path

# --- CONFIG: update these two paths for run ------ 
DATA_DIR   = Path("/mount/britton/Jose/Jobs/Jason/2025-07-29-novel-isolate-closed-genome-nextflow-run/data")
OUTPUT_CSV = Path("/mount/britton/Jose/Jobs/Jason/2025-12-15-Phylogram-Tree-Genus-Nextflow-Run/workflow/nextflow_readfile.csv")

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

    samples = {}  # sample -> {'ont': Path|None}

    # Strategy: For each ONT file, infer a candidate sample and attach if that sample exists.
    for p in files: 
        name = p.name 
        m = ONT_RE.match(name)
        if not m:
            continue 

        sid = m.group("sample")

        # Initialize the sample if we haven't seen it yet 
        if sid not in samples:
            samples[sid] = {"ont": None}

        # Keep largest file if duplicate exists 
        if samples[sid]["ont"] is None or p.stat().st_size > samples[sid]["ont"].stat().st_size:
            samples[sid]["ont"] = p 

    # Write CSV in the schema Nextflow expects: sample,r1,r2,ont
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(OUTPUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "ont"])
        for sid, paths in sorted(samples.items()):
            ont = paths["ont"]
            if not ont:
                # Print a clear message but keep going for other samples
                missing = ",".join(k for k, v in paths.items() if v is None)
                print(f"SKIP {sid}: missing {missing}", file=sys.stderr)
                continue
            w.writerow([sid, abspath(ont)])
            written += 1

    print(f"Wrote {written} row(s) to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()

