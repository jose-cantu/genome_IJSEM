#!/usr/bin/env python
# Emits: "<genus> <family>" from the GTDB-Tk TSV 'classification' column (first data row)

import sys, csv

def extract_tokens(classification):
    classification = (classification or "").replace("\r", "")
    if not classification:
        return "", ""
    gen = fam = ""
    for tok in (t.strip() for t in classification.split(";")):
        if not gen and tok.startswith("g__"):
            gen = tok[3:]
        if not fam and tok.startswith("f__"):
            fam = tok[3:]
        if gen and fam:
            break
    return gen, fam

def main():
    if len(sys.argv) < 2:
        sys.stdout.write("\n")
        return 0

    fn = sys.argv[1]
    # Py2: binary file for csv; Py3: text is fine (no newline kwarg for Py2)
    is_py2 = (sys.version_info[0] < 3)
    fh = open(fn, "rb" if is_py2 else "r")
    try:
        r = csv.DictReader(fh, delimiter="\t")
        gen = fam = ""
        for row in r:
            # header-robust: find the key equal to 'classification' ignoring case/whitespace
            key = None
            for k in row.keys():
                if k and k.strip().lower() == "classification":
                    key = k
                    break
            c = row.get(key, "") if key else ""
            gen, fam = extract_tokens(c)
            break  # first data row suffices
        sys.stdout.write("%s %s\n" % (gen, fam))   # Py2/3-safe
    finally:
        fh.close()

if __name__ == "__main__":
    main()

