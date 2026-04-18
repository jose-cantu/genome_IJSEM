#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import json
import os
import re
import shutil
import sys

def _read_jsonl(path):
    """Yield JSON objects from a .jsonl file (skip blank/bad lines)."""
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except ValueError:
                continue

def _coerce_assembly_record(rec):
    """
    Normalize record shape across Datasets variants:
      - assembly_data_report.json(l): top-level keys 'accession', 'organism'
      - summary outputs: may be nested under 'assembly'
    """
    if isinstance(rec, dict) and 'assembly' in rec and isinstance(rec['assembly'], dict):
        rec = rec['assembly']
    return rec if isinstance(rec, dict) else {}

def _get_accession(rec):
    # Try common keys
    return (rec.get('accession') or
            rec.get('currentAccession') or
            rec.get('assembly_accession') or
            rec.get('pairedAccession'))

def _get_organism(rec):
    org_obj = rec.get('organism') or {}
    name = org_obj.get('organismName') or org_obj.get('organism_name')
    if not name:
        # summary-style fallback
        org2 = rec.get('org') or {}
        name = org2.get('sci_name') or org2.get('title')
    return name or ''

def load_metadata(meta_path):
    """Return {accession: organism_name} mapping from .jsonl or .json."""
    mapping = {}
    base, ext = os.path.splitext(meta_path)
    ext = ext.lower()

    if ext == '.jsonl':
        itr = _read_jsonl(meta_path)
    else:
        try:
            with open(meta_path, 'r') as fh:
                obj = json.load(fh)
        except IOError:
            obj = []
        if isinstance(obj, list):
            itr = obj
        elif isinstance(obj, dict) and 'assemblies' in obj:
            itr = obj['assemblies']
        else:
            itr = [obj]

    for rec in itr:
        rec = _coerce_assembly_record(rec)
        acc = _get_accession(rec)
        org = _get_organism(rec)
        if acc and org and acc not in mapping:
            mapping[acc] = org
    return mapping

def sanitize(label):
    # Replace filesystem-unfriendly characters with underscores
    return re.sub(r'[ /()\[\],:;]', '_', label or '')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--refs-dir', required=True)
    p.add_argument('--out-dir',  required=True)
    p.add_argument('--query',    required=True)
    p.add_argument('--sample',   required=True)
    args = p.parse_args()

    refs_dir = os.path.abspath(args.refs_dir)
    out_dir  = os.path.abspath(args.out_dir)
    query    = os.path.abspath(args.query)
    sample   = args.sample

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # 1) Stage query FASTA
    qname = sanitize(sample) + '.fna'
    shutil.copy2(query, os.path.join(out_dir, qname))

    # 2) Resolve metadata file (.jsonl preferred, fallback to .json)
    meta_jsonl = os.path.join(refs_dir, 'ncbi_dataset', 'data', 'assembly_data_report.jsonl')
    meta_json  = os.path.join(refs_dir, 'ncbi_dataset', 'data', 'assembly_data_report.json')
    meta_path  = meta_jsonl if os.path.exists(meta_jsonl) else (meta_json if os.path.exists(meta_json) else None)
    mapping    = load_metadata(meta_path) if meta_path else {}

    # 3) Copy each FASTA with sanitized label "<label>__<accession>.fna"
    acc_re = re.compile(r'(GC[AF]_\d+\.\d+)')
    for root, _, files in os.walk(refs_dir):
        for fn in files:
            if fn.endswith('.fna') or fn.endswith('_genomic.fna'):
                pth = os.path.join(root, fn)
                m = acc_re.search(pth)
                if not m:
                    continue
                acc = m.group(1)
                label = sanitize(mapping.get(acc, acc))
                dst = os.path.join(out_dir, '%s__%s.fna' % (label, acc))
                shutil.copy2(pth, dst)

if __name__ == '__main__':
    sys.exit(main())

