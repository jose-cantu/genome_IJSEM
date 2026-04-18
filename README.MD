# genome_IJSEM

Genome-taxonomy and manuscript-support repository for **Strain Pseudo**, a candidate novel taxon within **Oscillospiraceae**.

## Project goal

Generate a publication-grade genome assembly for Strain Pseudo and evaluate whether it represents a **novel species** and potentially a **novel genus** using genome-based taxonomy workflows aligned with IJSEM expectations.

## Current status

### Best assembly currently used for taxonomy comparisons

The current assembly used for downstream taxonomy comparisons is `assembly.polished.fasta`.

**Selected taxonomy assembly**
- Selected as the current best assembly for taxonomy comparisons
- Size: **2,617,187 bp**
- GC: **59.52%**
- TYGS proteins: **2,593**

**QC summary for the selected assembly**
- Contigs: **6**
- Largest contig: **1,329,041 bp**
- N50: **1,329,041**
- N90: **1,018,847**
- auN: **1,098,850.5**
- L50: **1**
- L90: **2**
- Ns per 100 kbp: **0.00**
- CheckM completeness: **98.66%**
- CheckM contamination: **0.67%**

**Important assembly provenance note**
This repository uses the `2025-09-22` workflow as the **assembly/QC/taxonomy workflow backbone**, but the final assembly selected for taxonomy comparisons is documented as the **Flye-derived** `assembly.polished.fasta` based on overall downstream use and QC. Raven produced a strong alternative assembly with fewer contigs in exploratory comparison, but the current downstream taxonomic summaries are anchored to the Flye-selected draft.

## Taxonomic interpretation summary

### Placement
- GTDB-Tk places Strain Pseudo in **Oscillospiraceae**
- GTDB-Tk places it near placeholder genus **`g__Marseille-P3106`**

### Species-level evidence
TYGS GBDP/dDDH comparisons to type strains are all far below the usual species boundary (~70%), supporting **candidate novel species** status.

Examples:
- **26.4% dDDH(d4)** vs *Flintibacter porci* P01025
- **20.8% dDDH(d4)** vs *Vermiculatibacterium agrestimuris* DSM 112226T

### Genus-level supporting evidence
EzAAI currently shows the strongest proteome-level similarity to *Vermiculatibacterium agrestimuris*:
- **AAI = 74.109%**
- **Proteome coverage = 0.5537**
- Approximate POCP context: **55.37%**

### 16S context
Barrnap-extracted 16S support files are included for context only:
- **2 copies**
- **1,529 bp each**
- supportive, but not treated as the primary delimiter for species/genus claims

## ANI provenance

An initial FastANI screen against GTDB-derived neighbor genomes recovered the same Marseille-P3106 neighborhood indicated by GTDB-Tk, with a top ANI of approximately **97.44** in that exploratory neighbor panel.

However, this ANI result is retained as **supporting neighborhood evidence**, not as the primary taxonomic anchor, for two reasons:

1. the older ANI workflow produced a raw `fastani.tsv` against a GTDB-derived neighbor set rather than a finalized curated nomenclatural panel;
2. a later stricter ANI workflow that attempted to enforce type-material neighbors failed during panel construction and is retained in the repository for provenance only.

In addition, manual inspection of the downloaded metadata for the top Marseille-P3106 RefSeq neighbor indicated that the record was flagged as **contaminated/suppressed**, which is important context for why the exploratory ANI screen is not centered in the final taxonomic argument.

As a result, final taxonomic interpretation in this repository is based primarily on:
- **TYGS/dDDH**
- **EzAAI**
- **phylogenomic reconstruction**

rather than the initial Marseille-based ANI screen alone.

## Genome closure status

## Genome closure status

The genome is **not fully closed**, but the current assembly is still a **high-quality draft suitable for taxonomy-focused analysis**.

Current evidence indicates that closure is being limited by **long repeat structure** rather than by poor overall genome quality. Based on graph inspection and downstream review, two major contigs appear to be separated by repeat regions of approximately **~70 kb** and **~140 kb**. Resolving those repeats would require multiple long ONT reads that span the repeat plus unique flanking sequence on both sides.

### Why closure is difficult in this dataset

Although the ONT dataset has high total depth, the read-length distribution is not ideal for confidently bridging these repeat regions:

- Reads: **610,420**
- Total ONT bases: **1.432 Gb** (~550× coverage relative to a ~2.6 Mb genome)
- Read length N50: **3,389 bp**
- Longest read: **177,575 bp**
- Second longest read: **145,754 bp**

In practice, this means overall depth is high, but **most reads are much shorter than the repeat regions that need to be resolved**. As a result, repeat-spanning evidence is weak and can be confounded by ambiguous repeat mapping.

Bridging diagnostics were also not strongly supportive of confident closure:
- **20 kb end-window comparison:** **0 shared reads** across orientations
- **50–100 kb windows:** some shared reads were observed, but the signal was interpreted as **weak/ambiguous** and likely influenced by repeat-associated mappings rather than unique bridging evidence

### Assembly context

Multiple assembly strategies were explored:

- **Unicycler hybrid**: fragmented (**16 contigs**)
- **Unicycler bold**: improved but still fragmented (**12 contigs**)
- **Raven**: strongest contiguity in exploratory comparison (**5 contigs**)
- **Flye**: selected as the current best assembly for taxonomy comparisons (**6 contigs**)
- **Autocycler/Trycycler-style workflow**: limited benefit in this dataset because the input assemblies were not close enough to complete and the repeat structure still exceeded typical read span

### Why a high-quality draft is acceptable here

The selected assembly used for taxonomy comparisons remains strong by the metrics most relevant to IJSEM-style genome quality assessment:

- Assembly size: **2,617,187 bp**
- GC: **59.52%**
- Contigs: **6**
- Largest contig: **1,329,041 bp**
- N50: **1,329,041**
- L50: **1**
- CheckM completeness: **98.66%**
- CheckM contamination: **0.67%**
- Barrnap-extracted 16S copies: **2**, each **1,529 bp**

So while the genome is not closed, the current evidence suggests that this is a case of the read length being too short in compareison to the repeat regions rather than an assembly that is unusable for taxonomic inference. For that reason:

- a **high-quality draft** is already sufficient for the current taxonomy work
- closure-support experiments are preserved in this repository as **support history**

## Repository layout

```text
.
|-- README.MD
|-- assembly_qc
|   |-- README_2025-09-22.md
|   |-- results
|   |   `-- 2025-09-22
|   |-- support_history
|   |   `-- 2025-12-15-insilico
|   `-- workflow_backbone
|-- docs
|   `-- PROVENANCE.md
|-- genus_identity
|   |-- aai
|   |   |-- results
|   |   `-- workflow
|   `-- phylogeny
|       |-- results
|       `-- workflow
`-- species_identity
    `-- ani
        |-- README.md
        |-- archive
        |   |-- current_failed_type_enforced
        |   `-- old_success
        `-- results
            `-- 2025-09-30
```

## Module guide

### `assembly_qc/`
Main assembly/QC workflow backbone and the selected summary artifacts used for manuscript-relevant assembly evaluation.
- `workflow_backbone/` = main workflow scripts from the IJSEM-oriented assembly/QC backbone
- `results/2025-09-22/` = selected summary artifacts (QUAST, CheckM, BUSCO, Yak, 16S)
- `support_history/2025-12-15-insilico/` = closure-support experiments and exploratory repeat-resolution work

### `species_identity/ani/`
ANI provenance and workflow history.
- `archive/old_success/` = older ANI workflow that produced the raw `fastani.tsv`
- `archive/current_failed_type_enforced/` = stricter later ANI prototype that failed during type-material panel construction
- `results/2025-09-30/` = raw ANI result, refs list, and metadata panel

### `genus_identity/aai/`
EzAAI workflow and outputs.
- `workflow/` = AAI pipeline files
- `results/` = `Pseudo_vs_refs.aai.tsv`, `aai.tsv`, `aai.nwk`

### `genus_identity/phylogeny/`
Genus-placement phylogeny workflow and outputs.
- `workflow/` = bac120 / IQ-TREE genus phylogeny workflow files
- `results/` = current FastTree / IQ-TREE tree artifacts

### `docs/`
Repository-level provenance notes.

## What is finalized vs still being curated

### Finalized enough for repo provenance
- workflow backbone organization
- selected taxonomy assembly summary
- TYGS/dDDH evidence
- EzAAI module
- ANI provenance split into old-success vs current-failed
- genus phylogeny workflow/results staging

### Still under active scientific curation
- final curated **genus tree sampling set**
- final **outgroup choice / rooting presentation**
- final manuscript wording for genus placement
- final decision on whether Raven/Flye comparison needs a dedicated summary table in the paper

## Provenance principle

This repository is organized around the **functional modules that produced the manuscript-relevant results**, not as a direct mirror of all original dated run directories. Original run provenance is preserved in `docs/PROVENANCE.md` and in archived workflow subdirectories where needed.

