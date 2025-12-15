Date: 12/15/2025
@AUTHOR: José Cantú 

# Goal in this Project
Generation of a publication-grade genome assembly for Strain Pseudo and determine whether it represents a novel species/genus within Oscillospiraceae. 

## Current Best Assembly 
- Assembly used for taxonomy comparisons: `assembly.polished.fasta` 
- Size: 2,617,187 bp 
- GC: 59.52%
- Proteins (TYGS): 2,593
### From Quast using Flye since best overall stats from all assemblies tried:
* contigs (>= 0 bp)         6             
* contigs (>= 1000 bp)      6             
* contigs (>= 5000 bp)      5             
* contigs (>= 10000 bp)     4             
* contigs (>= 25000 bp)     4             
* contigs (>= 50000 bp)     4             
* Total length (>= 0 bp)      2575873       
* Total length (>= 1000 bp)   2575873       
* Total length (>= 5000 bp)   2574659       
* Total length (>= 10000 bp)  2564710       
* Total length (>= 25000 bp)  2564710       
* Total length (>= 50000 bp)  2564710       
* contigs                   6             
* Largest contig              1329041       
* Total length                2575873       
* GC (%)                      59.68         
* N50                         1329041       
* N90                         1018847       
* auN                         1098850.5     
* L50                         1             
* L90                         2             
* N's per 100 kbp           0.00

## Genome Quality (CheckM) 
- 98.66% Completeness 
- 0.67% Contamination 

## 16S rRNA 
- Barrnap-extracted 16s rRNA genes: 2 copies, each 1,529 bp `16S.fna` 
- To be used as support context only (not primary species/genus factor or delimiter) 

## Novelty evidence Summary 
**GTDB-TK** places strain Psuedo in family Oscillospiraceae, near place holder genus `g__Marseille-P3106`.

**TYGS (GDBP -> dDDH)** all dDDH values to type strains are far below the ~70% species boundary. 
- Best hit here Flintibacter porci P01025 dDDH(d4)=26.4%
- Vermiculatibacterium agrestimuris DSM 112226T dDDH(d4)=20.8%
See `Type_Strain_Genome_Server_dDDH_values_comparison_with_Psuedo_assemly_polish.*` 

**EzAAI** strongest genus level similarity is to Vermiculatibacterium agrestimuris.
- AAI or Average Amino Acid Identity = 74.109% 
- Proteome Cov. = 0.5537 (POCP = 55.37%) 
- To settle genus placement I plan to do the following: use the bac120 summary from gtdb-TK and build a phylogenomic tree with the curated sampling list provided from TYGS + boot strap support if low percentage that would settle novel genus. 

## Genome Closure Status 
Not fully closed due to long repeats (~70kb & 140kb) between major contigs. With ONT reads length N50 = 3.4 kb, repeat resolution is limited w/o substantiall longer ont reads spanning repeat regions. 

## Assembly Closure Summary 

### Why Closure is difficult for this prokaryotic genome 
Bandage/graph inspection indicates two large contigs are separated by long repeat regions (~70 kb and ~140 kb). Resolution requires multiple long reads spanning each repeat plus unique flanks.

### ONT read QC (NanoPlot) 
- Reads: 610,420
- Total Bases: 1.432 Gb (~550x coverage vs ~2.6 Mb genome) 
- Read length N50: 3,389 bp 
- Longest read: 177,575 bp (next: 145,754 bp) 
My interpretation: despite the high coverage, most reads are too short to robustly span the long repeats. I will try an insiloco filtering of the longest read and 2nd one to see if they are enough to close the genome assuming unique flanks are found to bridge the long repeat regions. 

### Bridging Diagnositcs | ONT mapped to assembly 
Shared read IDs between end-windows of the two largest contigs 
- 20 kb windows 0 shared reads (all orientations) 
- 50 - 100 kb window several shared reads depending on orientation 
My interpreatation from getting these stats: weak/ambigious bridging signal here, likely dominated by repeat mapping rather than unique spanning evidence. 

### Assemblies attempted 
- Unicycler hybrid and bold: bold mode improved contiguity but remained fragmented (16 contigs from hybrid normal and 12 from bold mode). 
- Raven: The best assembly I got (5 contigs).
- Flye: (Seq Coast Run rereplicated - 6 contigs). 
- Autocycler: limited benefit because input assemblies need to be near complete ~2/3 contigs; input assemblies some were completely fragmented best one was 12 contigs. Repeat structure exceeded typical read span. 

### Practical Next Steps for closure experiments (Low cost effort on my end) 
1) In-silico filtering of ONT reads keeping only long ONT reads at moderate coverage + reassembly with Raven/Flye and recheck graph - if closed then check if assembly is accurate 
2) Evaluate whether the longest read you have can bridge the repeats uniquely (Will need to check unique anchors) 
3) If closure is truly impossible ~ likely given the dataset thus far then a solid option would be to proceed with the high quality draft of 5/6 contigs since they meet the criteria that the Journal species - so then comes a question is briding truly necessary if the requirements are already met for the Journal? A high quaity draft is better than using the last resort of trying to bridge the repeat regions with a separate organism since that bias all the reads already given and risks it being a chimeric genome. 

### Decision NOte 
Closing the genome is not stricly required given the genome already meets the criteria specified thus far. I will include the rest of the values and figures IJSEM specifies. In IJSEM own words: "While the ideal is to have a closed genome, in most circumstances this is not
possible or cost-efficient, resulting in contigs of various lengths. Very short contigs are usually excluded from the final assemblies.
However, as there is not a clear standard on how to select contigs, the absolute number of contigs is not a really good quality
indicator. Indices like the N50 (length of the shortest contig that accumulatively show 50% or more of the genome size) or L50
(smallest number of contigs that sums half of genome size) give better assessment of the genome quality."  

## Taxonomy & Novelty Evidence (Based on IJSEM Criteria) 

### Taxon Content 
Strain Psuedo is placed within 
- o__Oscillospirales; f__Oscillospiraceae 
- GTDB-Tk places it near placeholder genus `g__Marseille-P3106`.

### Species-level delimitation (primary evidence)
TYGS GBDP/dDDH comparisons to type strains show all values far below the ~70% species boundary:
- Best dDDH(d4): 26.4% vs Flintibacter porci P01025
- dDDH(d4): 20.8% vs Vermiculatibacterium agrestimuris DSM 112226T
Conclusion: no species-level match to included type strains -> candidate novel species.


### Genus-level similarity (supporting OGRI)
EzAAI comparisons against reference/type genomes:
- Best AAI: 74.109% to Vermiculatibacterium agrestimuris DSM 112226T
- “Proteome cov.” 0.5537 corresponds approximately to POCP 55.37%
Interpretation: strongest proteome-level affinity is toward Vermiculatibacterium.

### Trees included (how to interpret)
- TYGS genome GBDP phylogram (`genome_GBDP_phylogram.svg`): whole-genome distance tree; bootstrap-supported; strongest phylogenetic evidence among current trees.
- TYGS 16S GBDP phylogram (`16S_GDBP_phylogram.svg`): supportive only; 16S can lack genus-level resolution and branch lengths can be biased if sequences are incomplete (not the case here; 16S ~1529 bp).
- EzAAI dendrogram (`aai.nwk`): hierarchical clustering on AAI-derived distances; no bootstrap support; use as similarity visualization, not as a formal phylogenomic tree.
- Will include a genus phylogram here.... 

### Analysis to settle if novel genus? 
Build a marker gene phylogenomic tree using the curated list from TYGS 
- Use GTDB-TK bac120 marker gene alignment and have gtbd-TK build the marker gene tree 
- Pass that marker gene tree over to IQ-TREE 2/3 and report the support bootstrap values to see evolutionary distances to settle if novel or not. 
This will then serve as the primary genus placement anchor, with dDDH values already in place for the TYGS tree tree and AAI/POCP as supporting OGRI.

### Guide lines for support of a novel genus when looking at the created phylogenetic tree 
When claiming a novel genus, the phylogenomic tree should include type strains and use single copy genes. The isolate must form a distinct branch, with strong support (e.g., bootstrap ≥90). If the genus is not novel, the isolate should nest within an existing genus, supported by high bootstrap values.

# In‑silico length filtering + reassembly

Goal: see if dataset already contains enough very long reads to bridge the ~70 kb and ~140 kb repeats.

Make read subsets like: >30 kb, >50 kb, >80 kb, >120 kb and cap depth (e.g., keep ~60–100×)

Reassemble with Raven and Flye on the filtered sets

Check graph in Bandage + see if contigs collapse

# “Do we actually have bridging reads?” test with adjustments in mind 

My bridging-summary output already hints at this kind of analysis, but to make it repeat relevant.

Only consider long reads (e.g., ≥120 kb for the 140 kb repeat; ≥80 kb for the 70 kb repeat)

Require each candidate bridging read to have unique anchors on both sides (high MAPQ, or unique k-mers), otherwise repeats will create false “bridges” which will reduce overall accuracy of accembly and QV score from Yak. 

If I can't find even one read that anchors uniquely on both sides of the repeat, closure is impossible with current data, regardless of assembler.

# Other options for bridging the genome?
I've listed long ont as an option given the recommendations from Ryan its the standard option for long repeat regions since it costs to much this is a no go. 
Cas9 bridging being to technical in designing gRNAs guides that would act as unique anchors so sequencing enriches and span the repeat regions to bridge on both sides { A lot can go wrong here } 
reference guided briding using a close reference based closed genome - however as I mentioned this will bias the results if the selling point here of the paper is something novel ~ risks include chimerism of the genome and bias and changed all the stats thus far; hide what was novel about it to begin with {e.g. such as the repeat arrangment}; it can inflant the dDDH + AAI/POCP scoring towards the reference genome being used. 

# Why does IJSEM using Type Strains as a reference point when settling for novelty for genus/species?

IJSEM is a journal that records and validates prokaryotic names under the rules of nomenclature provided by the international code of nomeclature of prokaryotes ICNP, recommending authors to follow this format when naming, and the name of a spcies is anchored to its nomenclatural type for bacteria, that is usually a living strain which is the type strain that is perserved and accessible to the community to use. That's why comparisons to type strains or the type species' type strains when making genus level decisions if its novel or not is the community standard for "what this name actually refers to." Having a novel species/genus claim that are framed as "compared to the closest type strain" because that is metric used that fixes the name in a reproducible way. 

































































