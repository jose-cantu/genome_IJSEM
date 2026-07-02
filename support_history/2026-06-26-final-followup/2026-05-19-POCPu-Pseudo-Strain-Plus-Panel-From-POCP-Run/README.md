@AUTHOR: JC 
Date Updated: 05/24/2026 

# Project Overview 

This project runs POCPu on strain PSEUDO and the same reference accession set used in the previous POCP workflow. Reusing the same genome/accession panel keeps the comparison controlled, so POCP and POCPu results can be interpreted side-by-side.

The goal is to test whether the protein-content similarity between PSEUDO and nearby taxa remains high when conserved proteins are counted using the POCPu approach, which emphasizes unique protein matches and reduces possible inflation from duplicated or paralogous proteins.

This analysis is being used to support genus-level placement of PSEUDO by comparing it against the closest AAI/genus-tree relatives, especially Vermiculatibacterium and Lawsonibacter. POCPu results will be interpreted alongside AAI, phylogenomics, TYGS/dDDH, morphology, and previous POCP results; they are not treated as standalone taxonomic evidence.

Update: 

One thing I did notice it uses `/share/apps/anaconda3/2023.03-1/bin/nextflow` instead of the conda env I have set up so will need to update PATH to 
`export PATH="$CONDA_PREFIX/bin:$JAVA_HOME/bin:$PATH"`

so I see `/mount/britton/Jose/conda-envs/consprot-dev/bin/nextflow`.

Keep in mind for future scripts. 
