Date: 2025/12/15
@AUTHOR: JC 

# Project 
Here I will be doing three things. 
1) I will be doing in-silico filtering to see if I can produce more long reads using `Psuedo` wgs ONT for assemblies using Raven and Flye to see if enough multiple reads can
bridge the repeat regions (140k & 70k) since more than half of the assemblies are made up of short reads. {2 best assemblies have scores of N50 being 3.4 kbp). 

2) I will be doing a patch bridge test for each of the two assemblies to require unique anchors or flanks for the repeat regions using MAPQ of 20 to see if any unique kmers are 
on either side of the repeat regions. 

3) I will be figuring out what kind of repeats those 2 repeat regions are. Typically they are mobile genetic elements such as prophages which can be tens of kb or even longer. 
Two similar prophages or one prophage duplicated can create long near-identical regions. I can use minimap2 and awk or samtools or filtering data in bakta. Fastest is using minimap2. 

