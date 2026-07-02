import os

# Define the directory where the genome assemblies are located
# I copied them from the POCP_genomes directory from my pocp work 
genome_dir = "/mount/britton/Jose/Jobs/Jason/2026-05-19-BAKTA-Annotation-Pseudo-Strain/data" 

# Define the CSV file to update
csv_file_path = "/mount/britton/Jose/Jobs/Jason/2026-05-19-BAKTA-Annotation-Pseudo-Strain/workflow/nextflow_readfile.csv"

# Open the CSV file in write mode
with open(csv_file_path, 'w') as csv_file:
    # Write the header
    csv_file.write("sample_id,assembly\n")
    
    # Traverse the output directory for assembly.fasta files
    for fasta_name in os.listdir(genome_dir):
        fasta_file = os.path.join(genome_dir, fasta_name)
        
        # Only include files ending in .fasta 
        if os.path.isfile(fasta_file) and fasta_name.endswith(".fasta"): 
           # Use the fasta filename without the extension as the sample ID 
           sample_id = os.path.splitext(fasta_name)[0] 

           # Write the line to the CSV file 
           csv_file.write(f"{sample_id},{fasta_file}\n")

print(f"CSV file updated: {csv_file_path}")
