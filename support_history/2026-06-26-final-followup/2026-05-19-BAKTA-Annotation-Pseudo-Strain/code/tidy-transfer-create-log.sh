#!/bin/bash

source "$(dirname "$0")/UPDATE_ME.sh"

input_directory="${main_dir}/workflow-outputs"
output_file="tidy-transfer-log-file.txt"

# Find symbolic links recursively and write to the output file
find "$input_directory" -type l -print | while read symlink; do
    target=$(readlink -f "$symlink")
    echo "Symbolic Link: $symlink" >> "$output_file"
    echo "Target File : $target" >> "$output_file"
    echo "--------------------------" >> "$output_file"
done

echo "Symbolic link information written to $output_file"
