#!/bin/bash

input_file="tidy-transfer-log-file.txt"
counter=0

# Check if the tidy-transfer-log-file.txt file exists
if [ ! -f "$input_file" ]; then
    echo "Error: $input_file does not exist."
    exit 1
fi

# Read the tidy-transfer-log-file.txt file and print corresponding files
while IFS= read -r line; do
    #define variables
    if [[ "$line" == "Symbolic Link:"* ]]; then
        X="${line#*: }"
    elif [[ "$line" == "Target File :"* ]]; then
        Y="${line#*: }"
        # safety checks
        if [ -n "$X" ] && [ -n "$Y" ] && [ -e "$X" ] && [ -e "$Y" ]; then
            #echo "Removing symbolic link at $X"
            rm -rf "$X"
            #echo "Moving target files from $Y to $X"
            mv "$Y" "$X"
            #echo "Creating new link from $X to $Y"
            ln -s "$X" "$Y"
            #echo "--------------------------"
            ((counter=counter+1))
        else
            echo "Error: Invalid paths or variables."
        fi
    fi
done < "$input_file"

echo "Transfer complete, $counter files moved"