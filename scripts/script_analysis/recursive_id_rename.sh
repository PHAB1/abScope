#!/bin/bash

# Define the pattern to find an asterisk followed by a number in the first column
pattern_asterisk="^([^\t]*)\\*([0-9]+)(.*?\t)"
replacement_asterisk="\1_\2\3"

# Define the pattern to find the pipe character
pattern_pipe="\\|"
replacement_pipe="_"

# Loop through all .tsv files in the current directory
for file in *.tsv; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file" .tsv)
    echo "Processing file: $file"
    # Use sed to perform the asterisk substitution
    sed -i -E "s/$pattern_asterisk/$replacement_asterisk/g" "$file"
    # Use sed again to perform the pipe substitution
    sed -i -E "s/$pattern_pipe/$replacement_pipe/g" "$file"
  fi
done

echo "Finished processing all .tsv files."
