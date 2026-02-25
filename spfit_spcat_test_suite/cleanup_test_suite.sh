#!/bin/bash
# Clean up test suite by removing output directories, keeping only references

cd /home/briggsjt/WORKING/Pickett/spfit_spcat_test_suite

# Find all molecule directories (containing v2008_results or subdirectories)
for category_dir in */; do
    category_dir="${category_dir%/}"
    
    for mol_dir in "$category_dir"/*/; do
        mol_dir="${mol_dir%/}"
        mol_name=$(basename "$mol_dir")
        
        # Skip if not a directory
        [ ! -d "$mol_dir" ] && continue
        
        # Remove all subdirectories except v2008_results and reference_outputs
        for subdir in "$mol_dir"/*/; do
            [ ! -d "$subdir" ] && continue
            subdir_name=$(basename "$subdir")
            
            if [ "$subdir_name" != "v2008_results" ] && [ "$subdir_name" != "reference_outputs" ]; then
                echo "Removing: $subdir"
                trash-put "$subdir"
            fi
        done
    done
done

echo "Cleanup complete!"
