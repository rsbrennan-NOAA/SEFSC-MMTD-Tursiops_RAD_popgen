#!/bin/bash
#
# Script to submit all stairwayPlot batch scripts
#
echo "Submitting all fourpop stairwayPlot runs"

# Ccounter
count=0

# Loop through all stairwayPlot_fourpop_*.sh files
for script in stairwayPlot/stairwayPlot_fourpop_*.sh; do
    if [ -f "$script" ]; then
        echo "Submitting $script"
        # Submit the script to Slurm
        sbatch "$script"
        
        # Check if submission was successful
        if [ $? -eq 0 ]; then
            echo "Successfully submitted $script"
            ((count++))
        else
            echo "Failed to submit $script"
        fi
        sleep 1
    fi
done

echo "Submission complete. Submitted $count stairwayPlot jobs."

