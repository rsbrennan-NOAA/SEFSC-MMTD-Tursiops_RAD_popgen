#!/usr/bin/env python3
"""
Moments demographic inference for Tursiops populations
"""

import moments
import numpy as np
import matplotlib.pylab as plt
import demes
import demesdraw
import os
import matplotlib.cm as cm
import time
from datetime import datetime
import csv
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run demographic inference for Tursiops populations.")
parser.add_argument("model_number", help="Model number to specify the model (e.g., '01', '02').")
parser.add_argument("rep_number", help="Repetition number for the model run.")  # REPNUMBER as string

args = parser.parse_args()

MODEL_NUMBER = args.model_number 
REPNUMBER = args.rep_number 

# define file paths
PERTURB = 1
MAXITER = 10
METHOD = "powell"

DEME_GRAPH_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/fourpop_{MODEL_NUMBER}.yaml"
OPTIONS_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/scripts/demography/moments/options_fourpop_{MODEL_NUMBER}_Simple.yaml"
OUTPUT_PATH = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/model_{MODEL_NUMBER}-bestfit_{METHOD}.yaml"

OUTPUT_RESULTS_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/optimization_model_{MODEL_NUMBER}_{REPNUMBER}_{METHOD}_results.txt"  # Include REPNUMBER in output file
PLOT_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/figures/moments/optimization_model_{MODEL_NUMBER}_{REPNUMBER}_{METHOD}.png"  # Include REPNUMBER in plot file name
SUMMARY_FILE = f"/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/initialRuns/model_{MODEL_NUMBER}_summary.csv"  # Summary file per model


# Print file paths for verification
print(f"Deme graph file: {DEME_GRAPH_PATH}")
print(f"Options file: {OPTIONS_PATH}")
print(f"Output plot: {OUTPUT_PATH}")

# record start time
start_time = time.time()
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"[{current_time}] Starting moments analysis for model {MODEL_NUMBER}, rep {REPNUMBER}")

# load the sfs 
print(f"Loading SFS")
fs = moments.Spectrum.from_file("/home/rbrennan/Tursiops-RAD-popgen/analysis/moments/fourpop_sfs/dadi/Coastal_Atlantic-Coastal_Gulf-Intermediate-Offshore.sfs")
    
# pop ids are in wrong order- fix them
new_pop_ids = ['CAtlantic', 'CGulf', 'Intermediate', 'Offshore']
fs = moments.Spectrum(fs, pop_ids=new_pop_ids)
print(f"Population IDs: {fs.pop_ids}")

# Run optimization
print(f"Starting demographic model optimization for model {MODEL_NUMBER}, rep {REPNUMBER}")
ret = moments.Demes.Inference.optimize(
    DEME_GRAPH_PATH,
    OPTIONS_PATH,
    fs,
    method=METHOD,
    output=OUTPUT_PATH,
    overwrite=True,
    verbose=1,
    perturb=PERTURB,
    maxiter=MAXITER
)


# Process results
param_names, opt_params, LL = ret
print(f"Log-likelihood: {-LL}")
print("Best fit parameters:")

# save to file
with open(OUTPUT_RESULTS_FILE, "w") as f:
    f.write(f"Log-likelihood: {-LL}\n")
    f.write("Best fit parameters\n")
    for n, p in zip(param_names, opt_params):
        line = f"{n}\t{p}"
        f.write(line + "\n")
        print(line)

# Create a dictionary of results to save to the CSV summary file
results_dict = {
    'Run': REPNUMBER,
    'LogLikelihood': -LL
}

# Add all parameters to the dictionary
for n, p in zip(param_names, opt_params):
    results_dict[n] = p

# Check if summary file exists and create/append appropriately
file_exists = os.path.isfile(SUMMARY_FILE)

with open(SUMMARY_FILE, 'a', newline='') as csvfile:
    fieldnames = ['Run', 'LogLikelihood'] + list(param_names)
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    # Write header only if the file doesn't exist yet
    if not file_exists:
        writer.writeheader()
    
    writer.writerow(results_dict)

print(f"Results appended to summary file: {SUMMARY_FILE}")

end_time = time.time()
end_time_print = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

elapsed_time = end_time - start_time
print(f"Analysis for model {MODEL_NUMBER}, rep {REPNUMBER} ended at {end_time_print}")
print(f"Analysis completed for model {MODEL_NUMBER}, rep {REPNUMBER} in {elapsed_time / 3600:.2f} hours")


# plot
try:
    graph = demes.load(OUTPUT_PATH)
    fig = plt.figure(figsize=(10, 8))
    demesdraw.tubes(graph, ax=fig.add_subplot(111))
    plt.savefig(PLOT_FILE, dpi=300)
    print(f"Saved plot to {PLOT_FILE}")
except Exception as e:
    print(f"Failed to generate plot for model {MODEL_NUMBER}, rep {REPNUMBER}: {e}")

